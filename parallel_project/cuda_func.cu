#include "cuda_func.hpp"
#include <cuda.h>
#include <cuda_runtime.h>

__device__ CROSS_RESULT cuda_cross(const Letter *main_letter, const Letter *other_letter)
{
    if (main_letter->letter == MUTANT_CHAR || other_letter->letter == MUTANT_CHAR)
        return Not_to_count;
    if (main_letter->letter == other_letter->letter)
        return Identical;
    int8_t i;
    for (i = CROSS_RESULT_SPECIALS; i < types_of_groups + CROSS_RESULT_SPECIALS; i++)
        if (main_letter->in_groups[i - CROSS_RESULT_SPECIALS] & (other_letter->in_groups[i - CROSS_RESULT_SPECIALS]))
            return CROSS_RESULT(i);
    return None; // according to CROSS_RESULT
}

__global__ void cross_kernel(const Letter **main_seq, int main_size, const Letter **other_seq, int other_size, int *weights,
    int offset, int mutation, int *ans)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= other_size)
        return;
    int my_weight;
    if (i >= mutation)
        my_weight = (int)weights[cuda_cross((main_seq)[i + offset + 1], (other_seq)[i])];
    else
        my_weight = (int)weights[cuda_cross((main_seq)[i + offset], (other_seq)[i])];

    atomicAdd(ans, my_weight);
}

void seq_cross(const Letter **main_seq, int main_size, const Letter **other_seq, int other_size, int *weights, int offset, int mutation, int *cuda_reduced_ans_temp)
{
    int compare_until = other_size;

    unsigned int blocks = DIV_ROUND_UP(compare_until, THREADS_PER_BLOCK);
    cudaMemset(cuda_reduced_ans_temp, 0, sizeof(int));

    cross_kernel<<<blocks, THREADS_PER_BLOCK>>>(main_seq, main_size, other_seq, other_size, weights, offset, mutation, cuda_reduced_ans_temp);

    cudaDeviceSynchronize();
}

__global__ void seq_to_cuda_kernel(const Letter *const *letter, int size, const Letter *cpu_alphabet_first_letter, const Letter *cuda_alphabet_first_letter, Letter **target)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size)
        target[i] = (Letter *)(cuda_alphabet_first_letter + (letter[i] - cpu_alphabet_first_letter)); // Calculate the address of the letter in the cuda_alphabet
}

Letter **seq_to_cuda(const Alphabet *cuda_alphabet, const Sequence *seq)
{
    Letter **cuda_seq, **cuda_host_sequence;
    cudaMalloc(&cuda_seq, sizeof(Letter *) * seq->size());
    cudaMalloc(&cuda_host_sequence, sizeof(Letter *) * seq->size());
    cudaMemcpy(cuda_host_sequence, seq->sequence.data(), sizeof(Letter *) * seq->size(), cudaMemcpyHostToDevice);

    unsigned int blocks = DIV_ROUND_UP(seq->size(), THREADS_PER_BLOCK);
    seq_to_cuda_kernel<<<blocks, THREADS_PER_BLOCK>>>(cuda_host_sequence, seq->size(), seq->alphabet->letters, cuda_alphabet->letters, cuda_seq);
    cudaDeviceSynchronize();

    return cuda_seq;
}

Fit_result seq_fit(const Sequence *main_seq, const Sequence *other_seq, array<int, CROSS_RESULT_COUNT> *weights)
{

    int max_offset = 0, max_mutation = 0;
    int max_score = numeric_limits<int>::min(), score;


    Alphabet *cuda_alphabet;
    cudaMalloc(&cuda_alphabet, sizeof(Alphabet));
    cudaMemcpy(cuda_alphabet, main_seq->alphabet, sizeof(Alphabet), cudaMemcpyHostToDevice);

    Letter **cuda_main_seq = seq_to_cuda(cuda_alphabet, main_seq);
    Letter **cuda_other_seq = seq_to_cuda(cuda_alphabet, other_seq);

    int *cuda_weights;
    cudaMalloc(&cuda_weights, sizeof(int) * CROSS_RESULT_COUNT);
    cudaMemcpy(cuda_weights, weights, sizeof(int) * CROSS_RESULT_COUNT, cudaMemcpyHostToDevice);

    int *cuda_reduced_ans_temp;
    cudaMalloc(&cuda_reduced_ans_temp, sizeof(int));

    #pragma omp parallel for
    for (int offset = 0; offset < (int)(main_seq->size()) - other_seq->size(); offset++)
        for (int mutation = 1; mutation <= (int)(other_seq->size()); mutation++)
        {
            seq_cross((const Letter **)cuda_main_seq, main_seq->size(), (const Letter **)cuda_other_seq, other_seq->size(), cuda_weights, offset, mutation, cuda_reduced_ans_temp);
            cudaMemcpy(&score, cuda_reduced_ans_temp, sizeof(int), cudaMemcpyDeviceToHost);
            
#ifdef LONG_DEBUG
            printf("\toffset = %zd\tmutation = %zd\tscore = %d\n", offset, mutation, score);
#endif
            #pragma omp critical
            if (score > max_score)
            {
                max_score = score;
                max_offset = offset;
                max_mutation = mutation;
            }
        }
    seq_cross((const Letter **)cuda_main_seq, main_seq->size(), (const Letter **)cuda_other_seq, other_seq->size(), cuda_weights, main_seq->size() - other_seq->size(), other_seq->size(), cuda_reduced_ans_temp);
    cudaMemcpy(&score, cuda_reduced_ans_temp, sizeof(int), cudaMemcpyDeviceToHost);
    if (score > max_score)
    {
        max_score = score;
        max_offset = main_seq->size() - other_seq->size();
        max_mutation = other_seq->size();
    }

    cudaFree(cuda_main_seq);
    cudaFree(cuda_other_seq);
    cudaFree(cuda_weights);
    cudaFree(cuda_alphabet);

    return {max_offset, max_mutation, max_score};
}