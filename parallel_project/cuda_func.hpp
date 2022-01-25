#ifndef _CUDA_FUNC_HPP
#define _CUDA_FUNC_HPP

#include "sequence.hpp"

#define WARP_SIZE 32
#define THREADS_PER_BLOCK 1024
#define FULL_32_BIT_MASK 0xFFFFFFFF
#define DIV_ROUND_UP(x,y) (((x) + ((y) - 1)) / (y))

void seq_cross(const Letter **main_seq, int main_size, const Letter **other_seq, int other_size, int *weights, int offset, int mutation, int *cuda_reduced_ans_temp);

Fit_result seq_fit(const Sequence *main_seq, const Sequence *other_seq, array<int, CROSS_RESULT_COUNT> *weights);

Letter **seq_to_cuda(Alphabet *cuda_alphabet, Sequence *seq);

#endif