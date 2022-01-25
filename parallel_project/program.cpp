#include "cuda_func.hpp"
#include <bits/stdc++.h>

#ifdef DEBUG
#include <cstdio>
#endif

#ifdef SHOW_TIME
//#include <ctime>
#include <chrono>
#endif

#define ROOT 0

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    int np, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

#ifdef SHOW_TIME
    int64_t start_time;
#endif

    const Alphabet alphabet = {{'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'}, '-'};

    int main_size;
    string main_seq_str;
    Sequence *main_sequence;
    vector<Sequence*> sequences;
    array<int, CROSS_RESULT_COUNT> weights;
    vector<string> sequences_str;
    vector<int> seq_partitions_count(np);
    vector<int> seq_partitions_size(np);
    vector<int> seq_partitions_index(np);
    vector<int> seq_partitions_gather_index(np);
    string all_sequences_in_one_string;
    string my_sequences_in_one_string;
    int my_long_string_length;
    vector<Sequence*> my_sequences;
    vector<Fit_result> fitres;
    vector<Fit_result> fitres_all;
    int number_of_sequences;

#ifdef SHOW_TIME
    start_time = std::chrono::system_clock::now().time_since_epoch().count();
#endif

    if (rank == ROOT)
    {

        
        weights[Not_to_count] = 0;
        cin >> weights[Identical];
        for (int i = CROSS_RESULT_READ_FROM; i < CROSS_RESULT_COUNT; i++)
        {
            cin >> weights[i];
            weights[i] = -weights[i];
        }
        cin.ignore(3000, '\n');

        getline(cin, main_seq_str);
        main_size = main_seq_str.size();
        
        cin >> number_of_sequences;
        cin.ignore(3000, '\n');
        cin.sync();
        sequences_str.resize(number_of_sequences);

        // To scatter the sequences we concatinate them into one big std::string and
        // use scatterv, then we split them by the space that we put between them

        int before_mod = DIV_ROUND_UP(number_of_sequences, np);
        int after_mod = number_of_sequences / np;
        int mod = number_of_sequences % np;

        int i, counter = 0;

        for (i = 0; i < mod; i++)
        {
            seq_partitions_count[i] = before_mod;
            seq_partitions_gather_index[i] = counter;
            counter += before_mod;
        }
        for (; i < np; i++)
        {
            seq_partitions_count[i] = after_mod;
            seq_partitions_gather_index[i] = counter;
            counter += after_mod;
        }

        fill(seq_partitions_size.begin(), seq_partitions_size.end(), 0);
        seq_partitions_index[0] = 0;

        int local_counter = 0, process_counter = 0;

        for (i = 0; i < number_of_sequences;)
        {
            string temp_seq_str;
            cin >> temp_seq_str;
            cin.ignore(3000, '\n');
            cin.sync();

            if (temp_seq_str.length() != 0)
            {
                if (local_counter == (int)seq_partitions_count[process_counter])
                {
                    process_counter++;
                    local_counter = 0;
                    seq_partitions_index[process_counter] = all_sequences_in_one_string.length();
                }

                seq_partitions_size[process_counter] += temp_seq_str.length() + 1;
                all_sequences_in_one_string.append(temp_seq_str);
                all_sequences_in_one_string.append(" ");

                i++;
                local_counter++;
            }
        }

        fitres_all.resize(number_of_sequences);
    }

    MPI_Bcast(weights.data(), weights.size(), MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&main_size, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    if (rank != ROOT)
    {
        main_seq_str.resize(main_size);
    }

    MPI_Bcast(main_seq_str.data(), main_size + 1, MPI_CHAR, ROOT, MPI_COMM_WORLD);
    main_sequence = new Sequence(&alphabet, &main_seq_str);

    int my_num_of_sequences;

    MPI_Scatter(seq_partitions_count.data(), 1, MPI_INT, &my_num_of_sequences, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    vector<string> my_sequences_str(my_num_of_sequences);
    fitres.resize(my_num_of_sequences);

    MPI_Scatter(seq_partitions_size.data(), 1, MPI_INT, &my_long_string_length, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    my_sequences_in_one_string.resize(my_long_string_length);

    MPI_Scatterv(all_sequences_in_one_string.data(), seq_partitions_size.data(), seq_partitions_index.data(), MPI_CHAR, my_sequences_in_one_string.data(), my_long_string_length, MPI_CHAR, ROOT, MPI_COMM_WORLD);

    my_sequences.resize(my_num_of_sequences);
    stringstream splitter(my_sequences_in_one_string);
    for (int i = 0; i < (int)my_num_of_sequences; i++)
    {
        splitter >> my_sequences_str[i];
        my_sequences[i] = new Sequence(&alphabet, &(my_sequences_str[i]));
    }

#ifdef DEBUG
    printf("Process %d has %zd sequences\n", rank, my_sequences.size());
#endif

    for (int i = 0; i < (int)my_sequences.size(); i++)
    {
        fitres[i] = seq_fit(main_sequence, my_sequences[i], &weights);
        delete my_sequences[i];

#ifdef DEBUG
        printf("Process %d calculated the fitting results no.%d.\n%s\nn = %d\tk = %d\tscore = %d\n", rank+1, i, my_sequences_str[i].c_str(), fitres[i].offset, fitres[i].mutation, fitres[i].weight);
#endif //DEBUG

    }

    delete main_sequence;

    Fit_result dummy_fitres;
    MPI::Datatype MPI_Fitres_type;
    MPI::Aint base_address = MPI::Get_address(&dummy_fitres);
    int mpi_fitres_block_length[] = {1,1,1};
    MPI::Aint mpi_fitres_displacements[] = {MPI::Get_address(&dummy_fitres.offset) - base_address, MPI::Get_address(&dummy_fitres.mutation) - base_address, MPI::Get_address(&dummy_fitres.weight) - base_address};
    MPI::Datatype mpi_fitres_types[] = {MPI_INT, MPI_INT, MPI_INT};
    MPI_Fitres_type = MPI::Datatype::Create_struct(3, mpi_fitres_block_length, mpi_fitres_displacements, mpi_fitres_types);
    MPI_Fitres_type.Commit();

#if defined(DEBUG) && defined(SHOW_TIME)
    printf("time until Gatherv = %lf seconds\n", ((double)(std::chrono::system_clock::now().time_since_epoch().count() - start_time)) / 1000000000.0);
#endif
    MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
    printf("disp = {%ld, %ld}\n", mpi_fitres_displacements[0], mpi_fitres_displacements[1]);

#endif

    int gather_stat = MPI_Gatherv(fitres.data(), my_num_of_sequences, MPI_Fitres_type, fitres_all.data(), seq_partitions_count.data(), seq_partitions_gather_index.data(), MPI_Fitres_type, ROOT, MPI_COMM_WORLD);

#ifdef DEBUG
    printf("gather_stat = %d\n", gather_stat);
#endif

    if (rank == ROOT)
    {
        for (int i = 0; i < number_of_sequences; i++)
        {
            cout << "n = " << fitres_all[i].offset << "\tk = " << fitres_all[i].mutation;
#ifdef SHOW_WEIGHTS
            cout << "\tscore = " << fitres_all[i].weight;
#endif
            cout << endl;
        }

#ifdef SHOW_TIME
        cout << "time = " << ((double)(std::chrono::system_clock::now().time_since_epoch().count() - start_time)) / 1000000000.0 << " seconds" << endl;
#endif

    }
    MPI_Finalize();
    return EXIT_SUCCESS;
}