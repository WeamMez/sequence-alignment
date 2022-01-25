
#ifndef _LETTER_HPP
#define _LETTER_HPP

#include <iostream>
#include <array>
#include <string>
#include <vector>
#include <limits>
#include <mpich/mpi.h>

//#define DEBUG
//#define LONG_DEBUG
//#define SHOW_WEIGHTS
//#define DEBUG_NO_MUTANT
//#define SHOW_TIME

#define MUTANT_CHAR '-'
#define CROSS_RESULT_COUNT 5
#define CROSS_RESULT_READ_FROM 2
enum CROSS_RESULT {Not_to_count, Identical, Type1, Type2, None};

using namespace std;

#define CROSS_RESULT_SPECIALS 2
constexpr uint8_t num_of_groups[] = {9, 11};
constexpr uint8_t types_of_groups = sizeof(num_of_groups) / sizeof(num_of_groups[0]);

const array<array<string, num_of_groups[1]>, types_of_groups> ALPHABET_GROUPS =
{{{"NDEQ", "NEQK", "STA", "MILV", "QHRK",
"NHQK", "FYW", "HY", "MILF"},
{"SAG", "ATV", "CSA", "SGND", "STPA", "STNK",
"NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"}}};


class Letter
{
    public:
        char letter;
        uint32_t in_groups[types_of_groups];

        Letter(char letter);
        //Letter();
        void init(char letter);

        friend ostream& operator<<(ostream& os, const Letter &l);
};

#endif