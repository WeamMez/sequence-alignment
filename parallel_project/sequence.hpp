
#ifndef _SEQUENCE_HPP
#define _SEQUENCE_HPP

#include "letter.hpp"
#include <algorithm>
#include <vector>

#ifdef CROSS_MUTANT_RESULT
#define SIZE_ERROR {-2, 0};
#else
#define SIZE_ERROR -2
#endif

using namespace std;

constexpr uint8_t LETTER_COUNT = 'Z' - 'A' + 1;

struct Alphabet
{
    Letter letters[LETTER_COUNT];
    Letter mutant;
};

struct Fit_result
{
    int offset;
    int mutation;
    int weight;
};

#ifdef CROSS_MUTANT_RESULT
struct Cross_mutant_result
{
    int maxi;
    int maxv;
};
#else
typedef int Cross_mutant_result;
#endif 

class Sequence
{
    private:
        vector<const Letter*> sequence;
        const Alphabet *alphabet;
    
    public:
        Sequence(const Alphabet *alphabet, vector<const Letter*> sequence);
        Sequence(const Alphabet *alphabet, string *str);
        Sequence(const Alphabet *alphabet);
        Sequence(const Sequence& other);
        ~Sequence();
        int size() const;
        Sequence& operator+=(const string *str);
        Sequence& operator+=(const Letter *letter);
        Sequence& operator+=(char ch);
        const Letter *operator[](int index) const;
        const vector<const Letter *> get_sequence() const;

        friend Letter **seq_to_cuda(const Alphabet *cuda_alphabet, const Sequence *seq);
        friend Fit_result seq_fit(const Sequence *main_seq, const Sequence *other_seq, array<int, CROSS_RESULT_COUNT> *weights);

        friend ostream& operator<<(ostream& out, const Sequence& seq);

};

#endif // _SEQUENCE_HPP