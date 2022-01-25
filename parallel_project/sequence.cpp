#include "sequence.hpp"

Sequence::Sequence(const Alphabet *alphabet, vector<const Letter*> sequence) : Sequence(alphabet)
{
    this->sequence = sequence;
}

Sequence::Sequence(const Alphabet *alphabet, string *str) : Sequence(alphabet)
{
    *this += str;
}

Sequence::Sequence(const Alphabet *alphabet)
{
    this->alphabet = alphabet;
}

Sequence::Sequence(const Sequence &other) : Sequence(other.alphabet, other.sequence) {}

Sequence::~Sequence()
{}

int Sequence::size() const
{
    return sequence.size();
}

Sequence &Sequence::operator+=(const string *str)
{
    for (int i = 0; i < (int)str->length(); ++i)
        *this += (*str)[i];
    return *this;
}

Sequence &Sequence::operator+=(const Letter *letter)
{
    this->sequence.push_back(letter);
    return *this;
}

Sequence &Sequence::operator+=(char ch)
{
    if (ch >= 'a' && ch <= 'z')
        *this += alphabet->letters + ch - 'a';
    else if (ch >= 'A' && ch <= 'Z')
        *this += alphabet->letters + ch - 'A';
    else if (MUTANT_CHAR == ch)
        *this += &(alphabet->mutant);
    return *this;
}

const Letter *Sequence::operator[](int index) const
{
    return (this->sequence[index]);
}

const vector<const Letter *> Sequence::get_sequence() const
{
    return this->sequence;
}

ostream& operator<<(ostream& out, const Sequence& seq)
{
    for (int i = 0; i < (int)seq.sequence.size(); ++i)
        out << *(seq.sequence[i]);
    return out;
}