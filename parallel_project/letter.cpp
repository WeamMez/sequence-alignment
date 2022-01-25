#include "letter.hpp"

Letter::Letter(char letter)
{
    init(letter);
}

//Letter::Letter(){}

void Letter::init(char letter)
{
    this->letter = letter;

/*
    if (letter == MUTANT_CHAR)
    {
        for (uint8_t i = 0; i < types_of_groups; i++)
            in_groups[i] = 0;
        return;
    }
*/

    for (uint8_t i = 0; i < types_of_groups; i++)
    { // Loop over all types of ALPHABET_GROUPS
        in_groups[i] = 0;
        for (uint8_t j = 0; j < num_of_groups[i]; j++)
        { // Loop over all ALPHABET_GROUPS in a specific type
            in_groups[i] |= ALPHABET_GROUPS[i][j].find(letter) == string::npos ? 0 : (1 << j);
            // if letter in group add 1 in place
        }
    }
}

ostream& operator<<(ostream& os, const Letter& l)
{
    os << l.letter;
    return os;
}