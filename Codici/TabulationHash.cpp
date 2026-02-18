#ifndef TABULATION_HASH_H
#define TABULATION_HASH_H

#include <cstdint>
#include <random>
typedef std::mt19937 RNG; // the Mersenne Twister with a popular choice of parameters

using namespace std;

class TabulationHash
{
private:
    uint32_t table[8][16];
    uint32_t U;

public:
    TabulationHash(uint32_t U) : U(U)
    {
        std::random_device rd; // Obtain a random seed from the hardware
        RNG rng(rd());
        std::uniform_int_distribution<uint32_t> uint_dist(0, U); // by default range [0, MAX]
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 16; j++)
                table[i][j] = uint_dist(rng);
    }

    TabulationHash()
    {
        new (this) TabulationHash(UINT32_MAX);
    }

    ~TabulationHash() {}

    uint32_t operator()(uint32_t x)
    {
        uint32_t res = 0;
        for (int i = 0; i < 8; i++)
            res ^= table[i][(uint8_t)((x >> 4 * i) & 0b1111)]; // added the end with binary 00001111 to slip the 4 least significant bits
        return res;
    }
};

#endif