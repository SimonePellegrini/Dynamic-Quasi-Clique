#ifndef DYNAMIC_MINHASH_H
#define DYNAMIC_MINHASH_H

#include <cstdint>
#include "TabulationHash.cpp"

class DynamicMinHash
{
private:
    /**
     * Number of buffers (i.e., size of the MinHash signature)
     */
    size_t k;

    /**
     * Maximum size of each buffer
     */
    size_t buffer_max_size;

    /**
     * Buffers to store the hashed values
     */
    uint32_t **buffers;

    /**
     * Current sizes of each buffer
     */
    size_t *buffer_sizes;

    /**
     * Current thresholds for each buffer
     */
    uint32_t *thresholds;

    /**
     * Hash functions for each buffer
     */
    TabulationHash **hashes;

    /**
     * MinHash signature.
     */
    uint32_t *signature;

public:

    /**
     * Constructor
     * @param k Number of buffers (i.e., size of the MinHash signature)
     * @param buffer_max_size Maximum size of each buffer
     * @param hashes Hash functions for each buffer. They must be provided externally to ensure that the same hash functions are used across different instances of DynamicMinHash.
     */
    DynamicMinHash(size_t k, size_t buffer_max_size, TabulationHash **hashes) : k(k), buffer_max_size(buffer_max_size), hashes(hashes)
    {
        buffers = new uint32_t *[k];
        buffer_sizes = new size_t[k];
        thresholds = new uint32_t[k];
        signature = new uint32_t[k];

        for (size_t i = 0; i < k; i++)
        {
            buffers[i] = new uint32_t[buffer_max_size];
            buffer_sizes[i] = 0;
            thresholds[i] = UINT32_MAX;
            signature[i] = UINT32_MAX;

            for (size_t j = 0; j < buffer_max_size; j++)
                buffers[i][j] = UINT32_MAX;
        }
    }

    ~DynamicMinHash()
    {
        for (size_t i = 0; i < k; i++)
        {
            delete[] buffers[i];
            delete hashes[i];
        }
        delete[] buffers;
        delete[] buffer_sizes;
        delete[] thresholds;
        delete[] hashes;
        delete[] signature;
    }

    static TabulationHash **createHashFunctions(size_t k)
    {
        TabulationHash **hashes = new TabulationHash *[k];
        for (size_t i = 0; i < k; i++)
        {
            hashes[i] = new TabulationHash();
        }
        return hashes;
    }

    /**
     * Hash the input value u using the i-th hash function.
     */
    uint32_t hash(size_t i, uint32_t u)
    {
        return (*hashes[i])(u);
    }

    /**
     * Insert the element u into the Dynamic MinHash structure.
     */
    void insert(uint32_t u)
    {
        for (size_t i = 0; i < k; i++)
        {
            uint32_t h = hash(i, u);
            if (h <= thresholds[i])
            {
                if (buffer_sizes[i] < buffer_max_size)
                {
                    buffers[i][buffer_sizes[i]++] = h;
                }
                else
                {
                    // Replace the maximum value in the buffer
                    for (size_t j = 0; j < buffer_max_size; j++)
                    {
                        if (buffers[i][j] == thresholds[i])
                        {
                            buffers[i][j] = h;
                            break;
                        }
                    }
                }

                // Update threshold
                if (buffer_sizes[i] == buffer_max_size)
                {
                    uint32_t max_val = 0;
                    for (size_t j = 0; j < buffer_max_size; j++)
                    {
                        if (buffers[i][j] > max_val)
                        {
                            max_val = buffers[i][j];
                        }
                    }
                    thresholds[i] = max_val;
                }

                // Update signature
                if (h < signature[i])
                {
                    signature[i] = h;
                }
            }
        }
    }

    /**
     * Remove the element u from the Dynamic MinHash structure.
     * If one of the buffers becomes empty after removal, the entire data structure is reset, and true is returned.
     * Otherwise, false is returned.
     */
    bool remove(uint32_t u)
    {
        for (size_t i = 0; i < k; i++)
        {
            uint32_t h = hash(i, u);
            if (h > this->thresholds[i])
                continue;

            // find the element to remove
            int index_to_remove = -1;
            for (int j = 0; j < this->buffer_sizes[i]; j++)
            {
                if (this->buffers[i][j] == h)
                {
                    index_to_remove = j;
                    break;
                }
            }

            // element not found, continue to the next buffer
            if (index_to_remove == -1)
                continue;

            // remove the element
            this->buffers[i][index_to_remove] = this->buffers[i][this->buffer_sizes[i] - 1];
            this->buffer_sizes[i]--;

            // if the buffer is empty reset all the data structure
            if (this->buffer_sizes[i] == 0)
            {
                this->resetBuffer();
                return true;
            }

            // update the minimum value if deleted
            if (this->signature[i] == h)
            {
                uint32_t min_val = UINT32_MAX;
                for (int j = 0; j < this->buffer_sizes[i]; j++)
                {
                    if (this->buffers[i][j] < min_val)
                    {
                        min_val = this->buffers[i][j];
                    }
                }
                this->signature[i] = min_val;
            }
        }

        return false;
    }

    void resetBuffer()
    {
        for (size_t i = 0; i < k; i++)
        {
            buffer_sizes[i] = 0;
            thresholds[i] = UINT32_MAX;
            signature[i] = UINT32_MAX;

            for (size_t j = 0; j < buffer_max_size; j++)
                buffers[i][j] = UINT32_MAX;
        }
    }

    virtual uint32_t *getSignature()
    {
        return signature;
    };

    static double similarity(DynamicMinHash *A, DynamicMinHash *B)
    {
        uint32_t *sigA = A->getSignature();
        uint32_t *sigB = B->getSignature();

        int k = A->k;
        double c = .0;
        for (int i = 0; i < k; i++)
            c += sigA[i] == sigB[i];
        return c / static_cast<double>(k);
    }
};

#endif