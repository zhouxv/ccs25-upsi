#pragma once

#include "coproto/Socket/Socket.h"
#include <cstdint>
#include <stddef.h>
#include <vector>
#include <array>

using Proto = coproto::task<void>;

template<typename T>
using vector = std::vector<T>;

template<typename T, size_t N>
using array = std::array<T, N>; 

namespace sparse_comp::sp_l2 {

    template<size_t tr, size_t ts, size_t d, uint8_t delta, uint8_t ssp>
    class Sender {

        PRNG* prng;
        AES* aes;
            
        public:
            Sender(PRNG& prng, AES& aes) {
                this->prng = &prng;
                this->aes = &aes;
            }

            Proto send(Socket& sock, array<point,ts>& ordIndexSet, array<array<uint32_t,d>,ts>& in_values);
    };

    template<size_t ts, size_t tr, size_t d, uint8_t delta, uint8_t ssp>
    class Receiver {

        PRNG* prng;
        AES* aes;
            
        public:
            Receiver(PRNG& prng, AES& aes) {
                this->prng = &prng;
                this->aes = &aes;
            }
            
            Proto receive(Socket& sock, array<point,tr>& ordIndexSet, array<array<uint32_t,d>,tr>& in_values, vector<size_t>& intersec_pos);
    };

}

#include "./SpL2.cpp"
