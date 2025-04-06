#pragma once

#include "coproto/Socket/Socket.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Crypto/AES.h"
#include <cstdint>
#include <stddef.h>
#include <vector>
#include <array>

using Proto = coproto::task<void>;

template<typename T>
using vector = std::vector<T>;

template<typename T, size_t N>
using array = std::array<T, N>; 

namespace sparse_comp::sp_l1 {

    template<size_t tr, size_t ts, size_t d, uint8_t delta, uint8_t ssp>
    class Sender {

        osuCrypto::PRNG* prng;
        osuCrypto::AES* aes;
            
        public:
            Sender(osuCrypto::PRNG& prng, osuCrypto::AES& aes) {
                this->prng = &prng;
                this->aes = &aes;
            }

            Proto send(coproto::Socket& sock, array<point,ts>& ordIndexSet, array<array<uint32_t,d>,ts>& in_values, array<array<block,1>,ts>& z_vec_shares);
    };

    template<size_t ts, size_t tr, size_t d, uint8_t delta, uint8_t ssp>
    class Receiver {

        osuCrypto::PRNG* prng;
        osuCrypto::AES* aes;
            
        public:
            Receiver(osuCrypto::PRNG& prng, osuCrypto::AES& aes) {
                this->prng = &prng;
                this->aes = &aes;
            }
            
            Proto receive(coproto::Socket& sock, array<point,tr>& ordIndexSet, array<array<uint32_t,d>,tr>& in_values, array<array<block,1>,tr>& z_vec_shares);
    };

}

#include "./SpL1.cpp"
