#pragma once

#include "coproto/Socket/Socket.h"
#include <cstdint>
#include <stddef.h>
#include <vector>
#include <array>

namespace sparse_comp::sp_l2 {

    template<size_t tr, size_t ts, uint8_t delta, uint8_t ssp>
    class Sender {

        osuCrypto::PRNG* prng;
        osuCrypto::AES* aes;
            
        public:
            Sender(osuCrypto::PRNG& prng, osuCrypto::AES& aes) {
                this->prng = &prng;
                this->aes = &aes;
            }

            coproto::task<void> send(coproto::Socket& sock, std::array<point,ts>& ordIndexSet, std::array<std::array<uint32_t,2>,ts>& in_values);
    };

    template<size_t ts, size_t tr, uint8_t delta, uint8_t ssp>
    class Receiver {

        osuCrypto::PRNG* prng;
        osuCrypto::AES* aes;
            
        public:
            Receiver(osuCrypto::PRNG& prng, osuCrypto::AES& aes) {
                this->prng = &prng;
                this->aes = &aes;
            }
            
            coproto::task<void> receive(coproto::Socket& sock, std::array<point,tr>& ordIndexSet, std::array<std::array<uint32_t,2>,tr>& in_values, std::vector<size_t>& intersec_pos);
    };

}

#include "./SpL2.cpp"
