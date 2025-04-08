#pragma once

#include "coproto/Socket/Socket.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Crypto/AES.h"
#include <cstdint>
#include <stddef.h>
#include <vector>
#include <array>

namespace sparse_comp::sp_linf {

    template<size_t tr, size_t t, size_t d, uint8_t delta, uint8_t ssp>
    class Sender {

        osuCrypto::PRNG* prng;
        osuCrypto::AES* aes;
            
        public:
            Sender(osuCrypto::PRNG& prng, osuCrypto::AES& aes) {
                this->prng = &prng;
                this->aes = &aes;
            }

            coproto::task<void> send(coproto::Socket& sock, std::vector<block>& ordIndexHashSet, std::array<std::array<uint32_t,d>,t>& in_values, std::array<std::array<block,1>,t>& out_vec_shares);
    };

    template<size_t ts, size_t t, size_t d, uint8_t delta, uint8_t ssp>
    class Receiver {

        osuCrypto::PRNG* prng;
        osuCrypto::AES* aes;
            
        public:
            Receiver(osuCrypto::PRNG& prng, osuCrypto::AES& aes) {
                this->prng = &prng;
                this->aes = &aes;
            }
            
            coproto::task<void> receive(coproto::Socket& sock, std::vector<osuCrypto::block>& ordIndexHashSet, std::array<std::array<uint32_t,d>,t>& in_values, std::array<std::array<block,1>,t>& z_vec_shares);
    };

}

#include "./SpLInf.cpp"
