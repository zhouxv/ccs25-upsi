#pragma once

#include "volePSI/RsOprf.h"
#include "coproto/Socket/Socket.h"
#include "cryptoTools/Crypto/AES.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/BitVector.h"
#include <cstdint>
#include <vector>
#include <array>

using PRNG = osuCrypto::PRNG;
using AES = osuCrypto::AES;
using block = osuCrypto::block;
using BitVector = osuCrypto::BitVector;

using Proto = coproto::task<void>;

namespace sparse_comp::multi_oprf {

    const size_t ell = 128;
    const size_t comp_sec_param = 128;

    class Sender {
        private:
            std::vector<block>* randSetupOtMsgs = nullptr;
            block s;
            size_t query_num = 0;
            std::vector<block>* okvs = new std::vector<block>();
            AES aes = AES(block(13133210048402866,17132091720387928));

        public:
            ~Sender();

            static Proto setup(coproto::Socket& sock, PRNG& prng, size_t num_instances, std::vector<Sender*>& senders);

            Proto send(coproto::Socket& sock, size_t query_num);
            void eval(std::vector<block>& idxs, std::vector<block>& vals);
    };

    class Receiver {
        private:
            std::vector<std::array<block, 2>>* randSetupOtMsgs = nullptr;
            AES aes = AES(block(13133210048402866,17132091720387928));
       
        public:
            ~Receiver();

            static Proto setup(coproto::Socket& sock, PRNG& prng, size_t num_instances, std::vector<Receiver*>& receivers);

            Proto receive(coproto::Socket& sock, std::vector<block>& idxs, std::vector<block>& vals);

    };

    

};