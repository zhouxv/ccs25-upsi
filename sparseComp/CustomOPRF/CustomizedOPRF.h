#pragma once

#include "../Common/Common.h"
#include "../Common/VecMatrix.h"
#include "../Common/ZN.h"
#include "../MultiOPRF/MultiOPRF.h"
#include "coproto/Socket/Socket.h"
#include "cryptoTools/Crypto/AES.h"
#include "cryptoTools/Crypto/PRNG.h"
#include <cstdint>


using PRNG = osuCrypto::PRNG;
using AES = osuCrypto::AES;

using MultiOprfRecvr = sparse_comp::multi_oprf::Receiver;
using MultiOprfSender = sparse_comp::multi_oprf::Sender;

using Proto = coproto::task<void>;

namespace sparse_comp::custom_oprf {

    struct oprf_point {
        sparse_comp::point point;
        uint32_t sot_idx;
        uint16_t sot_choice_share;

        oprf_point() = default;        

        oprf_point(sparse_comp::point point, uint32_t sot_idx, uint16_t sot_choice_share) {
            this->point = point;
            this->sot_idx = sot_idx;
            this->sot_choice_share = sot_choice_share;
        }
    };

    class Sender {

        private:
            inline static AES aes = AES(block(7,7));
            MultiOprfSender* oprfSender;

            Sender(MultiOprfSender* oprfSender);
            
        public:
            static Proto setup(coproto::Socket& sock, PRNG& prng, size_t num_instances, std::vector<Sender*>& senders);
            ~Sender();

            Proto send(coproto::Socket& sock, uint_fast32_t n);
            void eval(sparse_comp::point& point, size_t k, size_t n, VecMatrix<block>& out);

    };

    class Receiver {

        private:
            inline static AES aes = AES(block(7,7));
            MultiOprfRecvr* oprfRecvr;

            Receiver(MultiOprfRecvr* oprfRecvr);

        public:
            static Proto setup(coproto::Socket& sock, PRNG& prng, size_t num_instances, std::vector<Receiver*>& receivers);
            ~Receiver();

            Proto receive(coproto::Socket& sock, uint32_t n, std::vector<oprf_point>& points, std::vector<block>& outs);

    };

    

};