#include "./Common.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Crypto/AES.h"

using block = osuCrypto::block;
using AES = osuCrypto::AES;
using point = sparse_comp::point;

namespace sparse_comp {

    block hash_point(const AES& aes, point point);
    block hash_point(const AES& aes, point point, size_t sot_idx, size_t msg_vec_idx);

};