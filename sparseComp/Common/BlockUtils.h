#pragma once

#include "cryptoTools/Common/block.h"
#include <vector>

using block = osuCrypto::block;

template <typename T>
using vec = std::vector<T>;

namespace sparse_comp {

    // Let va and vb be two vectors of blocks. This function returns a vector of blocks vc, where vc[i] = va[i] xor vb[i]
    void block_vec_xor(std::vector<block>& va, std::vector<block>& vb, std::vector<block>& vc);

};