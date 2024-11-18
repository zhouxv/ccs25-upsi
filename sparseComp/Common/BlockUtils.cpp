#include "./BlockUtils.h"
#include <cstdint>
#include <vector>

using namespace std;

template <typename T>
using vec = std::vector<T>;

void sparse_comp::block_vec_xor(vec<block>& va, vec<block>& vb, vec<block>& vc) {
    for(uint32_t i=0;i < vc.size();i++) {
        vc[i] = va[i] ^ vb[i];
    }
}