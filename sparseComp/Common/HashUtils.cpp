#include "./HashUtils.h"
#include <iostream> 

block sparse_comp::hash_point(const AES& aes, point pot) {

    const size_t POINT_MAX_DIM = point::MAX_DIM;

    const uint8_t padded_coords_len = 4*((POINT_MAX_DIM/4) + 1);

    // Instantiate a vector of length divisible by 4 and bigger than POINT_MAX_DIM.
    uint32_t padded_coords[padded_coords_len] = {}; // All zero initialization.
    block blocksToHash[(POINT_MAX_DIM/4) + 1];

    for (uint8_t i=0;i < POINT_MAX_DIM;i++) {
        padded_coords[i] = pot.coords[i];
    }

    for(uint8_t i=0;i < padded_coords_len;i = i + 4) {
        
        uint64_t bp1 = (((uint64_t) padded_coords[i]) << 32) + ((uint64_t) padded_coords[i+1]);
        
        uint64_t bp2 = (((uint64_t) padded_coords[i+2]) << 32) + ((uint64_t) padded_coords[i+3]);
        
        blocksToHash[i/4] = block(bp1,bp2);

    }

    block hashes[(POINT_MAX_DIM/4) + 1];

    aes.hashBlocks(blocksToHash,(POINT_MAX_DIM/4) + 1, hashes);

    block digest = block(0,0);

    for(uint8_t i=0;i < (POINT_MAX_DIM/4) + 1;i++) {
        digest = digest ^ hashes[i]; 
    }

    return digest;

}

block sparse_comp::hash_point(const AES& aes, point point, size_t sot_idx, size_t msg_vec_idx) {

    block extra_block = block(sot_idx,msg_vec_idx);

    block extra_block_hash = aes.hashBlock(extra_block);

    block digest = extra_block_hash ^ hash_point(aes, point);

    return digest;

}