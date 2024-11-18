#pragma once

#include <vector>

template <typename T>
using vector = std::vector<T>;

namespace sparse_comp {

    template <typename T>
    class VecMatrix {
        
        private:
            size_t m_row_count;
            size_t m_col_count;
            vector<vector<T>*>* rows;

        public:
            VecMatrix(size_t row_count, size_t column_count);
            ~VecMatrix();

            vector<T>& row(size_t row_idx);
            vector<T>& operator [](size_t i);
            void cshift(size_t row_idx, size_t offset);
            void cshift(vector<size_t>& offsets);
            size_t row_count();
            size_t col_count();

            // Let M be a matrix and v a vector. This method returns a new matrix M' where M'[i][j] = M[i][j] - v[i].
            // VecMatrix<T>* subVec(std::vector<T> v);
    };

};

#include "./VecMatrix.cpp"