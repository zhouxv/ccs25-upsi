#include "./VecMatrix.h"
#include <vector>
#include <algorithm>

template <class T>
sparse_comp::VecMatrix<T>::VecMatrix(size_t row_count, size_t column_count) {
    this->m_row_count = row_count;
    this->m_col_count = column_count;
    this->rows = new vector<vector<T>*>(row_count);

    for(size_t i=0;i < row_count;i++) {
        this->rows->at(i) = new vector<T>(column_count);
    }
};

template <class T>
sparse_comp::VecMatrix<T>::~VecMatrix() {
    for(size_t i=0;i < this->m_row_count;i++) {
        delete (this->rows->at(i));
    }

    delete (this->rows);
};


template <class T>
size_t sparse_comp::VecMatrix<T>::row_count() {
    return this->m_row_count;
}

template <class T>
size_t sparse_comp::VecMatrix<T>::col_count() {
    return this->m_col_count;
}

template <typename T>
std::vector<T>& sparse_comp::VecMatrix<T>::row(size_t row_idx) {
    std::vector<T>* row = this->rows->at(row_idx);
    return (*row);
};

template <typename T>
std::vector<T>& sparse_comp::VecMatrix<T>::operator[](size_t i) {
    return this->row(i);
};

template <class T>
void sparse_comp::VecMatrix<T>::cshift(size_t row_idx, size_t offset) {
    vector<T>& row = *(this->rows->at(row_idx)); 

    std::rotate(row.begin(),row.begin()+offset,row.end());
    // std::rotate(row.rbegin(), row.rbegin() + offset, row.rend());
    // std::rotate(row.rbegin(), row.rbegin() - offset, row.rend());
};

template <class T>
void sparse_comp::VecMatrix<T>::cshift(vector<size_t>& offsets) {
    for(size_t i=0;i < offsets.size();i++) {
        size_t offset = offsets[i];
        vector<T>& row = *(this->rows[i]);

        return std::rotate(row.begin(),row.begin()+row.size()-offset,row.end());
    }
};

