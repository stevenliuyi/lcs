#pragma once

#include <vector>
#include <algorithm>
#include <tuple>
#include <cassert>
#include <stdexcept>
#include <memory>

namespace LCS {

// vector
template <typename T, unsigned Dim = 2>
struct Vector
{
    Vector(): x(), y() {}
    Vector(const T& x, const T& y): x(x), y(y) {}

    T x, y;
};

// tensor
// represent tensor (matrix) in a 1D vector
template <typename T, unsigned Dim = 2>
class Tensor
{
    public:
        // constructor 
        Tensor(unsigned nx, unsigned ny):
            nx_(nx), ny_(ny), data_(nx * ny, T()) {}

        // assignment
        inline void operator= (const Tensor<T, Dim>& t)
        {
            // need to implement size check
            data_ = t.GetAll();
        }

        inline T& operator() (unsigned i, unsigned j)
        {
            // need to implement size check
            return data_[i*ny_ + j];
        }

        inline T operator() (unsigned i, unsigned j) const
        {
            return data_[i*ny_ + j];
        }

        inline auto GetAll() const
        {
            return data_;
        }

    private:
        std::vector<T> data_;
        const unsigned nx_;
        const unsigned ny_;
};

}
