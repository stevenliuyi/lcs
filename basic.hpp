#pragma once

#include <vector>
#include <algorithm>
#include <tuple>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <string>
#include <fstream>

namespace LCS {

// vector
template <typename T, unsigned Size = 2>
struct Vector
{
    Vector(): x(), y() {}
    Vector(const T& x, const T& y): x(x), y(y) {}

    T x, y;
};

// scalar
template <typename T>
using Scalar = Vector<T, 1>;

template <typename T>
struct Vector<T, 1>
{
    Vector(): value() {}
    Vector(const T& x): value(x) {}

    T value;
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

        inline T& Get(unsigned i, unsigned j)
        {
            return data_[i*ny_ + j];
        }

        inline T Get(unsigned i, unsigned j) const
        {
            return data_[i*ny_ + j];
        }

        // for bool type
        inline T GetValue(unsigned i, unsigned j) const
        {
            return data_[i*ny_ + j];
        }

        inline auto& GetAll() const
        {
            return data_;
        }

        // obtain nearby points of a given point
        inline auto GetNearby(const unsigned i, const unsigned j) const
        {
            auto x_pre = (i != 0) ? Get(i-1,j) : Get(i,j);
            auto x_next = (i != nx_-1) ? Get(i+1,j) : Get(i,j);

            auto y_pre = (j != 0) ? Get(i,j-1) : Get(i,j);
            auto y_next = (j != ny_-1) ? Get(i,j+1) : Get(i,j);

            return std::make_tuple(x_pre, x_next, y_pre, y_next);

        }

        inline void SetValue(unsigned i, unsigned j, T value)
        {
            data_[i*ny_ + j] = value;
        }

        inline auto Size() const
        {
            return std::make_tuple(nx_, ny_);
        }

    private:
        std::vector<T> data_;
        const unsigned nx_;
        const unsigned ny_;
};

// interpolation
template <typename T>
inline T interpolate(T x1, T x2, T y1, T y2, T xm)
{
    return y1 + (xm - x1) * (y2 - y1) / (x2 - x1);
}

}
