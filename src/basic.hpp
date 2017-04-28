#pragma once

#include <vector>
#include <algorithm>
#include <tuple>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <string>
#include <fstream>
#include <functional>
#include <type_traits>
#include <chrono>
#include <iomanip>
#include <omp.h>

namespace LCS {

// elementwise addition for STL vectors
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::plus<T>());
    return result;
}

// elementwise substraction for STL vectors
template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::minus<T>());
    return result;
}

// multiplication of a scalar with STL vector
template <class T, class T2>
std::vector<T> operator* (const T2 c, std::vector<T> a)
{
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] = c * a[i];

    return a;
}

// vector
template <typename T, unsigned Size = 2>
struct Vector
{
    Vector(): x(), y() {}
    Vector(const T& x, const T& y): x(x), y(y) {}

    T x, y;
};

template <typename T, unsigned Size>
Vector<T, Size> operator+(const Vector<T, Size>& a, const Vector<T, Size>& b)
{
    Vector<T,Size> result(a.x+b.x, a.y+b.y);
    return result;
}

template <typename T, unsigned Size>
Vector<T, Size> operator-(const Vector<T, Size>& a, const Vector<T, Size>& b)
{
    Vector<T,Size> result(a.x-b.x, a.y-b.y);
    return result;
}

template <typename T, unsigned Size>
Vector<T, Size> operator*(const T c, Vector<T, Size> a)
{
    a.x *= c;
    a.y *= c;
    return a;
}

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
            data_(nx * ny, T()), nx_(nx), ny_(ny) {}

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

        inline void SetAll(std::vector<T>& data)
        {
            data_ = data;
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

template <typename Field, typename T>
void interpolate(T x1, T x2, Field& f1, Field& f2, T xm, Field& result)
{
    auto tensor_result = f1.GetAll();
    if (x1 != x2)
    {
        auto vec1 = f1.GetAll().GetAll();
        auto vec2 = f2.GetAll().GetAll();
        auto vec_result = vec1 + ((xm - x1)/(x2 - x1)) * (vec2 - vec1);

        tensor_result.SetAll(vec_result);
    }
    result.SetAll(tensor_result);
}

// clock for time recording
class Clock
{
    public:
        // constructor
        Clock():
            total_elapsed_time_(0), begin_time_(), end_time_(), started_(false) {}

        inline void Begin()
        {
            // TODO: make sure clock is not started
            if (!started_)
            {
                begin_time_ = std::chrono::system_clock::now();
                started_ = true;
            }
        }

        inline void End()
        {
            // TODO: make sure clock is started
            if (started_)
            {
                end_time_ = std::chrono::system_clock::now();
                started_ = false;

                // calculate elapsed time (in seconds) and add it to total_time_
                elapsed_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>
                    (end_time_ - begin_time_).count() / 1e9;
                total_elapsed_time_ += elapsed_time_;
            }
        }

        inline double GetElapsedTime() const
        {
            return elapsed_time_;
        }

        inline double GetTotalElapsedTime() const
        {
            return total_elapsed_time_;
        }

    private:
        double total_elapsed_time_, elapsed_time_;
        std::chrono::time_point<std::chrono::system_clock> begin_time_, end_time_;
        bool started_; // if the clock is started
};


}
