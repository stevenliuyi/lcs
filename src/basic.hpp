/**
    @file basic.hpp
    @brief Basic classes needed such as Vector and Tensor.
    
    This file contains the basic classes nedded in the project, including extensions of `std::vector` class to support vector arithemtic operations, functions to perfrom linear interpolations, and new Vecotr and Tensor classes to represent physical vectors (including scalars) and tensors (including matrices. It also contains elper class Clock to record program running time.
    */
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

/**
    @namespace LCS
    @brief Namespace for the whole project.

    All classes in this project are in this namespace.
    */
namespace LCS {

/** Elementwise addition for STL vectors.
    @tparam T Numeric data type of the elements.
    @param a The first STL vectors in the addition.
    @param b The second STL vectors in the addition.
    @return A new STL vector that is the elementwise addition of `a` and `b` (`a` + `b`).
    */
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

/** Elementwise subtraction for STL vectors.
    @tparam T Numeric data type of the elements.
    @param a The first STL vectors in the subtraction.
    @param b The second STL vectors in the subtraction.
    @return A new STL vector that is the elementwise subtraction of `a` and `b` (`a` - `b`).
    */
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

/** Multiplication of a scalar with a STL vector.
    @tparam T Numeric data type of the vector elements.
    @tparam T2 Numeric data tpye of the scalar value.
    @param c A scalar value.
    @param a A STL vector.
    @return A new STL vector that is the result of `c` * `a`.
    */
template <class T, class T2>
std::vector<T> operator* (const T2 c, std::vector<T> a)
{
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] = c * a[i];

    return a;
}

/** @brief Struct for vectors.

    This is a struct to represent vectors. It is for vectors in the physical sense, not the STL vectors which are basically dynamic size arrays.
    @tparam T Numeric data type of the vector elements.
    @tparam Size Length of the vector.
    */
template <typename T, unsigned Size = 2>
struct Vector
{
    /** Default constructor.  */
    Vector(): x(), y() {}

    /** Constructor for initializating the vector.
        @param x \f$x\f$-component of the vector.
        @param y \f$y\f$-component of the vector.
        */
    Vector(const T& x, const T& y): x(x), y(y) {}

    T x;/**<\f$x\f$-component of the vector.*/
    T y;/**<\f$y\f$-component of the vector.*/
};

/** Implementation of vector additions.
    @tparam T Numeric data type of the vector elements.
    @tparam Size Length of the vector.
    @param a The first Vector in the addition.
    @param b The second Vector in the addition.
    @return A Vector that is the result of `a` + `b`.
    */
template <typename T, unsigned Size>
Vector<T, Size> operator+(const Vector<T, Size>& a, const Vector<T, Size>& b)
{
    Vector<T,Size> result(a.x+b.x, a.y+b.y);
    return result;
}

/** Implementation of vector subtractions.
    @tparam T Numeric data type of the vector elements.
    @tparam Size Length of the vector.
    @param a The first Vector in the subtraction.
    @param b The second Vector in the subtraction.
    @return A Vector that is the result of `a` - `b`.
    */
template <typename T, unsigned Size>
Vector<T, Size> operator-(const Vector<T, Size>& a, const Vector<T, Size>& b)
{
    Vector<T,Size> result(a.x-b.x, a.y-b.y);
    return result;
}

/** Implementation of multiplication of a scalar and a vector.
    @tparam T Numeric data type of the vector elements.
    @tparam Size Length of the vector.
    @param c A scalar value.
    @param a A Vector.
    @return A Vector that is the result of `c` * `a`.
    */
template <typename T, unsigned Size>
Vector<T, Size> operator*(const T c, Vector<T, Size> a)
{
    a.x *= c;
    a.y *= c;
    return a;
}

/** Define a scalar as a Vector with one element.
    @tparam T Numeric data type of the element.
    */
template <typename T>
using Scalar = Vector<T, 1>;

/** @brief Vector with one element.

    This is a specialization of Vector with one element, which is actually a scalar.

    @tparam T Numeric data type of the element.
    */
template <typename T>
struct Vector<T, 1>
{
    /** Default constructor.  */
    Vector(): value() {}

    /** Constructor for initializating the scalar.
        @param x Value of the scalar.
        */
    Vector(const T& x): value(x) {}

    T value; /**<Scalar value.*/
};

/** @brief Class for tensors.
    
    This class is used for representing tensors (matrices) in a 1D vector.
    @tparam T Data type of each element in the tensor.
    @tparam Dim Dimension of the tensor.
    */
template <typename T, unsigned Dim = 2>
class Tensor
{
    public:
        /** Constructor for initializing the tensor.
            @param nx The number of grid points in \f$x\f$-direction.
            @param ny The number of grid points in \f$y\f$-direction.
            */
        Tensor(unsigned nx, unsigned ny):
            data_(nx * ny, T()), nx_(nx), ny_(ny) {}

        /** Class assignment.
            @param t A Tensor that needs to be assigned to this Tensor.
            */
        inline void operator= (const Tensor<T, Dim>& t)
        {
            // need to implement size check
            data_ = t.GetAll();
        }

        /** Get one element reference of the tensor using the () operator.
            @param i Index in \f$x\f$-coordinate.
            @param j Index in \f$y\f$-coordinate.
            @return Reference of the element at the given location.
            */
        inline T& operator() (unsigned i, unsigned j)
        {
            // need to implement size check
            return data_[i*ny_ + j];
        }

        /** Get one element of the tensor using the () operator.
            @param i Index in \f$x\f$-coordinate.
            @param j Index in \f$y\f$-coordinate.
            @return The element at the given location.
            */
        inline T operator() (unsigned i, unsigned j) const
        {
            return data_[i*ny_ + j];
        }

        /** Get one element reference of the tensor.
            @param i Index in \f$x\f$-coordinate.
            @param j Index in \f$y\f$-coordinate.
            @return Reference of the element at the given location.
            */
        inline T& Get(unsigned i, unsigned j)
        {
            return data_[i*ny_ + j];
        }

        /** Get one element of the tensor.
            @param i Index in \f$x\f$-coordinate.
            @param j Index in \f$y\f$-coordinate.
            @return The element at the given location.
            */
        inline T Get(unsigned i, unsigned j) const
        {
            return data_[i*ny_ + j];
        }

        /** An another implementation of getting one tensor element. The Get() function does not work well sometimes, perhaps due to the overloading.
            @param i Index in \f$x\f$-coordinate.
            @param j Index in \f$y\f$-coordinate.
            @return The element at the given location.
            */
        inline T GetValue(unsigned i, unsigned j) const
        {
            return data_[i*ny_ + j];
        }

        /** Get all elements of the tensor in STL vector format
            @return STL vector that contains raw data of the tensor
            */
        inline auto& GetAll() const
        {
            return data_;
        }

        /** Obtain nearby points of a given point
            @param i Index of the given point in \f$x\f$-coordinate.
            @param j Index of the given point in \f$y\f$-coordinate.
            @return Tuple that contains the \f$x\f$- and \f$y\f$-coordinates of the nearby points.
            */
        inline auto GetNearby(const unsigned i, const unsigned j) const
        {
            auto x_pre = (i != 0) ? Get(i-1,j) : Get(i,j);
            auto x_next = (i != nx_-1) ? Get(i+1,j) : Get(i,j);

            auto y_pre = (j != 0) ? Get(i,j-1) : Get(i,j);
            auto y_next = (j != ny_-1) ? Get(i,j+1) : Get(i,j);

            return std::make_tuple(x_pre, x_next, y_pre, y_next);

        }

        /** Set value for one tensor element.
            @param i Index in \f$x\f$-coordinate.
            @param j Index in \f$y\f$-coordinate.
            @param value The value to be assigned.
            */
        inline void SetValue(unsigned i, unsigned j, T value)
        {
            data_[i*ny_ + j] = value;
        }

        /** Set values of all tensor elements.
            @param data A vector contains all element values of the tensor.
            */
        inline void SetAll(std::vector<T>& data)
        {
            data_ = data;
        }

        /** Get the shape of the tensor.
            @return Tuple that contains the shape of the tensor.
            */
        inline auto Size() const
        {
            return std::make_tuple(nx_, ny_);
        }
    
    private:
        std::vector<T> data_; /**<STL vector that contains the raw values of the tensor.*/
        const unsigned nx_; /**The number of grid points in \f$x\f$-direction.*/
        const unsigned ny_; /**The number of grid points in \f$y\f$-direction.*/
};


/** Linear interpolation of two known points \f$(x_1,y_1)\f$ and \f$(x_2,y_2)\f$.
    @tparam T Data type of the values.
    @param x1 \f$x\f$-coordinate of the first known point.
    @param x2 \f$x\f$-coordinate of the second known point.
    @param y1 \f$y\f$-coordinate of the first known point.
    @param y2 \f$y\f$-coordinate of the second known point.
    @param xm \f$x\f$-coordinate of the interpolated point.
    @return \f$y\f$-coordinate of the interpolated point.
    */
template <typename T>
inline T interpolate(T x1, T x2, T y1, T y2, T xm)
{
    return y1 + (xm - x1) * (y2 - y1) / (x2 - x1);
}

/** Temporal linear interpolation of two Field.
    @tparam T Data type of the time.
    @param x1 First time point.
    @param x2 Second time point.
    @param f1 Field at time `x1`.
    @param f2 Field at time `x2`.
    @param xm Time that associated with the interpolated Field.
    @return Interpolated Field of `f1` and `f2`.
    */
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

/** @brief Class for time recording.

    This class utilizes `std::chrono::system_clock` to record the time. It is used for calculating the code running times.
    */
class Clock
{
    public:
        /** Constructor for initializating the clock.
            */
        Clock():
            total_elapsed_time_(0), begin_time_(), end_time_(), started_(false) {}

        /** Start the clock.  */
        inline void Begin()
        {
            // TODO: make sure clock is not started
            if (!started_)
            {
                begin_time_ = std::chrono::system_clock::now();
                started_ = true;
            }
        }

        /** Stop the clock.  */
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

        /** Get the elapsed time of the last run of the clock.
            @return elapsed time in seconds.
            */
        inline double GetElapsedTime() const
        {
            return elapsed_time_;
        }

        /** Get the total elapsed time of all runs of the clock.
            @return Total elapsed time in seconds.
            */
        inline double GetTotalElapsedTime() const
        {
            return total_elapsed_time_;
        }

    private:
        double total_elapsed_time_;/**<Total elapsed time in seconds.*/
        double elapsed_time_;/**<Elapsed time of the last run in seconds.*/
        std::chrono::time_point<std::chrono::system_clock> begin_time_;/**<Time when the clock begins.*/
        std::chrono::time_point<std::chrono::system_clock> end_time_;/**<Time when the clock stops.*/
        bool started_; /**<If the clock is started.*/
};


}
