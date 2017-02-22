#pragma once

#include <iostream>

namespace LCS {

// outputs

template <typename T, unsigned Dim = 2>
std::ostream& operator<< (std::ostream& os, const Vector<T, Dim>& vec)
{
    os << vec.x << std::endl;
    os << vec.y << std::endl;

    return os;
}

template <typename T, unsigned Dim = 2>
std::ostream& operator<< (std::ostream& os, const Tensor<T, Dim>& tensor)
{
    unsigned nx, ny;
    std::tie(nx, ny) = tensor.Size();

    for (unsigned i = 0; i < nx; ++i)
        for (unsigned j = 0; j < ny; ++j)
            os << tensor(i, j);

    return os;
}

template <typename T, unsigned Dim = 2>
std::ostream& operator<< (std::ostream& os, const Position<T, Dim>& pos)
{
    os << pos.GetAll();

    return os;
}

template <typename T, unsigned Dim = 2>
std::ostream& operator<< (std::ostream& os, const Velocity<T, Dim>& vel)
{
    os << vel.GetAll();

    return os;
}

// inputs

template <typename T, unsigned Dim = 2>
std::istream& operator>> (std::istream& is, Vector<T, Dim>& vec) {
    is >> vec.x;
    is >> vec.y;

    return is;
}

template <typename T, unsigned Dim = 2>
std::istream& operator>> (std::istream& is, Tensor<T, Dim>& tensor)
{
    unsigned nx, ny;
    std::tie(nx, ny) = tensor.Size();

    for (unsigned i = 0; i < nx; ++i)
        for (unsigned j = 0; j < ny; ++j)
            is >> tensor(i, j);

    return is;
}

template <typename T, unsigned Dim = 2>
std::istream& operator>> (std::istream& is, Position<T, Dim>& pos)
{
    is >> pos.GetAll();

    return is;
}

template <typename T, unsigned Dim = 2>
std::istream& operator>> (std::istream& is, Velocity<T, Dim>& vel)
{
    is >> vel.GetAll();

    return is;
}

}
