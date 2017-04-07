#pragma once

#include <iostream>

namespace LCS {

// outputs

template <typename T>
std::ostream& operator<< (std::ostream& os, const Vector<T, 2>& vec)
{
    os << vec.x << std::endl << vec.y;

    return os;
}

template <typename T>
std::ostream& operator<< (std::ostream& os, const Scalar<T>& scalar)
{
    os << scalar.value;

    return os;
}

template <typename T, unsigned Dim>
std::ostream& operator<< (std::ostream& os, const Tensor<T, Dim>& tensor)
{
    unsigned nx, ny;
    std::tie(nx, ny) = tensor.Size();

    for (unsigned i = 0; i < nx; ++i)
        for (unsigned j = 0; j < ny; ++j)
            os << tensor(i, j) << std::endl;

    return os;
}

template <typename T, unsigned Dim, unsigned Size>
std::ostream& operator<< (std::ostream& os, const Field<T, Dim, Size>& field)
{
    os << field.GetAll();

    return os;
}

// inputs

template <typename T>
std::istream& operator>> (std::istream& is, Vector<T, 2>& vec)
{
    is >> vec.x;
    is >> vec.y;

    return is;
}

template <typename T>
std::ostream& operator>> (std::ostream& is, Scalar<T>& scalar)
{
    is >> scalar.value;

    return is;
}

template <typename T, unsigned Dim>
std::istream& operator>> (std::istream& is, Tensor<T, Dim>& tensor)
{
    unsigned nx, ny;
    std::tie(nx, ny) = tensor.Size();

    for (unsigned i = 0; i < nx; ++i)
        for (unsigned j = 0; j < ny; ++j)
            is >> tensor(i, j);

    return is;
}

template <typename T, unsigned Dim, unsigned Size>
std::istream& operator>> (std::istream& is, Field<T, Dim>& field)
{
    is >> field.GetAll();

    return is;
}

}
