#pragma once

#include <iostream>

namespace LCS {

// outputs

/** Output a Vector.
    @param os Output stream object.
    @param vec Vector to be outputted.
    @return Updated output stream object.
    */
template <typename T>
std::ostream& operator<< (std::ostream& os, const Vector<T, 2>& vec)
{
    os << vec.x << std::endl << vec.y;

    return os;
}

/** Output a Scalar.
    @param os Output stream object.
    @param scalar Scalar to be outputted.
    @return Updated output stream object.
    */
template <typename T>
std::ostream& operator<< (std::ostream& os, const Scalar<T>& scalar)
{
    os << scalar.value;

    return os;
}

/** Output a Tensor.
    @param os Output stream object.
    @param tensor Tensor to be outputted.
    @return Updated output stream object.
    */
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

/** Output a Field.
    @param os Output stream object.
    @param field Field to be outputted.
    @return Updated output stream object.
    */
template <typename T, unsigned Dim, unsigned Size>
std::ostream& operator<< (std::ostream& os, const Field<T, Dim, Size>& field)
{
    os << field.GetAll();

    return os;
}

// inputs

/** Input a Vector.
    @param is Input stream object.
    @param vec Vector to be outputted.
    @return Updated input stream object.
    */
template <typename T>
std::istream& operator>> (std::istream& is, Vector<T, 2>& vec)
{
    is >> vec.x;
    is >> vec.y;

    return is;
}

/** Input a Scalar.
    @param is Input stream object.
    @param scalar Scalar to be outputted.
    @return Updated input stream object.
    */
template <typename T>
std::ostream& operator>> (std::ostream& is, Scalar<T>& scalar)
{
    is >> scalar.value;

    return is;
}

/** Input a Tensor.
    @param is Input stream object.
    @param tensor Tensor to be outputted.
    @return Updated input stream object.
    */
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

/** Input a Field.
    @param is Input stream object.
    @param field Field to be outputted.
    @return Updated input stream object.
    */
template <typename T, unsigned Dim, unsigned Size>
std::istream& operator>> (std::istream& is, Field<T, Dim>& field)
{
    is >> field.GetAll();

    return is;
}

}
