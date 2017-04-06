#pragma once

#include <cmath>
#include <tuple>
#include <cassert>

// analytic velocity function examples

namespace LCS {
namespace VelocityFunction {

// Bower model for meandering jet
// see "A Simple Kinematic Mechanism for Mixing Fluid Parcels across a Meandering Jet (Bower, 1991)
template <typename T>
struct BowerModel
{
    // constructor with default parameters
    BowerModel(): sc_(50), a_(50), l_(400), cx_(10), lambda_(40) {}

    // constructor with customized parameters
    BowerModel(std::vector<T>& p)
    {
        assert(p.size() == 5);
        sc_ = p[0]; a_ = p[1]; l_ = p[2]; cx_ = p[3]; lambda_ = p[4];
    }

    // steady flow, t doesn't matter
    inline auto operator() (const T x, const T y, const T t = 0) const
    {

        T phi0 = sc_ * lambda_; // scale factor
        T k = 8 * std::atan(1) / l_; // wave number (2*pi/L)

        T yc = a_ * std::sin(k * x);
        T dyc = a_ * k * std::cos(k * x);
        T alpha0 = lambda_ * std::sqrt(dyc*dyc + 1);

        T u = -cx_ + phi0 / std::pow(std::cosh((y-yc)/alpha0), 2) / alpha0;
        T v = -phi0 * ((yc*dyc*k*k*(y-yc)) / (lambda_*std::pow(dyc*dyc+1,1.5)) - dyc/alpha0)
            / std::pow(std::cosh((y-yc)/alpha0), 2);

        return std::make_tuple(u, v);
    }

    T sc_; // magnitude of downstream speed at the jet center (km/day)
    T a_; // wave amplitude (km)
    T l_; // wave length (km)
    T cx_; // jet phase speed (km/day)
    T lambda_; // scale width of the jet (40km)
};


// double gypre model that consists of a pair of counter-rotating gyres
template <typename T>
struct DoubleGyreModel
{
    // constructor with default parameters
    DoubleGyreModel(): epsilon_(0.1), a_(0.1), omega_(4*std::atan(1)/5) {}

    // constructor with customized parameters
    DoubleGyreModel(std::vector<T>& p)
    {
        assert(p.size() == 3);
        epsilon_ = p[0]; a_ = p[1]; omega_ = p[2];
    }

    inline auto operator() (const T x, const T y, const T t) const
    {
        T at = epsilon_ * std::sin(omega_ * t);
        T bt = 1 - 2 * epsilon_ * std::sin(omega_ * t);
        T f = at*x*x + bt*x;
        T dfdx = 2*at*x + bt;

        T pi = 4 * std::atan(1);
        T u = -pi * a_ * std::sin(pi*f) * std::cos(pi*y);
        T v = pi * a_ * std::cos(pi*f) * std::sin(pi*y) * dfdx;

        return std::make_tuple(u, v);
    }

    T epsilon_;
    T a_;
    T omega_;

};


}
}
