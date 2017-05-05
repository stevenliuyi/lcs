/**
    @file velocity_function.hpp
    @brief Demo analytic velocity fields.

    This file contains implementations of demo velocity fields with analytic velocity functions.
    */
#pragma once

#include <cmath>
#include <tuple>
#include <cassert>

namespace LCS {
/**
    @namespace LCS::VelocityFunction
    @brief Namespace for demo velocity functions.

    All demo analytic velocity functions are included in this namespace.
    */
namespace VelocityFunction {

/** @brief Bower model for meandering jet

    The velocity function of the Bower model for meandering jets is implemented here. @cite bower1991

    The streamfunction is defined as
    \f[
        \Psi(x,y,t)=\Psi_0\Big[1-\tanh\Big(\frac{y-y_c}{\lambda/\cos\alpha}\Big)\Big],
    \f]
    where \f$\Psi\f$ is the scale factor, \f$ y_c=A\sin[k(x-c_xt)]\f$ is the center streamline, \f$A\f$ is the wave amplitude, \f$k=2\pi/L\f$ is the wave number, \f$\lambda\f$ is the scale width of the jet, and \f$\alpha=\tan^{-1}\{Ak\cos[k(x-c_xt)]\}\f$ is the direction of current.

    Here, we implement the streamfunction in the moving frame, which is independent of time:
    \f[
        \Psi'(x',y')=\Psi_0\Big[1-\tanh\Big(\frac{y'-y_c'}{\lambda/\cos\alpha'}\Big)\Big],
    \f]
    where \f$y_c'=A\sin(kx')\f$, and \f$\alpha'=\tan^{-1}[Ak\cos(kx')]\f$.

    Then, the velocity in the moving frame can be obtained by calculating the partial derivatives of the streamfunction:
    \f[
        u'=-\frac{\partial\Psi'}{\partial y'},\quad v'=\frac{\partial\Psi'}{\partial x'}.
    \f]

    @tparam T Numeric data type of the values.
    */
template <typename T>
struct BowerModel
{
    /** Constructor with default parameters. */
    BowerModel(): sc_(50), a_(50), l_(400), cx_(10), lambda_(40) {}

    /** Constructor with customized parameters.
        @param p A vector of parameters.
        */
    BowerModel(std::vector<T>& p)
    {
        assert(p.size() == 5);
        sc_ = p[0]; a_ = p[1]; l_ = p[2]; cx_ = p[3]; lambda_ = p[4];
    }

    /** Get the velocity at a given point.
        @param x \f$x\f$-coordinate of the point.
        @param y \f$y\f$-coordinate of the point.
        @param t Time for velocity. Since the flow is steady (independent of time), `t` should be 0 here.
        */
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

    T sc_; /**<Magnitude of downstream speed at the jet center (km/day).*/
    T a_; /**<Wave amplitude (km).*/
    T l_; /**<Wave length (km).*/
    T cx_; /**<Jet phase speed (km/day).*/
    T lambda_; /**<Scale width of the jet (40km).*/
};


/** @brief Double-gyre model
    
    The velocity function of the double-gyre model is implemented here. The model  consists of a pair of counter-rotating gyres. @cite shadden2005

    The velocity for a given point \f$(x,y)\f$ at time \f$t\f$ is given by
    \f[
        u = -\pi A\sin(\pi f)\cos(\pi y),
    \f]
    and
    \f[
        v = \pi A\cos(\pi f)\sin(\pi y)\frac{\partial f}{\partial x},
    \f]
    where \f$f(x,t)=ax^2+bx\f$, \f$a(t)=\epsilon\sin(\omega t)\f$ and \f$b(t)=1-2\epsilon\sin(\omega t)\f$.

    @tparam T Numeric data type of the values.
    */
template <typename T>
struct DoubleGyreModel
{
    /** Constructor with default parameters. */
    DoubleGyreModel(): epsilon_(0.1), a_(0.1), omega_(4*std::atan(1)/5) {}

    /** Constructor with customized parameters.
        @param p A vector of parameters.
        */
    DoubleGyreModel(std::vector<T>& p)
    {
        assert(p.size() == 3);
        epsilon_ = p[0]; a_ = p[1]; omega_ = p[2];
    }

    /** Get the velocity for a given point at a given time
        @param x \f$x\f$-coordinate of the point.
        @param y \f$y\f$-coordinate of the point.
        @param t Time for velocity.
        */
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

    T epsilon_; /**<Parameter \f$\epsilon\f$ for the model.*/
    T a_; /**<Parameter \f$A\f$ for the model.*/
    T omega_; /**<Parameter \f$\omega\f$ for the model.*/

};


}
}
