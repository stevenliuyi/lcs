#pragma once

#include "flow.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace LCS {

/** @brief Class for finite-time Lyapunov exponent (FTLE) fields.

    This class is used for representing the finite-time Lyapunov exponent (FTLE) fields, which are scalar fields. It is a subclass of the Field class, and there is a FlowField associated with it.
    @tparam T Numeric data type of the FTLE values.
    @tparam Dim Dimension of the FTLE field.
    */
template <typename T, unsigned Dim = 2>
class FTLE : public Field<T, Dim, 1>
{
    public:
        /** Constructor for initializing the FTLE field.
            @param ff FlowField that is associated with the FTLE field.
            */
        FTLE(FlowField<T, Dim>& ff):
            Field<T, Dim, 1>(std::get<0>(ff.CurrentPosition().GetAll().Size()),
                std::get<1>(ff.CurrentPosition().GetAll().Size())),
            flow_field_(ff), 
            initial_time_(ff.InitialPosition().GetTime())
        {
            this->UpdateTime(ff.GetTime());
        }

        /** Get the FTLE value of the given point.
            @param i Index in \f$x\f$-coordinate.
            @param j Index in \f$y\f$-coordinate.
            @return FTLE value at the point.
            */
        inline T Get(const unsigned i, const unsigned j) const
        {
            return this->data_(i,j);
        }

        /** Calculate FTLE field from the flow map that is obtained from the flow field. OpenMP is used for calculation.*/
        void Calculate()
        {
            Clock clock;
            clock.Begin();
            auto initial_pos_data = flow_field_.InitialPosition().GetAll();
            auto current_pos_data = flow_field_.CurrentPosition().GetAll();
            auto dt = this->time_ - initial_time_;

            switch(flow_field_.GetDirection())
            {
                case Forward:
                    std::cout << "Forward FTLE calculation begins" << std::endl;
                    break;
                case Backward:
                    std::cout << "Backward FTLE calculation begins" << std::endl;
                    break;
                default: break;
            }

            #pragma omp parallel for
            for (unsigned i = 0; i < this->nx_; ++i)
            {
                Vector<T, Dim> x_pre, x_next, y_pre, y_next;
                Vector<T, Dim> x0_pre, x0_next, y0_pre, y0_next;

                Eigen::Matrix<T, 2, 2> deformation, cauchy_green;

                for (unsigned j = 0; j < this->ny_; ++j)
                {
                    std::tie(x0_pre, x0_next, y0_pre, y0_next) = initial_pos_data.GetNearby(i,j);
                    std::tie(x_pre, x_next, y_pre, y_next) = current_pos_data.GetNearby(i,j);

                    // deformation tensor 
                    deformation(0,0) = (x_next.x-x_pre.x) / (x0_next.x-x0_pre.x);
                    deformation(0,1) = (y_next.x-y_pre.x) / (y0_next.y-y0_pre.y);
                    deformation(1,0) = (x_next.y-x_pre.y) / (x0_next.x-x0_pre.x);
                    deformation(1,1) = (y_next.y-y_pre.y) / (y0_next.y-y0_pre.y);

                    cauchy_green = deformation.transpose() * deformation;

                    auto eivals = cauchy_green.template selfadjointView<Eigen::Lower>().eigenvalues();
                    this->data_(i,j).value = .5 * std::log(eivals.maxCoeff()) / dt;
                }
            }

            clock.End();
            switch(flow_field_.GetDirection())
            {
                case Forward:
                    std::cout << "Forward FTLE calculation ends";
                    break;
                case Backward:
                    std::cout << "Backward FTLE calculation ends";
                    break;
                default: break;
            }
            std::cout << " (Execution time: " <<
                std::setprecision(4) << clock.GetTotalElapsedTime() << "s)" << std::endl;
        }

    private:
        FlowField<T, Dim>& flow_field_;/**<FlowField that is associated with the FTLE field.*/
        T initial_time_;/**<Initial time for FTLE calculation.*/
};

}
