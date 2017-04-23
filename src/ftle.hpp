#pragma once

#include "flow.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace LCS {

// FTLE field
template <typename T, unsigned Dim = 2>
class FTLE : public Field<T, Dim, 1>
{
    public:
        // constructor
        FTLE(FlowField<T, Dim>& ff):
            Field<T, Dim, 1>(std::get<0>(ff.CurrentPosition().GetAll().Size()),
                std::get<1>(ff.CurrentPosition().GetAll().Size())),
            flow_field_(ff), 
            initial_time_(ff.InitialPosition().GetTime())
        {
            this->UpdateTime(ff.GetTime());
        }

        // getter
        inline T Get(const unsigned i, const unsigned j) const
        {
            return this->data_(i,j);
        }

        // calculate FTLE
        void Calculate()
        {
            auto initial_pos_data = flow_field_.InitialPosition().GetAll();
            auto current_pos_data = flow_field_.CurrentPosition().GetAll();
            auto dt = this->time_ - initial_time_;

            std::cout << "FTLE calculation begins" << std::endl;

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

            std::cout << "FTLE calculation ends" << std::endl;
        }

    private:
        FlowField<T, Dim>& flow_field_;
        T initial_time_;
};

}
