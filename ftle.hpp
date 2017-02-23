#pragma once

#include "field.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace LCS {

// FTLE field
template <typename T, unsigned Dim = 2>
class FTLE
{
    public:
        // constructor
        FTLE(FlowField<T, Dim>& ff): flow_field_(ff),
            nx_(std::get<0>(ff.CurrentPosition().GetAll().Size())),
            ny_(std::get<1>(ff.CurrentPosition().GetAll().Size())),
            data_(std::get<0>(ff.CurrentPosition().GetAll().Size()),
                  std::get<1>(ff.CurrentPosition().GetAll().Size())),
            initial_time_(ff.InitialPosition().GetTime()),
            current_time_(ff.GetTime()) {}

        inline void SetAll(LCS::Tensor<T, Dim>& data)
        {
            data_ = data;
        }

        // getter
        inline T Get(const unsigned i, const unsigned j) const
        {
            return data_(i,j);
        }

        inline auto& GetAll() const
        {
            return data_;
        }

        // calculate FTLE
        void Calculate()
        {
            Vector<T, Dim> x_pre, x_next, y_pre, y_next;
            Vector<T, Dim> x0_pre, x0_next, y0_pre, y0_next;

            auto initial_pos_data = flow_field_.InitialPosition().GetAll();
            auto current_pos_data = flow_field_.CurrentPosition().GetAll();
            auto dt = current_time_ - initial_time_;

            Eigen::Matrix<T, 2, 2> deformation, ftle_mat;

            for (unsigned i = 0; i < nx_; ++i)
            {
                for (unsigned j = 0; j < ny_; ++j)
                {
                    std::tie(x0_pre, x0_next, y0_pre, y0_next) = initial_pos_data.GetNearby(i,j);
                    std::tie(x_pre, x_next, y_pre, y_next) = current_pos_data.GetNearby(i,j);

                    // deformation tensor 
                    deformation(0,0) = (x_next.x-x_pre.x) / (x0_next.x-x0_pre.x);
                    deformation(0,1) = (y_next.x-y_pre.x) / (y0_next.y-y0_pre.y);
                    deformation(1,0) = (x_next.y-x_pre.y) / (x0_next.x-x0_pre.x);
                    deformation(1,1) = (y_next.y-y_pre.y) / (y0_next.y-y0_pre.y);

                    ftle_mat = deformation.transpose() * deformation;

                    auto eivals = ftle_mat.template selfadjointView<Eigen::Lower>().eigenvalues();
                    data_(i,j) = .5 * std::log(eivals.maxCoeff()) / dt;
                }
            }
        }

        inline void WriteToFile(const std::string& file_name) const
        {
            std::ofstream file;
            file.open(file_name);
            if (!file.is_open())
                throw std::runtime_error("file does not open correctly!");

            file.clear();
            file << nx_ << std::endl;
            file << ny_ << std::endl;
            file << current_time_ << std::endl;
            file << data_;
            file.close();
        }

    private:
        LCS::Tensor<T, Dim> data_;
        const unsigned nx_;
        const unsigned ny_;
        T initial_time_;
        T current_time_;
        FlowField<T, Dim>& flow_field_;
};

}
