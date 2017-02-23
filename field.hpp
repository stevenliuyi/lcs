#pragma once

#include "basic.hpp"

namespace LCS {

template <typename T, unsigned Dim>
class Position;

template <typename T, unsigned Dim>
class Velocity;

template <typename T, typename Func, unsigned Dim>
class AnalyticVelocity;

template <typename T, unsigned Dim>
class FlowField;

template <typename T, typename Func, unsigned Dim>
class AnalyticFlowField;

// flow field
template <typename T, unsigned Dim = 2>
class FlowField
{
    public:
        // constructor
        FlowField(unsigned nx, unsigned ny):
            nx_(nx), ny_(ny), delta_(), step_(), current_time_(),
            initial_pos_(new Position<T,Dim>(nx,ny)),
            current_pos_(new Position<T,Dim>(nx,ny)) {} 

        inline auto& InitialPosition()
        {
            return *initial_pos_;
        }

        inline auto& CurrentPosition()
        {
            return *current_pos_;
        }

        inline virtual Velocity<T, Dim>& CurrentVelocity()
        {
            if (current_vel_ == nullptr)
                throw std::invalid_argument("current velocity not set!");

            return *current_vel_;
        }

        inline auto GetTime()
        {
            return current_time_;
        }

        inline void CopyInitialPositionToCurrentPosition()
        {
            auto initial_pos_data = initial_pos_->GetAll();
            current_pos_->SetAll(initial_pos_data);
            current_pos_->UpdateTime(initial_pos_->GetTime());
        }

        inline void SetDelta(const T delta)
        {
            delta_ = delta;
        }

        inline void SetStep(const unsigned step)
        {
            step_ = step;
        }

        // overloaded in AnalyticFlowField class
        virtual void SetCurrentVelocity()
        {}

        // calculate the trajectories
        void Run()
        {
            CopyInitialPositionToCurrentPosition();

            for (unsigned i = 0; i < step_; ++i)
            {
                SetCurrentVelocity();
                current_pos_->Update(*current_vel_, delta_);

                // update time
                current_time_ += delta_;
                current_pos_->UpdateTime(current_time_);
            }
        }

        template <typename A, typename B, unsigned C>
        friend class AnalyticFlowField;

    private:
        std::unique_ptr<Position<T, Dim>> initial_pos_;
        std::unique_ptr<Position<T, Dim>> current_pos_;
        std::shared_ptr<Velocity<T, Dim>> current_vel_;
        const unsigned nx_;
        const unsigned ny_;
        T delta_;
        T current_time_;
        unsigned step_;
};

template <typename T, typename Func, unsigned Dim = 2>
class AnalyticFlowField : public FlowField<T, Dim>
{
    public:
        // constructor
        AnalyticFlowField(unsigned nx, unsigned ny):
            FlowField<T, Dim>(nx, ny), parameters_() {}

        AnalyticFlowField(unsigned nx, unsigned ny, std::vector<T>& parameters):
            FlowField<T, Dim>(nx, ny), parameters_(parameters) {}

        inline void SetCurrentVelocity()
        {
            // no parameters
            if (parameters_.size() == 0)
                current_analytic_vel_.reset(new AnalyticVelocity<T, Func, Dim>
                    (this->nx_, this->ny_, *(this->current_pos_)));
            else
            // there are parameters for analytic function 
                current_analytic_vel_.reset(new AnalyticVelocity<T, Func, Dim>
                    (this->nx_, this->ny_, *(this->current_pos_), parameters_));
            
            current_analytic_vel_->UpdateTime(this->current_time_);
            this->current_vel_ = current_analytic_vel_;
        }

        inline AnalyticVelocity<T, Func, Dim>& CurrentVelocity()
        {
            if (current_analytic_vel_ == nullptr)
                throw std::invalid_argument("current velocity not set!");

            return *current_analytic_vel_;
        }

    private:
        std::shared_ptr<AnalyticVelocity<T, Func, Dim>> current_analytic_vel_;
        std::vector<T> parameters_;
};

// general field
template <typename T, unsigned Dim = 2, unsigned Size = 2>
class Field
{
    public:
        Field(unsigned nx, unsigned ny): nx_(nx), ny_(ny), time_(),
            data_(nx, ny) {}

        // getter
        inline auto& GetAll() const
        {
            return data_;
        }

        inline auto GetTime() const
        {
            return time_;
        }

        // setter
        inline void SetAll(Tensor<Vector<T, Size>, Dim>& data)
        {
            data_ = data;
        }

        inline void UpdateTime(const T time)
        {
            time_ = time;
        }

    protected:
        Tensor<Vector<T, Size>, Dim> data_;
        const unsigned nx_;
        const unsigned ny_;
        T time_;
};

// poisiton field
template <typename T, unsigned Dim = 2>
class Position : public Field<T, Dim, Dim>
{
    public:
        using vec = LCS::Vector<T, Dim>;

        // constructor
        Position(unsigned nx, unsigned ny): Field<T, Dim, Dim>(nx, ny) {}

        // setter for all values in the field
        void SetAll(const std::vector<T>& xrange, const std::vector<T>& yrange)
        {
            // make sure sizes match
            if (xrange.size()!=this->nx_||yrange.size()!=this->ny_)
                throw std::domain_error("sizes do not match!");

            // fill in the values
            for (unsigned i = 0; i < this->nx_; ++i)
                for (unsigned j = 0; j < this->ny_; ++j)
                    this->data_(i,j) = vec(xrange[i], yrange[j]);
        }

        // overload, use end points to set values
        void SetAll(const T& xmin, const T& xmax, const T& ymin, const T& ymax)
        {
            std::vector<T> xrange(this->nx_, 0), yrange(this->ny_, 0);
            int i = 0, j = 0;
            // fill in the values with uniform step
            std::generate(xrange.begin(), xrange.end(),
                    [&]{ return xmin + (i++) * (xmax-xmin)/(this->nx_-1);});
            std::generate(yrange.begin(), yrange.end(),
                    [&]{ return ymin + (j++) * (ymax-ymin)/(this->ny_-1);});
            SetAll(xrange, yrange);
        }

        // call function in base class
        void SetAll(Tensor<vec, Dim>& data)
        {
            Field<T, Dim, Dim>::SetAll(data);
        }

        inline void ReadFromFile(const std::string& file_name)
        {
            std::ifstream file;
            file.open(file_name, std::ios::in);
            if (!file.is_open())
                throw std::runtime_error("file does not open correctly!");

            // check if sizes match
            unsigned nx, ny;
            file >> nx; file >> ny;
            if (this->nx_!=nx || this->ny_!=ny)
                throw std::domain_error("sizes do not match!");

            // read in time stamp
            file >> this->time_;

            // read in data
            file >> this->data_;

            file.close();
        }

        // getter
        inline auto Get(const unsigned i, const unsigned j) const
        {
            return std::make_tuple(this->data_(i,j).x, this->data_(i,j).y);
        }

        inline void WriteToFile(const std::string& file_name) const
        {
            std::ofstream file;
            file.open(file_name);
            if (!file.is_open())
                throw std::runtime_error("file does not open correctly!");

            file.clear();
            file << this->nx_ << std::endl;
            file << this->ny_ << std::endl;
            file << this->time_ << std::endl;
            file << this->data_;
            file.close();
        }
        
        // update the position using the velocity field
        void Update(Velocity<T, Dim>& vel, T delta)
        {
            for (unsigned i = 0; i < this->nx_; ++i)
            {
                for (unsigned j = 0; j < this->ny_; ++j)
                {
                    T vx, vy;
                    std::tie(vx, vy) = vel.Get(i,j);
                    this->data_(i,j).x += vx * delta;
                    this->data_(i,j).y += vy * delta;
                }
            }
        }
};

// velocity field
template <typename T, unsigned Dim = 2>
class Velocity : public Field<T, Dim, Dim>
{
    public:
        using vec = LCS::Vector<T, Dim>;

        // constructor
        Velocity(unsigned nx, unsigned ny, Position<T, Dim>& pos):
            Field<T, Dim, Dim>(nx, ny), pos_(pos) {}

        // getter
        inline auto Get(const unsigned i, const unsigned j) const
        {
            return std::make_tuple(this->data_(i,j).x, this->data_(i,j).y);
        }

        template <typename A, typename B, unsigned C>
        friend class AnalyticVelocity;

    private:
        Position<T, Dim>& pos_; // position field corresponding to this velocity field
};

// velocity field with analytic velocity function
template <typename T, typename Func, unsigned Dim>
class AnalyticVelocity : public Velocity<T, Dim>
{
    public:
        using vec = LCS::Vector<T, Dim>;

        // constructor without parameters
        AnalyticVelocity(unsigned nx, unsigned ny, Position<T, Dim>& pos):
            Velocity<T, Dim>(nx, ny, pos), f_()
        {
            SetAll();
        }

        // constructor with parameters
        AnalyticVelocity(unsigned nx, unsigned ny, Position<T, Dim>& pos,
                std::vector<T>& parameters): Velocity<T, Dim>(nx, ny, pos), f_(parameters)
        {
            SetAll();
        }

        // use analytic function to set all velocity values
        void SetAll()
        {
            for (unsigned i = 0; i < this->nx_; ++i)
            {
                for (unsigned j = 0; j < this->ny_; ++j)
                {
                    T x, y, vx, vy;
                    std::tie(x, y) = this->pos_.Get(i, j);
                    std::tie(vx, vy) = f_(x, y, this->time_);
                    this->data_(i,j) = vec(vx, vy);
                    auto a = this->Get(0, 0);
                }
            }
        }

        inline Func& Function()
        {
            return f_;
        }

    private:
        Func f_;
};

}
