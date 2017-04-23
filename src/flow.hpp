#pragma once

#include "field.hpp"
#include "basic.hpp"

namespace LCS {

template <typename T, unsigned Dim>
class FlowField;

template <typename T, typename Func, unsigned Dim>
class ContinuousFlowField;


// calculation direction
enum Direction
{
    Forward,
    Backward
};

// flow field
template <typename T, unsigned Dim = 2>
class FlowField
{
    public:
        // constructor
        FlowField(unsigned nx, unsigned ny):
            nx_(nx), ny_(ny), delta_(), initial_time_(), current_time_(), step_(), 
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

        virtual void CopyInitialPositionToCurrentPosition() {};

        inline void SetDelta(const T delta)
        {
            assert(delta > 0);
            delta_ = delta;
        }

        inline void SetStep(const unsigned step)
        {
            step_ = step;
        }

        // overloaded in ContinuousFlowField class
        virtual void SetCurrentVelocity()
        {}

        // set calculation direction
        virtual void SetDirection(const Direction direction)
        {}

        // set initial time
        inline void SetInitialTime(const T time)
        {
            initial_time_ = time;
            initial_pos_->UpdateTime(time);
            UpdateTime(time);
        }

        // update time
        inline void UpdateTime()
        {
            switch(direction_)
            {
                case Forward: current_time_ += delta_; break;
                case Backward: current_time_ -= delta_; break;
                default: break;
            }
            current_pos_->UpdateTime(current_time_);
            current_vel_->UpdateTime(current_time_);
        }

        inline void UpdateTime(const T time)
        {
            current_time_ = time;
            current_pos_->UpdateTime(current_time_);
            if (current_vel_ != nullptr) current_vel_->UpdateTime(current_time_);
        }

        // calculate the trajectories
        void Run()
        {
            std::cout << "Particle advection begins" << std::endl;
            CopyInitialPositionToCurrentPosition();

            Clock clock;
            for (unsigned i = 0; i < step_; ++i)
            {
                clock.Begin();
                std::cout << "Step " << i <<
                    " (time = " << current_time_ << ") begins" << std::endl;

                SetCurrentVelocity();

                switch(direction_)
                {
                    case Forward: current_pos_->Update(*current_vel_, delta_); break;
                    case Backward: current_pos_->Update(*current_vel_, -delta_); break;
                    default: break;
                }

                UpdateTime();

                std::cout << "Step " << i <<
                    " (time = " << current_time_ << ") ends" << std::endl;
                clock.End();
                std::cout << "Step executing time: " << clock.GetElapsedTime() <<
                    "s  Total executing time: " << clock.GetTotalElapsedTime() <<
                    "s" << std::endl;
            }
            std::cout << "Particle advection ends" << std::endl;
        }

    protected:
        const unsigned nx_;
        const unsigned ny_;
        T delta_; // time step for integration
        T initial_time_;
        T current_time_;
        unsigned step_;
        Direction direction_ = Forward;
        std::unique_ptr<Position<T, Dim>> initial_pos_;
        std::unique_ptr<Position<T, Dim>> current_pos_;
        std::shared_ptr<Velocity<T, Dim>> current_vel_;
};

template <typename T, unsigned Dim = 2>
class DiscreteFlowField : public FlowField<T, Dim>
{
    public:
        // constructor
        DiscreteFlowField(unsigned nx, unsigned ny, unsigned data_nx, unsigned data_ny):
            FlowField<T, Dim>(nx, ny), data_nx_(data_nx), data_ny_(data_ny),
            data_delta_(), current_data_time_(), begin_data_time_(), end_data_time_(),
            vel_file_name_prefix_(""), vel_file_name_suffix_(".txt"),
            data_pos_(new Position<T,Dim>(data_nx, data_ny)),
            previous_data_vel_(new Velocity<T,Dim>(data_nx, data_ny, *data_pos_)),
            next_data_vel_(new Velocity<T,Dim>(data_nx, data_ny, *data_pos_)),
            current_data_vel_(new Velocity<T,Dim>(data_nx, data_ny, *data_pos_)) {}

        // datainal data and calculation data have the same resolution
        DiscreteFlowField(unsigned nx, unsigned ny):
            DiscreteFlowField(nx, ny, nx, ny) {}

        inline void SetVelocityFileNamePrefix(const std::string prefix)
        {
            vel_file_name_prefix_ = prefix;
        }

        inline void SetVelocityFileNameSuffix(const std::string suffix)
        {
            vel_file_name_suffix_ = suffix;
        }

        inline void ReadDataVelocityFromFile(Velocity<T, Dim>& data_vel)
        {
            std::string file_name = vel_file_name_prefix_ +
                std::to_string(static_cast<int>(data_vel.GetTime())) + 
                vel_file_name_suffix_;

            data_vel.ReadFromFile(file_name);

            std::cout << "Read velocity data at time = " <<
                data_vel.GetTime() << " from " << file_name << std::endl;
        }

        inline void SetCurrentVelocity()
        {
            T signed_data_delta = (this->direction_ == Forward) ? data_delta_ : (-data_delta_);

            // initialize data velocities
            if (this->current_time_ == this->initial_time_)
            {
                current_data_time_ = begin_data_time_;

                switch(this->direction_)
                {
                    case Forward:
                        while (current_data_time_ < this->initial_time_)
                            current_data_time_ += signed_data_delta;
                        break;
                    case Backward:
                        while (current_data_time_ > this->initial_time_)
                            current_data_time_ -= signed_data_delta;
                        break;
                    default: break;
                }

                #pragma omp parallel
                #pragma omp single nowait
                {
                    #pragma omp task
                    previous_data_vel_->UpdateTime(current_data_time_);
                    ReadDataVelocityFromFile(*previous_data_vel_);
                    #pragma omp task
                    next_data_vel_->UpdateTime(current_data_time_ + signed_data_delta);
                    ReadDataVelocityFromFile(*next_data_vel_);
                }

            } else if (((this->current_time_ >= current_data_time_ + signed_data_delta) &&
                    (end_data_time_ > current_data_time_ + signed_data_delta) &&
                    this->direction_ == Forward) ||
                    ((this->current_time_ <= current_data_time_ + signed_data_delta) &&
                    (end_data_time_ < current_data_time_ + signed_data_delta) &&
                    this->direction_ == Backward))
            {
                // update data velocities
                current_data_time_ += signed_data_delta;
                #pragma omp parallel
                #pragma omp single nowait
                {
                    #pragma omp task
                    previous_data_vel_->UpdateTime(current_data_time_);
                    ReadDataVelocityFromFile(*previous_data_vel_);

                    #pragma omp task
                    next_data_vel_->UpdateTime(current_data_time_ + signed_data_delta);
                    ReadDataVelocityFromFile(*next_data_vel_);
                }
            }
            
            interpolate(previous_data_vel_->GetTime(),
                    next_data_vel_->GetTime(),
                    *previous_data_vel_, *next_data_vel_,
                    this->current_time_, *current_data_vel_);

            this->current_vel_->InterpolateFrom(*current_data_vel_);

        }

        inline void SetDirection(const Direction direction)
        {
            this->direction_ = direction;
            SetDataTimeRange(begin_data_time_, end_data_time_);
        }

        inline void CopyInitialPositionToCurrentPosition()
        {
            auto initial_pos_data = this->initial_pos_->GetAll();
            this->current_pos_->SetAll(initial_pos_data);
            this->current_pos_->GetRange(0) = this->initial_pos_->GetRange(0);
            this->current_pos_->GetRange(1) = this->initial_pos_->GetRange(1);
            this->current_pos_->UpdateTime(this->initial_pos_->GetTime());
            this->current_pos_->InitializeOutOfBoundTensor();


            // set out boundary
            T xmin, xmax, ymin, ymax;
            std::tie(xmin, ymin) = data_pos_->Get(0,0);
            std::tie(xmax, ymax) = data_pos_->Get(data_nx_-1,data_ny_-1);
            this->current_pos_->SetBound(xmin, xmax, ymin, ymax);

            // set corresponding velocity tensor
            this->current_vel_.reset(new Velocity<T, Dim>
                    (this->nx_, this->ny_, *(this->current_pos_)));
        }

        inline auto& DataPosition()
        {
            return *data_pos_;
        }

        inline auto& CurrentDataVelocity()
        {
            return *previous_data_vel_;
        }

        inline void SetDataDelta(const T delta)
        {
            data_delta_ = delta;
        }

        inline void SetDataTimeRange(const T t1, const T t2)
        {
            T begin_time = (t2 >= t1) ? t1 : t2;
            T end_time = (t2 >= t1) ? t2 : t2;

            switch(this->direction_)
            {
                case Forward:
                    begin_data_time_ = begin_time;
                    end_data_time_ = end_time;
                    break;
                case Backward:
                    begin_data_time_ = end_time;
                    end_data_time_ = begin_time;
                    break;
                default: break;
            }
        }

    private:
        const unsigned data_nx_;
        const unsigned data_ny_;
        T data_delta_; // time difference between two adjacent data files
        T current_data_time_;
        T begin_data_time_;
        T end_data_time_;
        std::string vel_file_name_prefix_;
        std::string vel_file_name_suffix_;
        std::unique_ptr<Position<T, Dim>> data_pos_;
        std::unique_ptr<Velocity<T, Dim>> previous_data_vel_;
        std::unique_ptr<Velocity<T, Dim>> next_data_vel_;
        std::unique_ptr<Velocity<T, Dim>> current_data_vel_; // time interpolation of two adjenct data file 
};


template <typename T, typename Func, unsigned Dim = 2>
class ContinuousFlowField : public FlowField<T, Dim>
{
    public:
        // constructor
        ContinuousFlowField(unsigned nx, unsigned ny):
            FlowField<T, Dim>(nx, ny), parameters_() {}

        ContinuousFlowField(unsigned nx, unsigned ny, std::vector<T>& parameters):
            FlowField<T, Dim>(nx, ny), parameters_(parameters) {}

        inline void SetCurrentVelocity()
        {
            // no parameters
            if (parameters_.size() == 0)
                current_continuous_vel_.reset(new ContinuousVelocity<T, Func, Dim>
                    (this->nx_, this->ny_, *(this->current_pos_),this->current_time_));
            else
            // there are parameters for continuous function 
                current_continuous_vel_.reset(new ContinuousVelocity<T, Func, Dim>
                    (this->nx_, this->ny_, *(this->current_pos_),
                     parameters_, this->current_time_));
            
            current_continuous_vel_->UpdateTime(this->current_time_);
            this->current_vel_ = current_continuous_vel_;
        }

        inline ContinuousVelocity<T, Func, Dim>& CurrentVelocity()
        {
            if (current_continuous_vel_ == nullptr)
                throw std::invalid_argument("current velocity not set!");

            return *current_continuous_vel_;
        }

        inline void CopyInitialPositionToCurrentPosition()
        {
            auto initial_pos_data = this->initial_pos_->GetAll();
            this->current_pos_->SetAll(initial_pos_data);
            this->current_pos_->GetRange(0) = this->initial_pos_->GetRange(0);
            this->current_pos_->GetRange(1) = this->initial_pos_->GetRange(1);
            this->current_pos_->UpdateTime(this->initial_pos_->GetTime());
        }

        inline void SetDirection(const Direction direction)
        {
            this->direction_ = direction;
        }

    private:
        std::shared_ptr<ContinuousVelocity<T, Func, Dim>> current_continuous_vel_;
        std::vector<T> parameters_;
};

}
