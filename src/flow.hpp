/**
    @file flow.hpp
    @brief Classes for flow fields.

    This file contains classes for reprensenting flow fields. It includes a general FlowField class, and its subclasses DiscreteFlowField and ContinuousFlowField, which are for discrete input data (e.g. from experiments and simulations) and continuous input data (when analytic velocity function is known) respectively.
    */
#pragma once

#include "field.hpp"
#include "basic.hpp"

namespace LCS {

template <typename T, unsigned Dim>
class FlowField;

template <typename T, typename Func, unsigned Dim>
class ContinuousFlowField;


/** @enum Direction
    Enum for particle advection direction.
    */
enum Direction
{
    Forward, /**< Forward in time. */
    Backward /**< Backward in time.*/
};

/** @brief Class for flow field.

    This class is used for representing flow fields, which contains many information (such as positions and velocities for all fluid particles) about the flow field.
    @tparam T Numeric data type of all the numeric values in the flow.
    @tparam Dim Dimension of the flow field (2 or 3).
    */
template <typename T, unsigned Dim = 2>
class FlowField
{
    public:
        /** Constructor for initializating the flow field.
            @param nx The number of grid points in \f$x\f$-direction.
            @param ny The number of grid points in \f$y\f$-direction.
            */
        FlowField(unsigned nx, unsigned ny):
            nx_(nx), ny_(ny), delta_(), initial_time_(), current_time_(), step_(), 
            initial_pos_(new Position<T,Dim>(nx,ny)),
            current_pos_(new Position<T,Dim>(nx,ny)) {}

        /** Get the initial positions of all particles.
            @return Position field associated with the initial time of the flow.
            */
        inline auto& InitialPosition()
        {
            return *initial_pos_;
        }

        /** Get the current positions of all particles.
            @return Position field associated with the current time of the flow.
            */
        inline auto& CurrentPosition()
        {
            return *current_pos_;
        }

        /** Get the current velocities of all particles.
            @return Velocity field associated with the current time of the flow.
            */
        inline virtual Velocity<T, Dim>& CurrentVelocity()
        {
            if (current_vel_ == nullptr)
                throw std::invalid_argument("current velocity not set!");

            return *current_vel_;
        }

        /** Get the current time of the flow.
            @return Current time of the flow.
            */
        inline auto GetTime()
        {
            return current_time_;
        }

        /** Get the direction associated with the flow field.
            @return Direction that shows that the flow advection is forward or backward.
            */
        inline auto GetDirection()
        {
            return direction_;
        }

        /** Copy the initial position to the current position. It is used before particle advection. The actual implementaions are in the subclasses.
            */
        virtual void CopyInitialPositionToCurrentPosition() {};

        /** Set the time step for particle advection.
            @param delta Time step for advection.
            */
        inline void SetDelta(const T delta)
        {
            assert(delta > 0);
            delta_ = delta;
        }

        /** Set the number of steps for particle advection.
            @param step Number of total advection steps.
            */
        inline void SetStep(const unsigned step)
        {
            step_ = step;
        }

        /** Set the current velocity. The actual implementations are in the subclasses.
            */
        virtual void SetCurrentVelocity()
        {}

        /** Set the advection direction. The actual implementations are in the subclesses.
            @param direction Direction for particle advection.
            */
        virtual void SetDirection(const Direction direction)
        {}

        /** Set the initial time.
            @param time Initial time of calculation.
            */
        inline void SetInitialTime(const T time)
        {
            initial_time_ = time;
            initial_pos_->UpdateTime(time);
            UpdateTime(time);
        }

        /** Update time of the flow field using the time step.
            */ 
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

        /** Update time of the flow field using the provided time.
            @param time New time for the flow field.
            */
        inline void UpdateTime(const T time)
        {
            current_time_ = time;
            current_pos_->UpdateTime(current_time_);
            if (current_vel_ != nullptr) current_vel_->UpdateTime(current_time_);
        }

        /** Particle advection (calculating particle trajectories).  */
        void Run()
        {
            switch(direction_)
            {
                case Forward:
                    std::cout << "Particle forward advection begins" << std::endl;
                    break;
                case Backward:
                    std::cout << "Particle backward advection begins" << std::endl;
                    break;
                default: break;
            }
            CopyInitialPositionToCurrentPosition();

            Clock clock;
            for (unsigned i = 0; i < step_; ++i)
            {
                clock.Begin();
                std::stringstream ss;
                ss << "[" << std::setw(std::to_string(step_).length()) << (i+1)
                    << "/" << step_ << "]" <<
                    " Simulation time: " << current_time_;

                SetCurrentVelocity();

                switch(direction_)
                {
                    case Forward: current_pos_->Update(*current_vel_, delta_); break;
                    case Backward: current_pos_->Update(*current_vel_, -delta_); break;
                    default: break;
                }

                UpdateTime();

                ss << " - " << current_time_;
                clock.End();
                ss << " (Execution time: " <<
                    std::setprecision(4) << clock.GetTotalElapsedTime() <<
                    "s)" << std::endl;
                std::cout << ss.str();
            }
            switch(direction_)
            {
                case Forward:
                    std::cout << "Particle forward advection ends" << std::endl;
                    break;
                case Backward:
                    std::cout << "Particle backward advection ends" << std::endl;
                    break;
                default: break;
            }
        }

    protected:
        const unsigned nx_;/**<The number of grid points in \f$x\f$-direction.*/
        const unsigned ny_;/**<The number of grid points in \f$y\f$-direction.*/
        T delta_; /**<Time step for integration.*/
        T initial_time_;/**<Initial time of the calculation.*/
        T current_time_;/**<Current time of the calculation.*/
        unsigned step_;/**<Number of calculation steps.*/
        Direction direction_ = Forward;/**<Direction of the calculation.*/
        std::unique_ptr<Position<T, Dim>> initial_pos_;/**<Pointer to the Position field at the initial time.*/
        std::unique_ptr<Position<T, Dim>> current_pos_;/**<Pointer to the Position field at the current time.*/
        std::shared_ptr<Velocity<T, Dim>> current_vel_;/**<Pointer to the Velocity field at the current time.*/
};

/** @brief Class for flow fields with discrete data.

    This class is used for representing flow fields with discrete data which typically come from experiments or CFD simulations.
    @tparam T Numeric data type of all the numeric values in the flow.
    @tparam Dim Dimension of the flow field (2 or 3).
    */
template <typename T, unsigned Dim = 2>
class DiscreteFlowField : public FlowField<T, Dim>
{
    public:
        /** Constructor for initializating the discrete flow field.
            @param nx The number of grid points in \f$x\f$-direction for calculation.
            @param ny The number of grid points in \f$y\f$-direction for calculation.
            @param data_nx The number of grid points in \f$x\f$-direction in the input flow data.
            @param data_ny The number of grid points in \f$y\f$-direction in the input flow data.
            */
        DiscreteFlowField(unsigned nx, unsigned ny, unsigned data_nx, unsigned data_ny):
            FlowField<T, Dim>(nx, ny), data_nx_(data_nx), data_ny_(data_ny),
            data_delta_(), current_data_time_(), begin_data_time_(), end_data_time_(),
            vel_file_name_prefix_(""), vel_file_name_suffix_(".txt"),
            data_pos_(new Position<T,Dim>(data_nx, data_ny)),
            previous_data_vel_(new Velocity<T,Dim>(data_nx, data_ny, *data_pos_)),
            next_data_vel_(new Velocity<T,Dim>(data_nx, data_ny, *data_pos_)),
            current_data_vel_(new Velocity<T,Dim>(data_nx, data_ny, *data_pos_)) {}

        /** Constructor for initializating the discrete flow field when the input data and calculation data have the same resolution.
            @param nx The number of grid points in \f$x\f$-direction.
            @param ny The number of grid points in \f$y\f$-direction.
            */
        DiscreteFlowField(unsigned nx, unsigned ny):
            DiscreteFlowField(nx, ny, nx, ny) {}

        /** Set the prefix of the file names of the input velocity data.
            @param prefix A string contains the prefix of file names.
            */
        inline void SetVelocityFileNamePrefix(const std::string prefix)
        {
            vel_file_name_prefix_ = prefix;
        }

        /** Set the suffix of the file names of the input velocity data.
            @param suffix A string contains the prefix of file names.
            */
        inline void SetVelocityFileNameSuffix(const std::string suffix)
        {
            vel_file_name_suffix_ = suffix;
        }

        /** Read velocity data from a file.
            @param data_vel Velocity fields that contains the input velocity data from a file.
            */
        inline void ReadDataVelocityFromFile(Velocity<T, Dim>& data_vel)
        {
            std::string file_name = vel_file_name_prefix_ +
                std::to_string(static_cast<int>(data_vel.GetTime())) + 
                vel_file_name_suffix_;

            data_vel.ReadFromFile(file_name);

            std::cout << "Read velocity data at time = " <<
                data_vel.GetTime() << " from " << file_name << std::endl;
        }

        /** Calculate and set the current velocities of the flow particles from temporal and spatical interpolations of the velocity data.
            */
        inline void SetCurrentVelocity()
        {
            T signed_data_delta = (this->direction_ == Forward) ? data_delta_ : (-data_delta_);

            // initialize data velocities
            if (this->current_time_ == this->initial_time_)
            {
                current_data_time_ = begin_data_time_;

                // find the data velocity field that is temporally closest to the initial time
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

                // read two temporally adjacent data field (for later interpolation)
                #pragma omp parallel
                #pragma omp single nowait
                {
                    #pragma omp task
                    {
                        previous_data_vel_->UpdateTime(current_data_time_);
                        ReadDataVelocityFromFile(*previous_data_vel_);
                    }
                    #pragma omp task
                    {
                        next_data_vel_->UpdateTime(current_data_time_ + signed_data_delta);
                        ReadDataVelocityFromFile(*next_data_vel_);
                    }
                }

            // check if the current time is outside the time
            // range of two temporally adjacent data field
            } else if (((this->current_time_ >= current_data_time_ + signed_data_delta) &&
                    (end_data_time_ > current_data_time_ + signed_data_delta) &&
                    this->direction_ == Forward) ||
                    ((this->current_time_ <= current_data_time_ + signed_data_delta) &&
                    (end_data_time_ < current_data_time_ + signed_data_delta) &&
                    this->direction_ == Backward))
            {
                // update data velocities

                current_data_time_ += signed_data_delta;

                // read the next two temporally adjacent data velocity field so that
                // the current time could be within the time range of the two data field
                #pragma omp parallel
                #pragma omp single nowait
                {
                    #pragma omp task
                    {
                        previous_data_vel_->UpdateTime(current_data_time_);
                        ReadDataVelocityFromFile(*previous_data_vel_);
                    }
                    #pragma omp task
                    {
                        next_data_vel_->UpdateTime(current_data_time_ + signed_data_delta);
                        ReadDataVelocityFromFile(*next_data_vel_);
                    }
                }
            }
            
            // temporal interpolation
            interpolate(previous_data_vel_->GetTime(),
                    next_data_vel_->GetTime(),
                    *previous_data_vel_, *next_data_vel_,
                    this->current_time_, *current_data_vel_);

            // spatial interpolation
            this->current_vel_->InterpolateFrom(*current_data_vel_);

        }

        /** Set the advection direction.
            @param direction Direction for particle advection.
            */
        inline void SetDirection(const Direction direction)
        {
            this->direction_ = direction;
            SetDataTimeRange(begin_data_time_, end_data_time_);
        }

        /** Copy the initial position to the current position.
            */
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

        /** Get the Position field of the input data.
            */
        inline auto& DataPosition()
        {
            return *data_pos_;
        }

        /** Get the current velocity field of the input data.
            */
        inline auto& CurrentDataVelocity()
        {
            return *previous_data_vel_;
        }

        /** Set the time difference between two adjacent data files.
            @param delta Time difference between two adjacent data files.
            */
        inline void SetDataDelta(const T delta)
        {
            data_delta_ = delta;
        }

        /** Set the time range of the data files.
            @param t1 One end point of the data time range.
            @param t2 The another end point of the data time range.
            */
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
        const unsigned data_nx_; /**<The number of grid points in \f$x\f$-direction in the input flow data.*/
        const unsigned data_ny_; /**<The number of grid points in \f$y\f$-direction in the input flow data.*/
        T data_delta_; /**<Time difference between two adjacent data files.*/
        T current_data_time_;/**<Time of the currently used velocity data file.*/
        T begin_data_time_;/**<Time of the data file that is used at the beginning of calculation.*/
        T end_data_time_;/**<Time of the data file that is used at the end of the calculation.*/
        std::string vel_file_name_prefix_;/**<Prefix of the data file names.*/
        std::string vel_file_name_suffix_;/**<Suffix of the data file names.*/
        std::unique_ptr<Position<T, Dim>> data_pos_;/**<Pointer to the Position field of the input data.*/
        std::unique_ptr<Velocity<T, Dim>> previous_data_vel_;/**<Pointer to the Velocity field of the input data whose time is earlier than the curren time.*/
        std::unique_ptr<Velocity<T, Dim>> next_data_vel_;/**<Pointer to the Velocity field of the input data whose time is later than the current time.*/
        std::unique_ptr<Velocity<T, Dim>> current_data_vel_; /**<Pointer to the Velocity filed that is the temporal interpolation of two adjenct data files.*/ 
};


/** @brief Class for flow fields with a continous velocity function.

    This class is used for representing flow fields with a given continous velocity function, which means that we could obtain the exact velocities of any point in the field.
    @tparam T Numeric data type of all the numeric values in the flow.
    @tparam Func Continous velocity function for this field.
    @tparam Dim Dimension of the flow field (2 or 3).
    */
template <typename T, typename Func, unsigned Dim = 2>
class ContinuousFlowField : public FlowField<T, Dim>
{
    public:
        /** Constructor for initializating the flow field.
            @param nx The number of grid points in \f$x\f$-direction.
            @param ny The number of grid points in \f$y\f$-direction.
            */
        ContinuousFlowField(unsigned nx, unsigned ny):
            FlowField<T, Dim>(nx, ny), parameters_() {}

        /** Constructor for initializating the flow field that the velocity function has parameters.
            @param nx The number of grid points in \f$x\f$-direction.
            @param ny The number of grid points in \f$y\f$-direction.
            @param parameters A vector that contains all the parameters of the velocity function.
            */
        ContinuousFlowField(unsigned nx, unsigned ny, std::vector<T>& parameters):
            FlowField<T, Dim>(nx, ny), parameters_(parameters) {}

        /** Calculate and set the current velocities of the flow particles using the velocity function.
            */
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

        /** Get the velocity field at the current time
            @return ContinuousVelocity field that contains the current velocities.
            */
        inline ContinuousVelocity<T, Func, Dim>& CurrentVelocity()
        {
            if (current_continuous_vel_ == nullptr)
                throw std::invalid_argument("current velocity not set!");

            return *current_continuous_vel_;
        }

        /** Copy the initial position to the current position.
            */
        inline void CopyInitialPositionToCurrentPosition()
        {
            auto initial_pos_data = this->initial_pos_->GetAll();
            this->current_pos_->SetAll(initial_pos_data);
            this->current_pos_->GetRange(0) = this->initial_pos_->GetRange(0);
            this->current_pos_->GetRange(1) = this->initial_pos_->GetRange(1);
            this->current_pos_->UpdateTime(this->initial_pos_->GetTime());
        }

        /** Set the advection direction.
            @param direction Direction for particle advection.
            */
        inline void SetDirection(const Direction direction)
        {
            this->direction_ = direction;
        }

    private:
        std::shared_ptr<ContinuousVelocity<T, Func, Dim>> current_continuous_vel_;/**<Pointer to the ContinuousVelocity field at the current time.*/
        std::vector<T> parameters_;/**<A vector of parameters for the velocity function associated with this flow field.*/
};

}
