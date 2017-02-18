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

        inline auto& CurrentVelocity()
        {
            return *current_vel_;
        }

        inline auto GetTime()
        {
            return current_time_;
        }

        inline void CopyInitialPositionToCurrentPosition()
        {
            current_pos_->SetAll(initial_pos_->GetAll());
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
        std::unique_ptr<Velocity<T, Dim>> current_vel_;
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
            FlowField<T, Dim>(nx, ny), f_() {}

        inline void SetCurrentVelocity()
        {
            this->current_vel_.reset(new AnalyticVelocity<T, Func, Dim>
                (this->nx_, this->ny_, *(this->current_pos_)));
            this->current_vel_->UpdateTime(this->current_time_);
        }

    private:
        Func f_;
};

// poisiton field
template <typename T, unsigned Dim = 2>
class Position
{
    public:
        using vec = LCS::Vector<T, Dim>;

        // constructor
        Position(unsigned nx, unsigned ny): nx_(nx), ny_(ny), time_(),
            data_(nx, ny) {
            data_(1,1) = vec(2,2);
            }

        // setter for all values in the field
        void SetAll(const std::vector<T>& xrange, const std::vector<T>& yrange)
        {
            // make sure dimensions match
            if (xrange.size()!=nx_||yrange.size()!=ny_)
                throw std::domain_error("dimensions do not match!");

            // fill in the values
            for (unsigned i = 0; i < nx_; ++i)
                for (unsigned j = 0; j < ny_; ++j)
                    data_(i,j) = vec(xrange[i], yrange[j]);
        }

        // overload, use end points to set values
        void SetAll(const T& xmin, const T& xmax, const T& ymin, const T& ymax)
        {
            std::vector<T> xrange(nx_, 0), yrange(ny_, 0);
            int i = 0, j = 0;
            // fill in the values with uniform step
            std::generate(xrange.begin(), xrange.end(),
                    [&]{ return xmin + (i++) * (xmax-xmin)/(nx_-1);});
            std::generate(yrange.begin(), yrange.end(),
                    [&]{ return ymin + (j++) * (ymax-ymin)/(ny_-1);});
            SetAll(xrange, yrange);
        }

        void SetAll(LCS::Tensor<vec, Dim> data)
        {
            data_ = data;
        }

        inline void UpdateTime(const T time)
        {
            time_ = time;
        }

        // getter
        inline auto Get(const unsigned i, const unsigned j) const
        {
            return std::make_tuple(data_(i,j).x, data_(i,j).y);
        }

        inline auto& GetAll() const
        {
            return data_;
        }

        inline auto GetTime() const
        {
            return time_;
        }
        
        // update the position using the velocity field
        void Update(Velocity<T, Dim>& vel, T delta)
        {
            for (unsigned i = 0; i < nx_; ++i)
            {
                for (unsigned j = 0; j < ny_; ++j)
                {
                    T vx, vy;
                    std::tie(vx, vy) = vel.Get(i,j);
                    data_(i,j).x += vx * delta;
                    data_(i,j).y += vy * delta;
                }
            }
        }

    private:
        //std::vector<std::vector<vec>> data_;
        LCS::Tensor<vec, Dim> data_;
        const unsigned nx_;
        const unsigned ny_;
        T time_;
};

// velocity field
template <typename T, unsigned Dim = 2>
class Velocity
{
    public:
        using vec = LCS::Vector<T, Dim>;

        // constructor
        Velocity(unsigned nx, unsigned ny, Position<T, Dim>& pos):
            nx_(nx), ny_(ny), time_(), data_(nx, ny), pos_(pos) {
                T x, y;
                std::tie(x, y) = pos_.Get(2,2);
                std::cout << x << std::endl;
            }

        // setter
        void SetAll(LCS::Tensor<vec, Dim> data)
        {
            data_ = data;
        }

        inline void UpdateTime(const T time)
        {
            time_ = time;
        }

        // getter
        inline auto Get(const unsigned i, const unsigned j) const
        {
            return std::make_tuple(data_(i,j).x, data_(i,j).y);
        }

        inline auto GetTime() const
        {
            return time_;
        }

        template <typename A, typename B, unsigned C>
        friend class AnalyticVelocity;

    private:
        LCS::Tensor<vec, Dim> data_;
        const unsigned nx_;
        const unsigned ny_;
        Position<T, Dim>& pos_; // position field corresponding to this velocity field
        T time_;
};

// velocity field with analytic velocity function
template <typename T, typename Func, unsigned Dim>
class AnalyticVelocity : public Velocity<T, Dim>
{
    public:
        using vec = LCS::Vector<T, Dim>;

        AnalyticVelocity(unsigned nx, unsigned ny, Position<T, Dim>& pos):
            Velocity<T, Dim>(nx, ny, pos), f_()
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

    private:
        Func f_;
};

}
