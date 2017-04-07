#pragma once

#include "basic.hpp"

namespace LCS {

template <typename T, unsigned Dim, unsigned Size>
struct FieldPolicy;

template <typename T, unsigned Dim>
class Position;

template <typename T, unsigned Dim>
class Velocity;

template <typename T, typename Func, unsigned Dim>
class ContinuousVelocity;


// general field
template <typename T, unsigned Dim = 2, unsigned Size = 2>
class Field
{
    public:
        Field(unsigned nx, unsigned ny): nx_(nx), ny_(ny), data_(nx, ny), time_() {}

        // getter
        inline auto& GetAll() const
        {
            return data_;
        }

        inline auto GetTime() const
        {
            return time_;
        }

        inline void ReadFromFile(const std::string& file_name)
        {
            FieldPolicy<T, Dim, Size>::ReadFromFile(*this, file_name);
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

        inline void WriteToFile(const std::string& file_name) const
        {
            FieldPolicy<T, Dim, Size>::WriteToFile(*this, file_name);
        }

        friend struct FieldPolicy<T, Dim, Size>;

    protected:
        const unsigned nx_;
        const unsigned ny_;
        Tensor<Vector<T, Size>, Dim> data_;
        T time_;
};


template <typename T, unsigned Dim, unsigned Size>
struct FieldPolicy {};

template <typename T, unsigned Size>
struct FieldPolicy<T, 2, Size>
{
    // write Field data to a file
    static inline void WriteToFile(const Field<T, 2, Size>& field, const std::string& file_name)
    {
        std::ofstream file;
        file.open(file_name);
        if (!file.is_open())
            throw std::runtime_error("file does not open correctly!");
    
        file.clear();
        file << field.nx_ << std::endl;
        file << field.ny_ << std::endl;
        file << field.time_ << std::endl;
        file << field.data_;
        file.close();
    }
    
    // read data form a file to Field
    static inline void ReadFromFile(Field<T, 2, Size>& field, const std::string& file_name)
    {
        std::ifstream file;
        file.open(file_name, std::ios::in);
        if (!file.is_open())
            throw std::runtime_error("file does not open correctly!");

        // check if sizes match
        unsigned nx, ny;
        file >> nx; file >> ny;
        if (field.nx_!=nx || field.ny_!=ny)
            throw std::domain_error("sizes do not match!");

        // read in time stamp
        file >> field.time_;

        // read in data
        file >> field.data_;

        file.close();
    }

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

            pos_xrange_ = xrange;
            pos_yrange_ = yrange;
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

        // getter
        inline auto Get(const unsigned i, const unsigned j) const
        {
            return std::make_tuple(this->data_(i,j).x, this->data_(i,j).y);
        }

        inline auto& GetRange(const unsigned axis)
        {
            // need to check the input and if ranges exist
            if (axis == 0)
                return pos_xrange_;
            else
                return pos_yrange_;
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

                    // check if the position is out of bound
                    // need to deal with NAN values
                    if (out_of_bound_ != nullptr)
                        if (this->data_(i,j).x < bound_xmin_ ||
                                this->data_(i,j).x > bound_xmax_ ||
                                this->data_(i,j).y < bound_ymin_ ||
                                this->data_(i,j).y > bound_ymax_)
                            out_of_bound_->SetValue(i, j, true);
                }
            }
        }

        // initialize out of bound tensor
        inline void InitializeOutOfBoundTensor()
        {
            // default value is false
            out_of_bound_.reset(new Tensor<bool, Dim>(this->nx_, this->ny_));
        }

        inline bool IsOutOfBound(const unsigned i, const unsigned j)
        {
            if (out_of_bound_ != nullptr)
                return out_of_bound_->GetValue(i,j);
            else
                return false; // default value
        }
        // set the out boundary 
        inline void SetBound(const T& xmin, const T& xmax, const T& ymin, const T& ymax)
        {
            bound_xmin_ = xmin;
            bound_xmax_ = xmax;
            bound_ymin_ = ymin;
            bound_ymax_ = ymax;
        }

    private:
        std::unique_ptr<Tensor<bool, Dim>> out_of_bound_;
        std::vector<T> pos_xrange_;
        std::vector<T> pos_yrange_;

        T bound_xmin_;
        T bound_xmax_;
        T bound_ymin_;
        T bound_ymax_;
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

        inline auto& GetPosition()
        {
            return pos_;
        }

        // interploate from another velocity field
        void InterpolateFrom(Velocity<T, Dim>& ref_vel)
        {
            // interpolation only works for orthogonal coordinates
            auto ref_pos_x = ref_vel.GetPosition().GetRange(0);
            auto ref_pos_y = ref_vel.GetPosition().GetRange(1);

            T pos_x, pos_y;
            Tensor<Vector<T, 2>, 2> ref_v(2,2);

            for (unsigned i = 0; i < this->nx_; ++i)
            {
                for (unsigned j = 0; j < this->ny_; ++j)
                {
                    // skip if the current position is out of bound
                    if (pos_.IsOutOfBound(i,j)) continue;

                    std::tie(pos_x, pos_y) = pos_.Get(i, j);

                    auto iter_x = std::upper_bound(ref_pos_x.begin(), ref_pos_x.end(), pos_x);
                    if (iter_x == ref_pos_x.end()) --iter_x;
                    if (iter_x == ref_pos_x.begin()) ++iter_x;
                    int i_next = std::distance(ref_pos_x.begin(), iter_x);
                    int i_pre = i_next - 1;

                    auto iter_y = std::upper_bound(ref_pos_y.begin(), ref_pos_y.end(), pos_y);
                    if (iter_y == ref_pos_y.end()) --iter_y;
                    if (iter_y == ref_pos_y.begin()) ++iter_y;
                    int j_next = std::distance(ref_pos_y.begin(), iter_y);
                    int j_pre = j_next - 1;

                    // add another getter for Tensor class?
                    std::tie(ref_v(0,0).x, ref_v(0,0).y) = ref_vel.Get(i_pre, j_pre);
                    std::tie(ref_v(0,1).x, ref_v(0,1).y) = ref_vel.Get(i_pre, j_next);
                    std::tie(ref_v(1,0).x, ref_v(1,0).y) = ref_vel.Get(i_next, j_pre);
                    std::tie(ref_v(1,1).x, ref_v(1,1).y) = ref_vel.Get(i_next, j_next);
                    
                    auto vx1 = interpolate(ref_pos_x[i_pre], ref_pos_x[i_next],
                            ref_v(0,0).x, ref_v(1,0).x, pos_x);
                    auto vx2 = interpolate(ref_pos_x[i_pre], ref_pos_x[i_next],
                            ref_v(0,1).x, ref_v(1,1).x, pos_x);
                    auto vxm = interpolate(ref_pos_y[j_pre], ref_pos_y[j_next],
                            vx1, vx2, pos_y);

                    auto vy1 = interpolate(ref_pos_x[i_pre], ref_pos_x[i_next],
                            ref_v(0,0).y, ref_v(1,0).y, pos_x);
                    auto vy2 = interpolate(ref_pos_x[i_pre], ref_pos_x[i_next],
                            ref_v(0,1).y, ref_v(1,1).y, pos_x);
                    auto vym = interpolate(ref_pos_y[j_pre], ref_pos_y[j_next],
                            vy1, vy2, pos_y);

                    this->data_(i,j).x = vxm;
                    this->data_(i,j).y = vym;
                }
            }
        }

    protected:
        Position<T, Dim>& pos_; // position field corresponding to this velocity field
};

// velocity field with continuous velocity function
template <typename T, typename Func, unsigned Dim>
class ContinuousVelocity : public Velocity<T, Dim>
{
    public:
        using vec = LCS::Vector<T, Dim>;

        // constructor without parameters
        ContinuousVelocity(unsigned nx, unsigned ny, Position<T, Dim>& pos, T time=0):
            Velocity<T, Dim>(nx, ny, pos), f_()
        {
            this->time_ = time;
            SetAll();
        }

        // constructor with parameters
        ContinuousVelocity(unsigned nx, unsigned ny, Position<T, Dim>& pos,
                std::vector<T>& parameters, T time=0): Velocity<T, Dim>(nx, ny, pos), f_(parameters)
        {
            this->time_ = time;
            SetAll();
        }

        // use continuous function to set all velocity values
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
