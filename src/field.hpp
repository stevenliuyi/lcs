#pragma once

#include "basic.hpp"
#include <omp.h>

namespace LCS {

template <typename T, unsigned Dim, unsigned Size>
struct FieldPolicy;

template <typename T, unsigned Dim>
class Position;

template <typename T, unsigned Dim>
class Velocity;

template <typename T, typename Func, unsigned Dim>
class ContinuousVelocity;


/** @brief Class for general physical fields.

    This class is used for representing scalar or vector fields, such as displacement and velocity fields.
    @tparam T Numeric data type of the field.
    @tparam Dim Dimension of the field (2 or 3).
    @tparam Size Length of the vector aassociated with each point in the field (e.g. 1 for scalar fields).
    */
template <typename T, unsigned Dim = 2, unsigned Size = 2>
class Field
{
    public:
        /** Constructor for initializating the Field.    
            @param nx The number of grid points in \f$x\f$-direction.
            @param ny The number of grid points in \f$y\f$-direction.
            */
        Field(unsigned nx, unsigned ny): nx_(nx), ny_(ny), data_(nx, ny), time_() {}

        /** Get all data values of the Field.
            @return A Tensor containing all data values of the Field.
            */
        inline auto& GetAll() const
        {
            return data_;
        }

        /** Get the time associated with the Field.
            @return Time of the Field.
            */
        inline auto GetTime() const
        {
            return time_;
        }

        /** Read data from a file.
            @param file_name Input data file name.
            */
        inline void ReadFromFile(const std::string& file_name)
        {
            FieldPolicy<T, Dim, Size>::ReadFromFile(*this, file_name);
        }

        /** Set all data vlues of the Field with an existing Tensor.
            @param data Input Tensor.
            */
        inline void SetAll(Tensor<Vector<T, Size>, Dim>& data)
        {
            data_ = data;
        }

        /** Update the time associated with the Field.
            @param time New time for the Field.
            */
        inline void UpdateTime(const T time)
        {
            time_ = time;
        }

        /** Write data to a file.
            @param file_name Output data file name.
            */
        inline void WriteToFile(const std::string& file_name) const
        {
            FieldPolicy<T, Dim, Size>::WriteToFile(*this, file_name);
        }

        friend struct FieldPolicy<T, Dim, Size>;

    protected:
        const unsigned nx_;/**<Number of grid points in \f$x\f$-direction.*/
        const unsigned ny_;/**<Number of grid points in \f$y\f$-direction.*/
        Tensor<Vector<T, Size>, Dim> data_;/**<Tensor containing all Field data values.*/
        T time_;/**<Time of the Field.*/
};

/** @brief Dimension-dependent implementations of Field class.
    */
template <typename T, unsigned Dim, unsigned Size>
struct FieldPolicy {};

/** @brief 2D implementations of Field class.
    
    This is a specialization of FieldPolicy template struct that for 2-dimension Field.
    */
template <typename T, unsigned Size>
struct FieldPolicy<T, 2, Size>
{
    /** Write Field data to a file.
        @param field Input field.
        @param file_name Output data file name.
        */
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
    
    /** Read data form a file to Field.
        @param field Field that needs to read the data.
        @param file_name Input data file name.
        */
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

/** @brief Field of particle positions.

    This class is used for representing positions (displacements) of particles. It is a subclass of Field.
    @tparam T Numeric data type of the field.
    @tparam Dim Dimension of the field (2 or 3).
    */
template <typename T, unsigned Dim = 2>
class Position : public Field<T, Dim, Dim>
{
    public:
        using vec = LCS::Vector<T, Dim>;

        /** Constructor for initializating the position field.    
            @param nx The number of points in \f$x\f$-direction.
            @param ny The number of points in \f$y\f$-direction.
            */
        Position(unsigned nx, unsigned ny): Field<T, Dim, Dim>(nx, ny) {}

        /** Set data values of the field using \f$x\f$- and \f$y\f$-coordinates. The grid is assumed to be a Cartesian grid when using this function. This is an overloaded function.
            @param xrange Vector that contains \f$x\f$-coordinates of the grid.
            @param yrange Vector that contains \f$y\f$-coordinates of the grid.
            */
        void SetAll(const std::vector<T>& xrange, const std::vector<T>& yrange)
        {
            // make sure sizes match
            if (xrange.size()!=this->nx_||yrange.size()!=this->ny_)
                throw std::domain_error("sizes do not match!");

            // fill in the values
            #pragma omp parallel for
            for (unsigned i = 0; i < this->nx_; ++i)
                for (unsigned j = 0; j < this->ny_; ++j)
                    this->data_(i,j) = vec(xrange[i], yrange[j]);

            pos_xrange_ = xrange;
            pos_yrange_ = yrange;
        }

        /** Set data values of the field using the coordinates of end points. The grid is assumed to be an uniform Cartesian grid when using this function. This is a overloaded function.
            @param xmin The minimum \f$x\f$-coordinate.
            @param xmax The maximum \f$x\f$-coordinate.
            @param ymin The minimum \f$y\f$-coordinate.
            @param ymax The maximum \f$y\f$-coordinate.
            */
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

        /** Set data values of the field. This function would class the function with the same name in the base class.
            */
        void SetAll(Tensor<vec, Dim>& data)
        {
            Field<T, Dim, Dim>::SetAll(data);
        }

        /** Get the positions of a grid point.
            @param i Index in \f$x\f$-coordinate.
            @param j Index in \f$y\f$-coordinate.
            @return Tuple contains \f$x\f$- and \f$y\f$-coordinates of a given grid point.
            */
        inline auto Get(const unsigned i, const unsigned j) const
        {
            return std::make_tuple(this->data_(i,j).x, this->data_(i,j).y);
        }

        /** Get the range of coordinates for a given axis.
            @param axis Integer that represents an axis (0 for \f$x\f$-axis, 1 for \f$y\f$-axis).
        */
        inline auto& GetRange(const unsigned axis)
        {
            // need to check the input and if ranges exist
            if (axis == 0)
                return pos_xrange_;
            else
                return pos_yrange_;
        }

        /** Update positions of each grid point using a velocity field. OpenMP is used here for accelerating the calculations.
            @param vel Velocity field contains the velocity information of each grid point.
            @param delta Time step for updating the positions.
            */
        void Update(Velocity<T, Dim>& vel, T delta)
        {
            #pragma omp parallel for
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

        /** Initialize the Tensor that contains out-of-bound information. All points are assumed to be inside the boundary.
            */
        inline void InitializeOutOfBoundTensor()
        {
            // default value is false
            out_of_bound_.reset(new Tensor<bool, Dim>(this->nx_, this->ny_));
        }

        /** Check if a point is out of the boundary.
            @param i Index in \f$x\f$-coordinate.
            @param j Index in \f$y\f$-coordinate.
            @return A bool that shows if the point is out of boundary.
            */
        inline bool IsOutOfBound(const unsigned i, const unsigned j)
        {
            if (out_of_bound_ != nullptr)
                return out_of_bound_->GetValue(i,j);
            else
                return false; // default value
        }

        /** Set the out boundary of the field.
            @param xmin The minimum \f$x\f$-coordinate of the out boundary.
            @param xmax The maximum \f$x\f$-coordinate of the out boundary.
            @param ymin The minimum \f$y\f$-coordinate of the out boundary.
            @param ymax The maximum \f$y\f$-coordinate of the out boundary.
            */
        inline void SetBound(const T& xmin, const T& xmax, const T& ymin, const T& ymax)
        {
            bound_xmin_ = xmin;
            bound_xmax_ = xmax;
            bound_ymin_ = ymin;
            bound_ymax_ = ymax;
        }

    private:
        std::unique_ptr<Tensor<bool, Dim>> out_of_bound_;/**<Pointer to a Tensor that shows if the points are out of boundary.*/
        std::vector<T> pos_xrange_;/**<Range of \f$x\$-coordinates.*/
        std::vector<T> pos_yrange_;/**<Range of \f$y\$-coordinates.*/

        T bound_xmin_;/**<The minimum \f$x\f$-coordinate of the out boundary.*/
        T bound_xmax_;/**<The maximum \f$x\f$-coordinate of the out boundary.*/
        T bound_ymin_;/**<The minimum \f$y\f$-coordinate of the out boundary.*/
        T bound_ymax_;/**<The maximum \f$y\f$-coordinate of the out boundary.*/
};

/** @brief Field of particle velocities.

    This class is used for representing velocities of particles. It is a subclass of Field.
    A Position field is associated with each Velocity field.
    @tparam T Numeric data type of the field.
    @tparam Dim Dimension of the field (2 or 3).
    */
template <typename T, unsigned Dim = 2>
class Velocity : public Field<T, Dim, Dim>
{
    public:
        using vec = LCS::Vector<T, Dim>;

        /** Constructor for initializating the velocity field.    
            @param nx The number of points in \f$x\f$-direction.
            @param ny The number of points in \f$y\f$-direction.
            @param pos The Position field that associated with the Velocity field.
            */
        Velocity(unsigned nx, unsigned ny, Position<T, Dim>& pos):
            Field<T, Dim, Dim>(nx, ny), pos_(pos) {}

        /** Get the velocities of a grid point.
            @param i Index in \f$x\f$-coordinate.
            @param j Index in \f$y\f$-coordinate.
            @return Tuple contains \f$x\f$- and \f$y\f$-coordinates of a given grid point.
            */
        inline auto Get(const unsigned i, const unsigned j) const
        {
            return std::make_tuple(this->data_(i,j).x, this->data_(i,j).y);
        }

        /** Get the Position field that associated with the Velocity field.
            @return Position field for this Velocity field.
            */
        inline auto& GetPosition()
        {
            return pos_;
        }

        /** Interploate another Velocity field to set the velocities of this field.
            @param ref_vel The reference Velocity field.
            */
        void InterpolateFrom(Velocity<T, Dim>& ref_vel)
        {
            // interpolation only works for orthogonal coordinates
            auto const ref_pos_x = ref_vel.GetPosition().GetRange(0);
            auto const ref_pos_y = ref_vel.GetPosition().GetRange(1);

            #pragma omp parallel for
            for (unsigned i = 0; i < this->nx_; ++i)
            {
                T pos_x, pos_y;
                Tensor<Vector<T, 2>, 2> ref_v(2,2);

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
        Position<T, Dim>& pos_; /**<Position field associated with this Velocity field*/
};

/** @brief Velocity field with continuous velocity function
    
    This class represents velocity fields that are defined by known continuous velocity function. It is a subclass of Velocity.
    @tparam T Numeric data type of the field.
    @tparam Func Continous velocity function for this field.
    @tparam Dim Dimension of the field (2 or 3).
    */
template <typename T, typename Func, unsigned Dim>
class ContinuousVelocity : public Velocity<T, Dim>
{
    public:
        using vec = LCS::Vector<T, Dim>;

        /** Constructor for initializating the continous velocity field.    
            @param nx The number of points in \f$x\f$-direction.
            @param ny The number of points in \f$y\f$-direction.
            @param pos The Position field that associated with the ContinuousVelocity field.
            @param time Initial time, default value is 0.
            */
        ContinuousVelocity(unsigned nx, unsigned ny, Position<T, Dim>& pos, T time=0):
            Velocity<T, Dim>(nx, ny, pos), f_()
        {
            this->time_ = time;
            SetAll();
        }

        /** Constructor for initializating the continous velocity field that the velocity function has parameters.
            @param nx The number of points in \f$x\f$-direction.
            @param ny The number of points in \f$y\f$-direction.
            @param pos The Position field that associated with the ContinuousVelocity field.
            @param parameters A vector that contains all the parameters of the velocity function.
            @param time Initial time, default value is 0.
            */
        ContinuousVelocity(unsigned nx, unsigned ny, Position<T, Dim>& pos,
                std::vector<T>& parameters, T time=0): Velocity<T, Dim>(nx, ny, pos), f_(parameters)
        {
            this->time_ = time;
            SetAll();
        }

        /** Use continuous function to set all velocity values.
            */
        void SetAll()
        {
            #pragma omp parallel for
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

        /** Get the velocity function of this field.
            @return Continuous velocity function of this field.
            */
        inline Func& Function()
        {
            return f_;
        }

    private:
        Func f_;/**<Continous velocity function of this field.*/
};

}
