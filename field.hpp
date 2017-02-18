#include <vector>
#include <algorithm>
#include <tuple>
#include <cassert>
#include <stdexcept>

// vector
template <typename T, unsigned Dim = 2>
struct Vector
{
    Vector(): x(), y() {}
    Vector(const T& x, const T& y): x(x), y(y) {}

    T x, y;
};

// poisiton field
template <typename T, unsigned Dim = 2>
class Position
{
    public:
        using vec = Vector<T, Dim>;

        // constructor
        Position(unsigned nx, unsigned ny): nx_(nx), ny_(ny),
            data_(nx, std::vector<vec>(ny, vec())) {}

        // setter for all values in the field
        void SetAll(const std::vector<T>& xrange, const std::vector<T>& yrange)
        {
            // make sure dimensions match
            if (xrange.size()!=nx_||yrange.size()!=ny_)
                throw std::domain_error("dimensions do not match!");

            // fill in the values
            auto itx = xrange.begin();
            std::generate(data_.begin(), data_.end(), [&]
                {
                    auto ity = yrange.begin();
                    std::vector<vec> v(ny_, vec());
                    std::generate(v.begin(), v.end(),
                        [&]{ return vec(*itx, *(ity++)); });
                    ++itx;
                    return v;
                });
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

        // getter
        auto Get(const unsigned i, const unsigned j) const
        {
            return std::make_tuple(data_[i][j].x, data_[i][j].y);
        }

    private:
        std::vector<std::vector<vec>> data_;
        const unsigned nx_;
        const unsigned ny_;
};

// velocity field
template <typename T, unsigned Dim>
class AnalyticVelocity;

template <typename T, unsigned Dim = 2>
class Velocity
{
    public:
        using vec = Vector<T, Dim>;

        // constructor
        Velocity(unsigned nx, unsigned ny): nx_(nx), ny_(ny),
            data_(nx, std::vector<vec>(ny, vec())) {}

        // setter
        void SetAll(std::vector<std::vector<vec>> data)
        {
            data_ = data;
        }

        // getter
        auto Get(const unsigned i, const unsigned j) const
        {
            return std::make_tuple(data_[i][j].x, data_[i][j].y);
        }

        friend class AnalyticVelocity<T, Dim>;

    private:
        std::vector<std::vector<vec>> data_;
        const unsigned nx_;
        const unsigned ny_;
};

// velocity field with analytic velocity function
template <typename T, unsigned Dim>
class AnalyticVelocity : public Velocity<T, Dim>
{
    public:
        using vec = Vector<T, Dim>;
        using func = std::tuple<T, T>(T, T);

        AnalyticVelocity(unsigned nx, unsigned ny, func& f):
            Velocity<T, Dim>(nx, ny), f_(f) {}

        // use analytic function to set all velocity values
        void SetAll(const Position<T, Dim>& pos)
        {
            for (unsigned i = 0; i < this->nx_; ++i)
            {
                for (unsigned j = 0; j < this->ny_; ++j)
                {
                    T x, y, vx, vy;
                    std::tie(x, y) = pos.Get(i, j);
                    std::tie(vx, vy) = f_(x, y);
                    this->data_[i][j] = vec(vx, vy);
                    auto a = this->Get(0, 0);
                }
            }
        }

    private:
        func& f_;
};
