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
            Set(xrange, yrange);
        }

        // getter
        auto Get(const unsigned i, const unsigned j) const
        {
            return std::make_tuple(data_[i][j].x, data_[i][j].y);
        }

    private:
        std::vector<std::vector<vec>> data_;
        unsigned nx_;
        unsigned ny_;
};
