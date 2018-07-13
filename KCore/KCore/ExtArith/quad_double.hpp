#ifndef _MATH_QUADDOUBLE_HPP_
#define _MATH_QUADDOUBLE_HPP_
# include <cmath>
namespace ExtendedArithmetics
{
    /**
     * @brief     Extended floating point type : nearly double precision floating numbers
     * 
     * Reference : A Floating-Point technique for Extending the Available Precision, Author T.J. Dekker, Numerical Mathematics, March 18th 1971
     *             Springer-Verlag
     */
    class quad_double {
    public:
        constexpr quad_double() : m_val0(0.), m_val1(0.) {  }
        constexpr quad_double(double x) : m_val0(x), m_val1(0.) { }
        constexpr quad_double(double x, double xs ): m_val0(x), m_val1(xs) {}
        constexpr quad_double( const quad_double& x ) : m_val0(x.m_val0), m_val1(x.m_val1) { }
        ~quad_double() = default;

#ifdef __INTEL_COMPILER        
#       pragma optimize("", off)
#endif
        quad_double operator + ( const quad_double& x ) const {
            double r = m_val0 + x.m_val0;
            double s = ( std::abs(m_val0) > std::abs(x.m_val0) ? m_val0 - r + x.m_val0 + x.m_val1 + m_val1 :
                                                                 x.m_val0 - r + m_val0 + m_val1 + x.m_val1 );
            double z = r + s;                                           
            return quad_double(z, r - z + s );
        }

        quad_double operator + ( double x ) const {
            return (*this) + quad_double(x);
        }

#ifdef __INTEL_COMPILER    
#       pragma optimize("", off)
#endif
        quad_double operator - ( const quad_double& x ) const {
            double r = m_val0 - x.m_val0;
            double s = ( std::abs(m_val0) > std::abs(x.m_val0) ? m_val0 -r - x.m_val0 - x.m_val1 + m_val1 :
                                                                -x.m_val0 -r + m_val0 + m_val1 - x.m_val1 );  
            double z = r + s;
            return quad_double(z, r-z+s);
        }
        quad_double operator - ( double x ) const {
            return (*this) - quad_double(x);
        }
        quad_double operator - () const {
            return quad_double(-m_val0,-m_val1);
        }
#ifdef __INTEL_COMPILER        
#       pragma optimize("", off)
#endif 
        quad_double operator * ( const quad_double& x ) const {
            quad_double c = mul_exact_(m_val0, x.m_val0);
            c.m_val1 = m_val0 * x.m_val1 + m_val1 * x.m_val0 + c.m_val1;
            double z = c.m_val0 + c.m_val1;
            return quad_double(z, c.m_val0 - z + c.m_val1 );        
        }
        quad_double operator * ( double x ) const {
            return (*this) * quad_double(x);
        }

#ifdef __INTEL_COMPILER
#       pragma optimize("", off)
#endif 
        quad_double operator / ( const quad_double& x ) const {
             double c = m_val0/x.m_val0;
             quad_double u = mul_exact_( c, x.m_val0);
             double cc = ( m_val0 - u.m_val0 - u.m_val1 + m_val1 - c * x.m_val1 )/x.m_val0;
             double z = c + cc;
             return quad_double( z , c - z + cc );
        }
        quad_double operator / ( double x ) const {
            return (*this) / quad_double(x);
        }

        quad_double sqrt() const {
            double c = std::sqrt(m_val0);
            if (m_val0 > 0)
            {
                 quad_double u = mul_exact_(c,c);
                 double cc = (m_val0-u.m_val0-u.m_val1+m_val1)*0.5/c;
                 double y = c + cc;
                 return quad_double( y, c - y + cc );
            } else return quad_double(0.,0.);
        }
                        
        double get_double_precision() const { return m_val0; }
        double& get_double_precision() { return m_val0; }
        double get_extended_precision() const { return m_val1; }
        double& get_extended_precision() { return m_val1; }
        
        operator double() const { return m_val0; }

    private:
        quad_double mul_exact_( double x, double y ) const
        {
            // 26 : nb bytes mantisse - (nb bytes mantisse / 2 )        
            const double c = (1L<<26)+1;
            double p  =  x * c;
            double hx = x - p + p, tx = x - hx;
            p = y * c;
            double hy = y - p + p, ty = y - hy;
            p = hx * hy;
            double q = hx * ty + tx * hy;
            double z = p + q;
            return quad_double( z, p-z+q+tx*ty );
        }

        double m_val0, m_val1;
    };
    inline quad_double sqrt( const quad_double& x ) {
        return x.sqrt();
    }
    inline quad_double operator + ( double x, const quad_double& qy ) {
        return quad_double(x) + qy;
    }
    inline quad_double operator - ( double x, const quad_double& qy ) {
        return quad_double(x) - qy;
    }
    inline quad_double operator * ( double x, const quad_double& qy ) {
        return quad_double(x) * qy;
    }
    inline quad_double operator / ( double x, const quad_double& qy ) {
        return quad_double(x) / qy;
    }
}

#endif
