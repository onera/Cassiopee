#ifndef _CONVERTER_EXPRESSION_MATH_FUNCTIONS_HPP_
#  define _CONVERTER_EXPRESSION_MATH_FUNCTIONS_HPP_
#  include <cmath>
#  include "Expression/function.hpp"

namespace Expression
{
	class sin;

	class cos : public function
	{
	public:
		cos( std::shared_ptr<node> r_stmt ) : function(r_stmt) {}
	private:
		virtual std::string name() const final { return "cos"; }
		virtual double apply_to(double x) const final { return std::cos(x); }
		virtual simd_vector_wrapper apply_to( const simd_vector_wrapper& u ) const  final
		{
			simd_vector_wrapper w;
			double* __restrict__  dw = w.data();
			const double* __restrict__ du = u.data();
#           pragma omp simd			
			for ( std::size_t i = 0; i < simd_vector_wrapper::max_size; ++i )
				dw[i] = std::cos(du[i]);
			return w;
		}

		virtual std::shared_ptr<ast::node> df() final;
    	std::shared_ptr<function> clone( std::shared_ptr<node>& right_stmt ) const
    	{
    		return std::make_shared<cos>( right_stmt );
    	}
        std::shared_ptr<node> simplify() final;

	};

	class sin : public function
	{
	public:
		sin( std::shared_ptr<node> r_stmt ) : function(r_stmt) {}
	private:
		virtual std::string name() const final { return "sin"; }
		virtual double apply_to(double x) const final { return std::sin(x); }
		virtual simd_vector_wrapper apply_to( const simd_vector_wrapper& u ) const   final
		{
			simd_vector_wrapper w;
			double* __restrict__ dw = w.data();
			const double* __restrict__ du = u.data();
# pragma omp simd
			for ( std::size_t i = 0; i < simd_vector_wrapper::max_size; ++i )
				dw[i] = std::sin(du[i]);
			return w;
		}
		virtual std::shared_ptr<ast::node> df()  final;
    	std::shared_ptr<function> clone( std::shared_ptr<node>& right_stmt ) const
    	{
    		return std::make_shared<sin>( right_stmt );
    	}
        std::shared_ptr<node> simplify() final;
	};

	class tan : public function
	{
	public:
		tan( std::shared_ptr<node> r_stmt ) : function(r_stmt) {}
	private:
		virtual std::string name() const final { return "tan"; }
		virtual double apply_to(double x) const final { return std::tan(x); }
		virtual simd_vector_wrapper apply_to( const simd_vector_wrapper& u ) const   final
		{
			simd_vector_wrapper w;
			double* __restrict__ dw = w.data();
			const double* __restrict__ du = u.data();
# pragma omp simd
			for ( std::size_t i = 0; i < simd_vector_wrapper::max_size; ++i )
				dw[i] = std::tan(du[i]);
			return w;
		}
		virtual std::shared_ptr<ast::node> df()  final;
    	std::shared_ptr<function> clone( std::shared_ptr<node>& right_stmt ) const
    	{
    		return std::make_shared<tan>( right_stmt );
    	}
        std::shared_ptr<node> simplify() final;
	};

	class ln : public function
	{
	public:
		ln( std::shared_ptr<node> r_stmt ) : function(r_stmt) {}
	private:
		virtual std::string name() const final { return "ln"; }
		virtual double apply_to(double x) const final { return std::log(x); }
		virtual simd_vector_wrapper apply_to( const simd_vector_wrapper& u ) const   final
		{
			simd_vector_wrapper w;
			double* __restrict__ dw = w.data();
			const double* __restrict__ du = u.data();
# pragma omp simd
			for ( std::size_t i = 0; i < simd_vector_wrapper::max_size; ++i )
				dw[i] = std::log(du[i]);
			return w;
		}
		virtual std::shared_ptr<ast::node> df()  final;
    	std::shared_ptr<function> clone( std::shared_ptr<node>& right_stmt ) const
    	{
    		return std::make_shared<ln>( right_stmt );
    	}
        std::shared_ptr<node> simplify() final;
	};

	class exp : public function
	{
	public:
		exp( std::shared_ptr<node> r_stmt ) : function(r_stmt) {}
	private:
		virtual std::string name() const final { return "exp"; }
		virtual double apply_to(double x) const final { return std::exp(x); }
		virtual simd_vector_wrapper apply_to( const simd_vector_wrapper& u ) const   final
		{
			simd_vector_wrapper w;
			double* __restrict__ dw = w.data();
			const double* __restrict__ du = u.data();
# pragma omp simd
			for ( std::size_t i = 0; i < simd_vector_wrapper::max_size; ++i )
				dw[i] = std::exp(du[i]);
			return w;
		}
		virtual std::shared_ptr<ast::node> df()  final
		{
			return std::make_shared<exp>(m_right_statement);
		}
    	std::shared_ptr<function> clone( std::shared_ptr<node>& right_stmt ) const
    	{
    		return std::make_shared<exp>( right_stmt );
    	}
        std::shared_ptr<node> simplify() final;
	};

	class sqrt : public function
	{
	public:
		sqrt( std::shared_ptr<node> r_stmt ) : function(r_stmt) {}
	private:
		virtual std::string name() const final { return "sqrt"; }
		virtual double apply_to(double x) const final { return std::sqrt(x); }
		virtual simd_vector_wrapper apply_to( const simd_vector_wrapper& u ) const   final
		{
			simd_vector_wrapper w;
			double* __restrict__ dw = w.data();
			const double* __restrict__ du = u.data();
# pragma omp simd
			for ( std::size_t i = 0; i < simd_vector_wrapper::max_size; ++i )
				dw[i] = std::sqrt(du[i]);
			return w;
		}
		virtual std::shared_ptr<ast::node> df()  final;

    	std::shared_ptr<function> clone( std::shared_ptr<node>& right_stmt ) const
    	{
    		return std::make_shared<sqrt>( right_stmt );
    	}
        std::shared_ptr<node> simplify() final;
	};

	class cosh : public function
	{
	public:
		cosh( std::shared_ptr<node> r_stmt ) : function(r_stmt) {}
	private:
		virtual std::string name() const final { return "cosh"; }
		virtual double apply_to(double x) const final { return std::cosh(x); }
		virtual simd_vector_wrapper apply_to( const simd_vector_wrapper& u ) const   final
		{
			simd_vector_wrapper w;
			double* __restrict__ dw = w.data();
			const double* __restrict__ du = u.data();
# pragma omp simd
			for ( std::size_t i = 0; i < simd_vector_wrapper::max_size; ++i )
				dw[i] = std::cosh(du[i]);
			return w;
		}
		virtual std::shared_ptr<ast::node> df()  final;

    	std::shared_ptr<function> clone( std::shared_ptr<node>& right_stmt ) const
    	{
    		return std::make_shared<cosh>( right_stmt );
    	}
        std::shared_ptr<node> simplify() final;
	};

	class sinh : public function
	{
	public:
		sinh( std::shared_ptr<node> r_stmt ) : function(r_stmt) {}
	private:
		virtual std::string name() const final { return "sinh"; }
		virtual double apply_to(double x) const final { return std::sinh(x); }
		virtual simd_vector_wrapper apply_to( const simd_vector_wrapper& u ) const   final
		{
			simd_vector_wrapper w;
			double* __restrict__ dw = w.data();
			const double* __restrict__ du = u.data();
# pragma omp simd
			for ( std::size_t i = 0; i < simd_vector_wrapper::max_size; ++i )
				dw[i] = std::sinh(du[i]);
			return w;
		}
		virtual std::shared_ptr<ast::node> df()  final;

    	std::shared_ptr<function> clone( std::shared_ptr<node>& right_stmt ) const
    	{
    		return std::make_shared<sinh>( right_stmt );
    	}
        std::shared_ptr<node> simplify() final;
	};

	class abs : public function
	{
	public:
		abs( std::shared_ptr<node> r_stmt ) : function(r_stmt) {}
	private:
		virtual std::string name() const final { return "abs"; }
		virtual double apply_to(double x) const final { return std::abs(x); }
		virtual simd_vector_wrapper apply_to( const simd_vector_wrapper& u ) const   final
		{
			simd_vector_wrapper w;
			double* __restrict__ dw = w.data();
			const double* __restrict__ du = u.data();
# pragma omp simd
			for ( std::size_t i = 0; i < simd_vector_wrapper::max_size; ++i )
				dw[i] = std::abs(du[i]);
			return w;
		}
		virtual std::shared_ptr<ast::node> df()  final;

    	std::shared_ptr<function> clone( std::shared_ptr<node>& right_stmt ) const
    	{
    		return std::make_shared<abs>( right_stmt );
    	}
        std::shared_ptr<node> simplify() final;
	};

	class sign : public function
	{
	public:
		sign( std::shared_ptr<node> r_stmt ) : function(r_stmt) {}
	private:
		virtual std::string name() const final { return "sign"; }
		virtual double apply_to(double x) const final { return (0. < x) - (x < 0.); }
		virtual simd_vector_wrapper apply_to( const simd_vector_wrapper& u ) const   
		{
			simd_vector_wrapper w;
			double* __restrict__ dw = w.data();
			const double* __restrict__ du = u.data();
# pragma omp simd
			for ( std::size_t i = 0; i < simd_vector_wrapper::max_size; ++i )
				dw[i] = (0. < du[i]) - (du[i] < 0.);
			return w;
		}
		virtual std::shared_ptr<ast::node> df()  ;

    	std::shared_ptr<function> clone( std::shared_ptr<node>& right_stmt ) const
    	{
    		return std::make_shared<sign>( right_stmt );
    	}
        std::shared_ptr<node> simplify() ;
	};

	void init_math_functions();
}

#endif



