#ifndef _CONVERTER_EXPRESSION_CONSTANTE_HPP_
#define _CONVERTER_EXPRESSION_CONSTANTE_HPP_
# include <string>
# include "Expression/node.hpp"
# include "Expression/simd_vector_wrapper.hpp"

namespace Expression
{
	class constante : public ast::node
	{
	public:
		constante(double val) : m_value(val)
		{}

		double operator() ( std::size_t i ) const final
		{
			return m_value;
		}

		simd_vector_wrapper eval_simd( std::size_t i ) const final
		{
			return {m_value};
		}

		explicit operator std::string() const final
		{
		  return std::to_string(m_value);
		  //char buffer[256];
		  //sprintf(buffer, "%g", m_value);
		  //return std::string(buffer);
		}

		virtual std::shared_ptr<node> derivate( ) final
		{
			return std::make_shared<constante>(0.);
		}

		std::shared_ptr<node> optimize() final
		{
			return std::make_shared<constante>(m_value);
		}
	private:
		double m_value;
	};
}

#endif
