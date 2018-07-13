#ifndef _CONVERTER_EXPRESSION_UNARY_OPERATOR_HPP_
# define _CONVERTER_EXPRESSION_UNARY_OPERATOR_HPP_
# include <memory>
# include <cmath>
# include <typeinfo>
# include "Expression/node.hpp"
# include "Expression/constante.hpp"

namespace Expression
{
	/**
	 * @brief Arithmetic binary operators
	 * @details Arithmetics operators having left and right statements.
	 * 
	 */
	class unary_operator : public ast::node
	{
	public:
		enum OP {
			ADD = 0, SUB
		};
		unary_operator( OP op, std::shared_ptr<ast::node> r_stmt ) :
			m_right_statement(r_stmt), m_op(op)
		{}

		double operator() ( std::size_t i ) const final
		{
			switch(m_op)
			{
			case ADD:
				return (*m_right_statement)(i);
				break;
			case SUB:
				return -(*m_right_statement)(i);
				break;
			}
			return 0.;
		}

		simd_vector_wrapper eval_simd ( std::size_t i ) const final
		{
			switch(m_op)
			{
			case ADD:
				return m_right_statement->eval_simd(i);
				break;
			case SUB:
				return -m_right_statement->eval_simd(i);
				break;
			}
			return simd_vector_wrapper{};
		}

		explicit operator std::string() const final
		{
			switch(m_op)
			{
			case ADD:
				return "+"+std::string(*m_right_statement);
				break;
			case SUB:
				return "-"+std::string(*m_right_statement);
				break;
			}
			return std::string("");
		}

		virtual std::shared_ptr<node> derivate( ) final
		{
			switch(m_op)
			{
			case ADD:
				return m_right_statement->derivate();
				break;
			case SUB:
				return std::make_shared<unary_operator>(SUB,m_right_statement->derivate());
				break;
			}
			return nullptr;
		}

		std::shared_ptr<node> optimize() final
		{
			m_right_statement = m_right_statement->optimize();
			// SI + unaire, on ne fait rien avec...
			if (m_op == ADD) return m_right_statement;
			auto& r = *m_right_statement;
			if (typeid(r).hash_code() == typeid(constante).hash_code() ) {
				return std::make_shared<constante>(-(*m_right_statement)(0));
			}
			else
				return std::make_shared<unary_operator>(m_op, m_right_statement);
		}
	private:
	    std::shared_ptr<ast::node> m_right_statement;		
	    OP m_op;
	};

}
#endif

