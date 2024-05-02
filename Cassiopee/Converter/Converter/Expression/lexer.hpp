#ifndef _CONVERTER_EXPRESSION_LEXER_HPP_
#define _CONVERTER_EXPRESSION_LEXER_HPP_
#include <utility>
#include <unordered_map>
#include "Expression/ast.hpp"

namespace Expression
{
	class ast::lexer
	{
	public:
		enum kind
		{
			LEFT_BRACKET = 0,
			RIGHT_BRACKET= 1,
			FUNCTION     = 2,
			OPERATOR     = 3,
			NUMBER       = 4,
			VARIABLE     = 5,
			UNARY_OPERATOR=6,
			PI_VALUE      =7,
			SEPARATOR     =8, // la virgule pour min et max
			MINMAXOP      =9,

			UNKNOWN      = -1			
		};
		using token=std::pair<std::string,kind>;
		using container=std::vector<token>;
		using iterator =container::iterator;
		using const_iterator = container::const_iterator;

		lexer() : m_tokens{} {}

		void analyze(const std::string& expr );

		iterator begin() { return m_tokens.begin(); }
		const_iterator cbegin() const { return m_tokens.cbegin(); }
		const_iterator begin() const { return cbegin(); }

		iterator end() { return m_tokens.end(); }
		const_iterator cend() const { return m_tokens.cend(); }
		const_iterator end() const { return cend(); }

		std::size_t size() const { return m_tokens.size(); }
		const token& operator [] ( std::size_t ind ) const { return m_tokens[ind]; }
	private:
		container m_tokens;
	};
}

#endif