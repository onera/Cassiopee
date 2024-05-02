#ifndef _CONVERTER_EXPRESSION_PARSER_HPP_
# define _CONVERTER_EXPRESSION_PARSER_HPP_
# include <vector>
# include <unordered_map>
# include <memory>

#include "Expression/lexer.hpp"

namespace Expression
{
	class ast::parser 
	{
	public:
		using token = lexer::token;
		using symbol_table_t = std::vector<std::string>;
		/**
		 * @brief      Arbre des tokens permettant de creer ensuite l'ast
		 */
		class tree {
		public:
			enum child_type
			{
				left_child = 0, right_child=1
			};

			tree( const token& tok ) : m_token(tok)
			{}

			std::shared_ptr<tree>& operator [] ( child_type ch )
			{
				return m_pt_children[ch];
			}

			const std::shared_ptr<tree>& operator [] ( child_type ch ) const
			{
				return m_pt_children[ch];
			}

			const token& operator *() const { return m_token; }

			std::string display_tree() const;

		private:
			std::shared_ptr<tree> m_pt_children[2];
			token m_token;
		};

		using tree_pointer = std::shared_ptr<tree>;
		using container=tree_pointer;

		parser( const lexer& lex );
		parser( const lexer::const_iterator& beg, const lexer::const_iterator& end );
		parser( const lexer::const_iterator& beg, std::size_t lgth = 0 );

		std::size_t size_of_symbol_table() const { return m_symbol_table.size(); }
		const symbol_table_t& get_symbol_table() { return m_symbol_table; }

		const tree_pointer& getTree() const { return m_pt_tree; }

		explicit operator std::string() const {
			std::string s("\nSymbol table: ");
			for ( auto& item : m_symbol_table ) s += item + " ";
			s += "\n";
			return m_pt_tree->display_tree() + s; 
		}

	private:
		void parse( const lexer::const_iterator& beg, const lexer::const_iterator& end, tree_pointer& t, bool is_neg = false );
		tree_pointer m_pt_tree;
		symbol_table_t m_symbol_table;
	};
}

# endif