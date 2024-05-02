#ifndef _EXPRESSIONTREE_NODE_HPP_
#  define _EXPRESSIONTREE_NODE_HPP_
#  include <string>
#  include <memory> 
#  include "Expression/ast.hpp"

namespace Expression
{
	class simd_vector_wrapper;
	class ast::node 
	{
	public:
		/**
		 * @brief Evaluate the node at index i
		 * @details evaluate the node at index i. If this node is not a terminal node
		 * of the abstract syntax tree, evaluate left ( if any ) and right child before
		 * applying the operation on the node ( function or operator ... )
		 * 
		 * @param i The index for vector variables
		 * @return The value of the expression
		 */
		virtual double operator() ( std::size_t i ) const = 0;

		/**
		 * @brief Evaluate the node for a pack of values ( a vector ) for vectorization
		 * @details Evaluate the node for a tiny vector issued from variables or from computation
		 *          coming from the children.
		 * 
		 * @param i The beginning index for the tiny vectors
		 * @return A tiny vector with temporary values ( or pointing on variable )
		 */
		virtual simd_vector_wrapper eval_simd( std::size_t i ) const = 0;
		/**
		 * @brief The representation in string of the node
		 * @details Return the representation in string of the node for debuggin or to display the abstract
		 * syntax tree in a human form. If the node is not a terminal node, call before the conversion of left
		 * child ( if any ) which is concataned with current node representation and call the conversion of
		 * the right child which is appended at the string.
		 * @return Return the representation of the node and his children...
		 */
		explicit virtual operator std::string() const = 0;

		/**
		 * @brief      Simplification of a node with these children
		 * 
		 */
		virtual std::shared_ptr<node> optimize() = 0;

		virtual std::shared_ptr<node> derivate( ) = 0;
	};
}

#endif

