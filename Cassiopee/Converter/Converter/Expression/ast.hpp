#ifndef _CONVERTER_EXPRESSION_AST_HPP_
#define _CONVERTER_EXPRESSION_AST_HPP_
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include "Memory/vector_view.hpp"
#include "Expression/symbol_table.hpp"
using K_MEMORY::vector_view;
namespace Expression
{
/**
 * @brief Abstract Syntax tree
 * @details Abstract syntax tree to evaluate litteral expression and computing derivate of these expressions
 * 
 * Usage:
 * ------
 *     std::vector<double> arr1={1.,2.,3.};
 *     std::vector<double> arr2={-2.,3.,5.};
 *     std::vector<double> zeros={0.,0.,0.};
 *     std::vector<double> ones ={1.,1.,1.};
 *     double scal1 = 1.5, scal2 = 3.2, zero = 0.;
 *     Expression::ast expr1("{a}*{x}+{b}*{y}");// The expression
 *     auto res = expr1({ {"x", arr1}, {"y", arr2}, {"scal1", {1,&scal1}}, {"scal2", {1,&scal2}}});// Evaluating
 *     for ( auto v : res ) std::cout << v << " ";
 *     std::cout << std::endl;
 *     auto d_expr1 = expr1.derivate();// Compute the derivate
 *     auto d_res = d_expr1({ {"x", arr1}, {"y", arr2}, {"scal1", {1,&scal1}}, {"scal2", {1,&scal2}},
 *                            {"d_x", ones}, {"d_y", zeros}, {"d_a", {1,&zero}}, {"d_b", {1,&zero}}});     
 *     for ( auto v : d_res ) std::cout << v << " ";
 *     std::cout << std::endl;
 *
 */
class ast
{
public:

    using symbol_table_t = Expression::symbol_table;

    class lexer;   // Le lexer pour decouper notre chaine de caractere en token
    class parser; // Le parser pour transformer la liste de token en arbre de token

    class node; // Noeud d'un arbre AST
    class tree; // Arbre de syntaxe abstraite utilisant le parser pour cela et ayant des node pour noeuds de l'arbre.

    /**
     * @brief Build the ast from the expression given in the string expr
     * @details Build the abstract syntax tree for the expression given in the
     *          string expr. If the expression contains unknown symbols or expression,
     *          return a logic_error exception.
     * 
     * @param expr The string containing the expression
     */
    ast(const std::string &expr);

    /**
     * @brief Build an ast object from a new abstract syntax tree
     * @details Build an ast object from an internal representation of the
     *          abstract syntax tree. The intent of this method is to be
     *          used only internal in this function.
     * 
     * @param tr The internal abstract syntax tree.
     */
    ast( std::shared_ptr<tree> tr );

    /**
     * @brief Evaluate expression for some values
     * @details Evaluate expression with some values defined in a dictionnary.
     * 
     * @param params The dictionnary of the values
     * @return The value which is the result of the expression evaluation. This value is always a vector,
     *         with dimension one if the result is a scalar.
     */
    std::vector<double> operator( )(const std::unordered_map<std::string, vector_view<double>> &params);
    /**
     * @brief Evaluate expression for some values
     * @details Evaluate expression with some values defined in a dictionnary and store the result in the
     * output parameter.
     * 
     * @param params The dictionnary of the values
     */

    void eval(const std::unordered_map<std::string, vector_view<double>> &params, vector_view<double>& output);
    /**
     * @brief Evaluate the expression
     * @details Evaluate the expression. The result must be written in the expression with the assign operator ( symbol = )
     *          None result is returned.
     * 
     * @param params The dictionnary of the values and variables used.
     */
    void eval(const std::unordered_map<std::string, vector_view<double>> &params);

    /**
     * @brief Return the derivate of the expression
     * @details Return the derivate of the expression, addind some new variables in the table of symbol
     *          for derivates of the symbols. The name of these new variables is the name of the variable
     *          prefixed with a "d_" characters. By example, from variable "x", this method create the
     *          new variable "d_x" if this variable doesn"t exist in the table of the symbol.
     * 
     * @return The abstract syntax tree of the derivate.
     */
    ast derivate( );

    /**
     * @brief Return a string with human readable expression
     * @details Return a string with the representation of the expression for
     *          human. 
     * @return The string which contains the representation of the ast object
     */
    std::string to_string() const;


    /**
     * @brief Convert an ast object in string which purpose is to be readable by human.
     * @details Convert an ast object in string. This string is a representation of the
     *          abstract syntax tree which can be read easily by a human.
     * @return The string representation of the ast object.
     */
    explicit operator std::string() const { return to_string(); }

private:
    std::shared_ptr<tree> m_pt_tree;
};

}

#endif





