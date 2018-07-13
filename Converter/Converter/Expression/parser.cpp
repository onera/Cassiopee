#include <algorithm>
#include <cassert>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include "Expression/parser.hpp"

namespace
{
Expression::ast::lexer::const_iterator find_closing_bracket(const Expression::ast::lexer::const_iterator &it,
                                                            const Expression::ast::lexer::const_iterator &end)
{
    using Expression::ast;
    int  bracket_level = 1;
    auto itC           = it;
    do {
        itC++;
        if (itC->second == ast::lexer::LEFT_BRACKET) bracket_level++;
        if (itC->second == ast::lexer::RIGHT_BRACKET) bracket_level--;
    } while ((bracket_level != 0) && (itC != end));
    return itC;
}

/**
 * Find a token in same level as beg and end iterator ( levels are defined by brackets ) */
Expression::ast::lexer::const_iterator find_statement(const Expression::ast::lexer::const_iterator &beg,
                                                      const Expression::ast::lexer::const_iterator &end,
                                                      const Expression::ast::lexer::token &         tok)
{
    using namespace Expression;
    ast::lexer::const_iterator itCur = beg;
    for (; itCur < end; ++itCur) {
        if (itCur->second == ast::lexer::LEFT_BRACKET) itCur = find_closing_bracket(itCur, end);
        if ((itCur->second == tok.second) && (itCur->first == tok.first)) break;
    }
    return itCur;
}

/**
 * Find a token in same level as beg and end iterator ( levels are defined by brackets^) */
Expression::ast::lexer::const_iterator find_next_function(const Expression::ast::lexer::const_iterator &beg,
                                                          const Expression::ast::lexer::const_iterator &end)
{
    using namespace Expression;
    ast::lexer::const_iterator itCur = beg;
    for (; itCur < end; ++itCur) {
        if (itCur->second == ast::lexer::LEFT_BRACKET) itCur = find_closing_bracket(itCur, end);
        if (itCur->second == ast::lexer::FUNCTION) break;
    }
    return itCur;
}
} // end Namespace

namespace Expression
{
std::string ast::parser::tree::display_tree( ) const
{
    std::string d = "[ ";
    if ((*this)[left_child]) d += (*this)[left_child]->display_tree( );
    d += m_token.first;
    if ((*this)[right_child]) d += (*this)[right_child]->display_tree( );
    d += "]";
    return d;
}

ast::parser::parser(const lexer &lex)
{
    // On va chercher les opérateurs/fonctions/etc. par ordre de précédence inverse :
    ast::lexer::const_iterator begLex = lex.begin( );
    ast::lexer::const_iterator endLex = lex.end( );
    parse(begLex, endLex, m_pt_tree);
    std::sort(m_symbol_table.begin(),m_symbol_table.end());
    auto last = std::unique(m_symbol_table.begin(), m_symbol_table.end());
    m_symbol_table.erase(last,m_symbol_table.end());
}

void ast::parser::parse(const lexer::const_iterator &beg, const lexer::const_iterator &end, tree_pointer &tr, bool is_neg)
{
	assert(beg != end);
    // Première chose à faire : voir si l'expression qu'on examine est entièrement encadrée par des parenthèses ou non :
    // Dans ce cas, on supprime cette paire de parenthèse pour analyser l'expression sous-jacente :
    if (beg->second == lexer::LEFT_BRACKET) {
        auto itRBracket = find_closing_bracket(beg, end);
        if (itRBracket == std::prev(end)) {
            parse(std::next(beg), std::prev(end), tr, is_neg);
            return;
        }
    }
    // Viens la virgule qui sépare deux statements gauches et droite ( pour le min et le max ! )
    auto it_separator = find_statement(beg, end, {",", lexer::SEPARATOR});
    if (it_separator != end) // Oui, il y a bien un separateur :
    {
        // Le pere est toujours l'operateur min ou max...
        if ((*(*tr)).second != lexer::MINMAXOP)
            throw std::runtime_error("Separator , is used only for min, max and logical functions nowaday !");
        parse(beg, it_separator, (*tr)[tree::left_child], is_neg);
        parse(std::next(it_separator), end, (*tr)[tree::right_child], is_neg);
        return;
    }
    // On cherche ensuite si l'expression possede un operateur d'assignement ( plus petite précédence ) qui n'est pas
    // inclu dans des parenthèses
    auto it_assignment = find_statement(beg, end, {"=", lexer::OPERATOR});
    if (it_assignment != end) // Oui, il y a bien un opérateur d'assignement :
    {
        // On vérifie qu'à gauche de =, c'est bien une variable prête à recevoir le résultat de l'expression à droite de
        // l'opérateur =
        if (std::prev(it_assignment)->second != lexer::VARIABLE)
            throw std::runtime_error("Variable waited at the left of the assignment operator !");
        tr = std::make_shared<tree>(*it_assignment); // Père
        // Parcourt fils gauche - fils droit
        (*tr)[tree::left_child] = std::make_shared<tree>(*std::prev(it_assignment));
        parse(std::next(it_assignment), end, (*tr)[tree::right_child]);
        return;
    }
    // Puis si elle possède ( en dehors des parenthèses ) un symbole == ou != ( test égalité )
    for (std::string s : {"==", "!=", ">", ">=", "<", "<=", "+", "-", "*", "/", "**"}) {
        auto it_sparse = find_statement(beg, end, {s, lexer::OPERATOR});
        if (it_sparse != end) { // Il y a bien un opérateur binaire de ce type
            {
                // Si opérateur + ou -, regarder si il n'est pas unaire :
                // Il est unaire si : en début d'expression, après une parenthèse ouvrante, après un autre opérateur
                if ((s == "+") || (s == "-")) {
                    // On vérifie si il est unaire et dans ce cas ignorer cet opérateur
                    if (it_sparse == beg) continue;                    
                    auto it_prev = std::prev(it_sparse);
                    if ((it_prev->second == lexer::LEFT_BRACKET) || (it_prev->second == lexer::OPERATOR)) continue;
                    // L'opérateur est binaire. Pour le -, il va falloir travailler les - successifs :

                }
                //tr = std::make_shared<tree>(*it_sparse); // Père
                if ( ( s == "-" ) and is_neg )
                    tr = std::make_shared<tree>(token{"+",lexer::OPERATOR}); // Père
                else
                    tr = std::make_shared<tree>(*it_sparse); // Père
                parse(beg, it_sparse, (*tr)[tree::left_child], is_neg);
                if ( s == "-" ) is_neg = true;
                parse(std::next(it_sparse), end, (*tr)[tree::right_child], is_neg);
                return;
            }
        }
    }
    // Traitement des + et - unaires :
    for (std::string s : {"+", "-"}) {
        auto it_sparse = find_statement(beg, end, {s, lexer::OPERATOR});
        if (it_sparse != end) { // Il y a bien un opérateur unaire de ce type
            {
                tr = std::make_shared<tree>(token{it_sparse->first, lexer::UNARY_OPERATOR}); // Père
                parse(std::next(it_sparse), end, (*tr)[tree::right_child]);
                return;
            }
        }
    }
    // Traitement des min et max :
    for ( std::string s : {"min", "max", "logical_and", "logical_or", "logical_xor"}) {
        auto it_minmax = find_statement(beg, end, {s, lexer::MINMAXOP});
        if (it_minmax != end) {
            tr = std::make_shared<tree>(*it_minmax);
            parse(std::next(it_minmax), end, tr);
            return;
        }
    }
    // Traitement des fonctions :
    auto it_fct = find_next_function(beg, end);
    if (it_fct != end) {
        tr = std::make_shared<tree>(*it_fct);
        parse(std::next(it_fct), end, (*tr)[tree::right_child]);
        return;
    }
    // Ne reste plus que les nombres, pi et les variables :
    assert(std::next(beg) == end); // Sinon, on a oublié un truc dans le sparser...
    if (beg->second == lexer::PI_VALUE)
    {
        tr = std::make_shared<tree>(token{std::string("3.1415926535897931"), lexer::NUMBER});
        return; 
    }
    tr = std::make_shared<tree>(*beg);
    if (beg->second == lexer::VARIABLE) m_symbol_table.push_back(beg->first);
    return;
}
}
//===============================================================================================================================