#include "Expression/math_function.hpp"
#include "Expression/binary_operator.hpp"
#include "Expression/constante.hpp"
#include "Expression/unary_operator.hpp"

namespace Expression
{
std::shared_ptr<ast::node> cos::df( )
{
    return std::make_shared<unary_operator>(unary_operator::SUB,
                                            std::make_shared<::Expression::sin>(m_right_statement));
}
// .................................................................................................
std::shared_ptr<ast::node> cos::simplify( )
{
    // Faudra faire le cas du acos avec le cos quand on aura rajouter acos comme fonction...
    // Sinon, ici, a part le coup de la constante deja traitee par fonction, bof...
    // Idem, simplification a faire avec tan si besoin
    return std::make_shared<::Expression::cos>(m_right_statement);
}
// -------------------------------------------------------------------------------------------------
std::shared_ptr<ast::node> sin::df( ) { return std::make_shared<::Expression::cos>(m_right_statement); }
// .................................................................................................
std::shared_ptr<ast::node> sin::simplify( )
{
    // Faudra faire le cas du asin avec le sin quand on aura rajouter acos comme fonction...
    // Sinon, ici, a part le coup de la constante deja traite par fonction, bof...
    // Besoin simplification avec sin / tan ?
    return std::make_shared<::Expression::sin>(m_right_statement);
}
// -------------------------------------------------------------------------------------------------
std::shared_ptr<ast::node> tan::df( ) { 
    return std::make_shared<binary_operator>(binary_operator::DIV, std::make_shared<constante>(1.),
                                                                   std::make_shared<binary_operator>(binary_operator::POW,
                            std::make_shared<::Expression::cos>(m_right_statement), std::make_shared<constante>(2.))); 
}
// .................................................................................................
std::shared_ptr<ast::node> tan::simplify( )
{
    // Faudra faire le cas du asin avec le sin quand on aura rajouter acos comme fonction...
    // Sinon, ici, a part le coup de la constante deja traite par fonction, bof...
    // Besoin simplification avec sin / tan ?
    return std::make_shared<::Expression::tan>(m_right_statement);
}
// -------------------------------------------------------------------------------------------------
std::shared_ptr<ast::node> ln::df( )
{
    return std::make_shared<binary_operator>(binary_operator::DIV, std::make_shared<constante>(1.), m_right_statement);
}
// .................................................................................................
std::shared_ptr<ast::node> ln::simplify( ) {
    auto& r = *m_right_statement;// Deja simplifie dans function...
    if (typeid(r).hash_code() == typeid(::Expression::exp).hash_code() ) {
        return m_right_statement;
    }
    return std::make_shared<::Expression::ln>(m_right_statement); 
}
// -------------------------------------------------------------------------------------------------
std::shared_ptr<ast::node> exp::simplify( ) {
    auto& r = *m_right_statement;// Deja simplifie dans function...
    if (typeid(r).hash_code() == typeid(::Expression::ln).hash_code() ) {
        return m_right_statement;
    }
    return std::make_shared<::Expression::exp>(m_right_statement); 
}
// -------------------------------------------------------------------------------------------------
std::shared_ptr<ast::node> sqrt::df( )
{
    return std::make_shared<binary_operator>(binary_operator::DIV, std::make_shared<constante>(0.5),std::make_shared<::Expression::sqrt>(m_right_statement));
}
// .................................................................................................
std::shared_ptr<ast::node> sqrt::simplify( ) {
    // Simplification avec puissance au carré : donne la valeur absolue ( pas encore mise ) ?
    return std::make_shared<::Expression::sqrt>(m_right_statement);
}
// -------------------------------------------------------------------------------------------------
std::shared_ptr<ast::node> cosh::df( )
{
    return std::make_shared<sinh>(m_right_statement);
}
// .................................................................................................
std::shared_ptr<ast::node> cosh::simplify( ) {
    // Simplification avec puissance au carré : donne la valeur absolue ( pas encore mise ) ?
    return std::make_shared<::Expression::cosh>(m_right_statement);
}
// -------------------------------------------------------------------------------------------------
std::shared_ptr<ast::node> sinh::df( )
{
    return std::make_shared<cosh>(m_right_statement);
}
// .................................................................................................
std::shared_ptr<ast::node> sinh::simplify( ) {
    // Simplification avec puissance au carré : donne la valeur absolue ( pas encore mise ) ?
    return std::make_shared<::Expression::sinh>(m_right_statement);
}
// -------------------------------------------------------------------------------------------------
std::shared_ptr<ast::node> abs::df( )
{
    return std::make_shared<sign>(m_right_statement);
    //return std::make_shared<constante>(0.);
}
// .................................................................................................
std::shared_ptr<ast::node> abs::simplify( ) {
    return std::make_shared<::Expression::abs>(m_right_statement);
}
// -------------------------------------------------------------------------------------------------
std::shared_ptr<ast::node> sign::df( )
{
    return std::make_shared<constante>(0.);
}
// .................................................................................................
std::shared_ptr<ast::node> sign::simplify( ) {
    return std::make_shared<::Expression::sign>(m_right_statement);
}
// #################################################################################################
void init_math_functions( )
{
    function::list_of_functions["cos"] = std::make_shared<cos>(std::make_shared<constante>(0.));
    function::list_of_functions["sin"] = std::make_shared<sin>(std::make_shared<constante>(0.));
    function::list_of_functions["tan"] = std::make_shared<tan>(std::make_shared<constante>(0.));
    function::list_of_functions["ln"]  = std::make_shared<ln>(std::make_shared<constante>(0.));
    function::list_of_functions["log"]  = std::make_shared<ln>(std::make_shared<constante>(0.));
    function::list_of_functions["exp"] = std::make_shared<exp>(std::make_shared<constante>(0.));
    function::list_of_functions["sqrt"] = std::make_shared<sqrt>(std::make_shared<constante>(0.));
    function::list_of_functions["cosh"] = std::make_shared<cosh>(std::make_shared<constante>(0.));
    function::list_of_functions["sinh"] = std::make_shared<sinh>(std::make_shared<constante>(0.));
    function::list_of_functions["abs"] = std::make_shared<abs>(std::make_shared<constante>(0.));
    function::list_of_functions["sign"] = std::make_shared<sign>(std::make_shared<constante>(0.));
}
}