# include "Expression/function.hpp"
# include "Expression/binary_operator.hpp"
# include "Expression/constante.hpp"

namespace Expression
{
    std::unordered_map<std::string, std::shared_ptr<function>> function::list_of_functions;

    std::shared_ptr<ast::node> function::derivate( )
    {
        return std::make_shared<binary_operator>(binary_operator::MULT, m_right_statement->derivate( ), df( ));
    }
    // ---------------------------------------------------------------------------------------------------------
    std::shared_ptr<ast::node> function::optimize()
    {
        m_right_statement = m_right_statement->optimize();
        auto& r = *m_right_statement;
        if (typeid(r).hash_code() == typeid(constante).hash_code() ) {
            return std::make_shared<constante>((*this)(r(0)));
        }
        return simplify();
    }
}