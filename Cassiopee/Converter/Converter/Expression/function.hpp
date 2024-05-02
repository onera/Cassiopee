#ifndef _CONVERTER_EXPRESSION_FUNCTION_HPP_
#define _CONVERTER_EXPRESSION_FUNCTION_HPP_
#include "Expression/node.hpp"
#include "Expression/simd_vector_wrapper.hpp"
#include <memory>
#include <unordered_map>

namespace Expression {
    class function : public ast::node {
      public:
        function(std::shared_ptr<node> &r_stmt) : m_right_statement(r_stmt) {}
        double operator()(std::size_t i) const final {
            return apply_to((*m_right_statement)(i));
        }
        simd_vector_wrapper eval_simd(std::size_t i) const final {
            return apply_to(m_right_statement->eval_simd(i));
        }

        explicit operator std::string() const final {
            return name() + "(" + std::string(*m_right_statement) + ")";
        }

        virtual std::shared_ptr<node> derivate() final;

        /**
         * @brief      Clone a function but modify the argument ( right
         * statement ) of the function
         *
         * @param      right_stmt  The right statement
         *
         * @return     Copy of this object.
         */
        virtual std::shared_ptr<function>
        clone(std::shared_ptr<ast::node> &right_stmt) const = 0;

        std::shared_ptr<node> optimize() final;

        static std::unordered_map<std::string, std::shared_ptr<function>>
            list_of_functions;

      protected:
        virtual std::string                name() const = 0;
        virtual std::shared_ptr<ast::node> df()         = 0;
        virtual std::shared_ptr<node>      simplify()   = 0;
        virtual simd_vector_wrapper
                       apply_to(const simd_vector_wrapper &x) const = 0;
        virtual double apply_to(double x) const                     = 0;
        // Pas de statement à gauche pour une fonction
        std::shared_ptr<ast::node> m_right_statement;
    };

    std::shared_ptr<function>
    make_function(const std::string &         name,
                  std::shared_ptr<ast::node> &right_stmt);
} // namespace Expression

#endif
