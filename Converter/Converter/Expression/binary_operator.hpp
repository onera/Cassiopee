#ifndef _CONVERTER_EXPRESSION_BINARY_OPERATOR_HPP_
#define _CONVERTER_EXPRESSION_BINARY_OPERATOR_HPP_
#include "Expression/constante.hpp"
#include "Expression/math_function.hpp"
#include "Expression/node.hpp"
#include "Expression/unary_operator.hpp"
#include "Expression/variable.hpp"
#include "Expression/logical_operator.hpp"
#include <algorithm>
#include <cmath>
#include <memory>

namespace Expression {
    /**
     * @brief Arithmetic binary operators
     * @details Arithmetics operators having left and right statements.
     *
     */
    class binary_operator : public ast::node {
      public:
        enum OP { ADD = 0, SUB, MULT, DIV, POW, ASSIGN, MIN, MAX };
        binary_operator(OP op, std::shared_ptr<ast::node> l_stmt, std::shared_ptr<ast::node> r_stmt)
            : m_left_statement(l_stmt), m_right_statement(r_stmt), m_op(op) {}

        double operator()(std::size_t i) const final {
            switch (m_op) {
            case ADD:
                return (*m_left_statement)(i) + (*m_right_statement)(i);
                break;
            case SUB:
                return (*m_left_statement)(i) - (*m_right_statement)(i);
                break;
            case MULT:
                return (*m_left_statement)(i) * (*m_right_statement)(i);
                break;
            case DIV:
                return (*m_left_statement)(i) / (*m_right_statement)(i);
                break;
            case POW:
                return std::pow((*m_left_statement)(i), (*m_right_statement)(i));
            case MIN:
                return std::min((*m_left_statement)(i), (*m_right_statement)(i));
            case MAX:
                return std::max((*m_left_statement)(i), (*m_right_statement)(i));
            case ASSIGN: {
                variable &var = dynamic_cast<variable &>(*m_left_statement);
                var[i]        = (*m_right_statement)(i);
                return var[i];
            }
            }
            return 0.;
        }

        simd_vector_wrapper eval_simd(std::size_t i) const final {
            switch (m_op) {
            case ADD: {
                return m_left_statement->eval_simd(i) + m_right_statement->eval_simd(i);
            } break;
            case SUB:
                return m_left_statement->eval_simd(i) - m_right_statement->eval_simd(i);
                break;
            case MULT:
                return m_left_statement->eval_simd(i) * m_right_statement->eval_simd(i);
                break;
            case DIV:
                return m_left_statement->eval_simd(i) / m_right_statement->eval_simd(i);
                break;
            case POW: {
                simd_vector_wrapper w;
                auto                u  = m_left_statement->eval_simd(i);
                auto                v  = m_right_statement->eval_simd(i);
                const double *      du = u.data(), *dv = v.data();
                double *            dw = w.data();
#pragma omp simd
                for (std::size_t ind = 0; ind < simd_vector_wrapper::max_size; ++ind)
                    dw[ind] = std::pow(du[ind], dv[ind]);
                return w;
            } break;
            case ASSIGN: {
                variable &var    = dynamic_cast<variable &>(*m_left_statement);
                auto      output = var.eval_simd(i);
                output.copy_in(m_right_statement->eval_simd(i));
                return output;
            }
            case MIN: {
                simd_vector_wrapper w;
                auto                u  = m_left_statement->eval_simd(i);
                auto                v  = m_right_statement->eval_simd(i);
                const double *      du = u.data(), *dv = v.data();
                double *            dw = w.data();
#pragma omp simd
                for (std::size_t ind = 0; ind < simd_vector_wrapper::max_size; ++ind)
                    dw[ind] = (du[ind] < dv[ind] ? du[ind] : dv[ind]);
                return w;
            }
            case MAX: {
                simd_vector_wrapper w;
                auto                u  = m_left_statement->eval_simd(i);
                auto                v  = m_right_statement->eval_simd(i);
                const double *      du = u.data(), *dv = v.data();
                double *            dw = w.data();
#pragma omp simd
                for (std::size_t ind = 0; ind < simd_vector_wrapper::max_size; ++ind)
                    dw[ind] = (du[ind] > dv[ind] ? du[ind] : dv[ind]);
                return w;
            }
            }
            return 0.;
        }

        explicit operator std::string() const final {
            std::string s;
            if ((m_op == MIN) or (m_op == MAX)) {
                s = std::string(m_op == MIN ? "min(" : "max(") + std::string(*m_left_statement) + ", ";
            } else {
                s = std::string("(") + std::string(*m_left_statement);
                switch (m_op) {
                case ADD:
                    s += "+";
                    break;
                case SUB:
                    s += "-";
                    break;
                case MULT:
                    s += "*";
                    break;
                case DIV:
                    s += "-";
                    break;
                case POW:
                    s += "**";
                    break;
                case ASSIGN:
                    s += "=";
                    break;
                }
            }
            return s + std::string(*m_right_statement) + ")";
        }

        virtual std::shared_ptr<node> derivate() final {
            switch (m_op) {
            case ASSIGN:
                return std::make_shared<binary_operator>(ASSIGN, m_left_statement->derivate(),
                                                         m_right_statement->derivate());
                break;
            case ADD:
                return std::make_shared<binary_operator>(ADD, m_left_statement->derivate(),
                                                         m_right_statement->derivate());
                break;
            case SUB:
                return std::make_shared<binary_operator>(SUB, m_left_statement->derivate(),
                                                         m_right_statement->derivate());
                break;
            case MULT:
                return std::make_shared<binary_operator>(
                    ADD, std::make_shared<binary_operator>(MULT, m_left_statement->derivate(), m_right_statement),
                    std::make_shared<binary_operator>(MULT, m_left_statement, m_right_statement->derivate()));
                break;
            case DIV:
                return std::make_shared<binary_operator>(
                    DIV,
                    std::make_shared<binary_operator>(
                        SUB, std::make_shared<binary_operator>(MULT, m_left_statement->derivate(), m_right_statement),
                        std::make_shared<binary_operator>(MULT, m_left_statement, m_right_statement->derivate())),
                    std::make_shared<binary_operator>(MULT, m_right_statement, m_right_statement));
                break;
            case MIN:
            	return std::make_shared<logical_operator>(logical_operator::LE, m_left_statement, m_right_statement);
            case MAX:
            	return std::make_shared<logical_operator>(logical_operator::GE, m_left_statement, m_right_statement);
            case POW: {
                auto &r             = *m_right_statement;
                bool  r_stmt_is_cte = typeid(r).hash_code() == typeid(constante).hash_code();
                if (r_stmt_is_cte) { // d(x^a) = a * dx * x^{a-1}
                    double value = r(0);
                    return std::make_shared<binary_operator>(
                        MULT,
                        std::make_shared<binary_operator>(MULT, std::make_shared<constante>(value),
                                                          m_left_statement->derivate()),
                        std::make_shared<binary_operator>(POW, m_left_statement,
                                                          std::make_shared<constante>(value - 1)));
                }
                // d(x^y) = (dy * ln(x) + y * dx / x ) * x^y
                return std::make_shared<binary_operator>(
                    MULT,
                    std::make_shared<binary_operator>(
                        ADD,
                        std::make_shared<binary_operator>(MULT, m_right_statement->derivate(),
                                                          std::make_shared<ln>(m_left_statement)),
                        std::make_shared<binary_operator>(
                            DIV,
                            std::make_shared<binary_operator>(MULT, m_right_statement, m_left_statement->derivate()),
                            m_left_statement)),
                    std::make_shared<binary_operator>(POW, m_left_statement, m_right_statement));
            } break;
            }
            return nullptr;
        }

        std::shared_ptr<node> optimize() final {
            // Cas ou les deux statements gauche et droite sont les memes => simplification
            // A faire avant simplification, car sinon on perd facilement l'égalité par adresse
            if (m_left_statement == m_right_statement) {
                if (m_op == ADD) { // x+x -> 2 * x
                    m_left_statement = m_left_statement->optimize();
                    return std::make_shared<binary_operator>(MULT, std::make_shared<constante>(2.), m_left_statement);
                }
                if (m_op == SUB) { // x-x -> 0
                    return std::make_shared<constante>(0.);
                }
                if (m_op == DIV) { // x / x -> 1
                    return std::make_shared<constante>(1.);
                }
                if (m_op == MIN) { // min(x,x) -> x
                    return m_left_statement;
                }
                if (m_op == MAX) { // max(x,x) -> x
                    return m_left_statement;
                }
            }

            m_left_statement    = m_left_statement->optimize();
            m_right_statement   = m_right_statement->optimize();
            auto &l             = *m_left_statement;
            auto &r             = *m_right_statement;
            bool  l_stmt_is_cte = typeid(l).hash_code() == typeid(constante).hash_code();
            bool  r_stmt_is_cte = typeid(r).hash_code() == typeid(constante).hash_code();
            // Si left stmt et right stmt sont deux constantes -> on evalue
            if (l_stmt_is_cte && r_stmt_is_cte) return std::make_shared<constante>((*this)(0));
            // Si left stmt est une constante :
            if (l_stmt_is_cte) {
                double value = l(0);
                // 0 + x -> x
                if ((m_op == ADD) && (value == 0)) return m_right_statement;
                // 0 - x -> -x
                if ((m_op == SUB) && (value == 0))
                    return std::make_shared<unary_operator>(unary_operator::SUB, m_right_statement);
                // 0 * x -> 0
                if ((m_op == MULT) && (value == 0)) return std::make_shared<constante>(0.);
                // 1 * x -> x
                if ((m_op == MULT) && (value == 1)) return m_right_statement;
                // 0 / x -> 0
                if ((m_op == DIV) && (value == 0)) return std::make_shared<constante>(0.);
                // 1 ** x -> 1
                if ((m_op == POW) && (value == 1)) return std::make_shared<constante>(1.);
            }
            if (r_stmt_is_cte) {
                double value = r(0);
                // x + 0 -> x
                if ((m_op == ADD) && (value == 0)) return m_left_statement;
                // x - 0 -> x
                if ((m_op == SUB) && (value == 0)) return m_left_statement;
                // x * 0 -> 0
                if ((m_op == MULT) && (value == 0)) return std::make_shared<constante>(0.);
                // x * 1 -> x
                if ((m_op == MULT) && (value == 1)) return m_left_statement;
                // x / 1 -> x
                if ((m_op == DIV) && (value == 1)) return m_left_statement;
                // x**1 -> x
                if ((m_op == POW) && (value == 1)) return m_left_statement;
                // x**(-1) -> 1/x
                if ((m_op == POW) && (value == -1))
                    return std::make_shared<binary_operator>(DIV, std::make_shared<constante>(1.), m_left_statement);
            }
            return std::make_shared<binary_operator>(m_op, m_left_statement, m_right_statement);
        }

      private:
        std::shared_ptr<ast::node> m_left_statement, m_right_statement;
        OP                         m_op;
    };

} // namespace Expression
#endif
