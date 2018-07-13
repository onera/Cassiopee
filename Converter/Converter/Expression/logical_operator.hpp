#ifndef _CONVERTER_EXPRESSION_LOGICAL_OPERATOR_HPP_
#define _CONVERTER_EXPRESSION_LOGICAL_OPERATOR_HPP_
#include "Expression/node.hpp"
#include "Expression/constante.hpp"

namespace Expression
{
    class logical_operator : public ast::node
    {
    public:
        enum OP {
            EQU = 0, NE, LT, LE, GT, GE, AND, OR, XOR, NOT
        };

        logical_operator( OP op, std::shared_ptr<ast::node> l_stmt, std::shared_ptr<ast::node> r_stmt ) :
            m_left_statement(l_stmt), m_right_statement(r_stmt), m_op(op)
        {}

        double operator() ( std::size_t i ) const final
        {
            switch(m_op)
            {
            case EQU:
                return double((*m_left_statement)(i) == (*m_right_statement)(i));
                break;
            case NE:
                return double((*m_left_statement)(i) != (*m_right_statement)(i));
                break;
            case LT:
                return double((*m_left_statement)(i) < (*m_right_statement)(i));
                break;
            case LE:
                return double((*m_left_statement)(i) <= (*m_right_statement)(i));
                break;          
            case GT:
                return double((*m_left_statement)(i) > (*m_right_statement)(i));
                break;
            case GE:
                return double((*m_left_statement)(i) >= (*m_right_statement)(i));
                break;
            case AND:
                return double(long((*m_left_statement)(i)) && long((*m_right_statement)(i))); 
                break;
            case OR:
                return double(long((*m_left_statement)(i)) || long((*m_right_statement)(i))); 
                break;
            case XOR:
                {
                    bool l = (*m_left_statement)(i)  != 0;
                    bool r = (*m_right_statement)(i) != 0;
                    return double(l != r);
                }
                break;
            }
            return 0.;
        }

        simd_vector_wrapper eval_simd ( std::size_t i ) const final
        {
            switch(m_op)
            {
            case EQU:
                return m_left_statement->eval_simd(i) == m_right_statement->eval_simd(i);
                break;
            case NE:
                return m_left_statement->eval_simd(i) != m_right_statement->eval_simd(i);
                break;
            case LT:
            	return m_left_statement->eval_simd(i) < m_right_statement->eval_simd(i);
                break;
            case LE:
            	return m_left_statement->eval_simd(i) <= m_right_statement->eval_simd(i);
                break;          
            case GT:
                return m_left_statement->eval_simd(i) > m_right_statement->eval_simd(i);
                break;
            case GE:
                return m_left_statement->eval_simd(i) >= m_right_statement->eval_simd(i);
                break;
            case AND:
                {
                    simd_vector_wrapper w;
                    auto                u  = m_left_statement->eval_simd(i);
                    auto                v  = m_right_statement->eval_simd(i);
                    const double *      du = u.data(), *dv = v.data();
                    double *            dw = w.data();

#                   pragma omp simd
                    for (std::size_t ind = 0; ind < simd_vector_wrapper::max_size; ++ind)
                        dw[ind] = double( (du[ind]!=0) and (dv[ind]!=0) );
                    return w;
                }
                break;
            case OR:
                {
                    simd_vector_wrapper w;
                    auto                u  = m_left_statement->eval_simd(i);
                    auto                v  = m_right_statement->eval_simd(i);
                    const double *      du = u.data(), *dv = v.data();
                    double *            dw = w.data();

#                   pragma omp simd
                    for (std::size_t ind = 0; ind < simd_vector_wrapper::max_size; ++ind)
                        dw[ind] = double( (du[ind]!=0) or (dv[ind]!=0) );
                    return w;
                }
                break;
            case XOR:
                {
                    simd_vector_wrapper w;
                    auto                u  = m_left_statement->eval_simd(i);
                    auto                v  = m_right_statement->eval_simd(i);
                    const double *      du = u.data(), *dv = v.data();
                    double *            dw = w.data();

#                   pragma omp simd
                    for (std::size_t ind = 0; ind < simd_vector_wrapper::max_size; ++ind)
                        dw[ind] = double( (du[ind]!=0) != (dv[ind]!=0) );
                    return w;
                }
                break;
            }
            return simd_vector_wrapper();
        }


        explicit operator std::string() const final
        {
            std::string s = std::string("[") + std::string(*m_left_statement);
            switch(m_op)
            {
            case EQU:
                s += "==";
                break;
            case NE:
                s += "!=";
                break;
            case LT:
                s += "<";
                break;
            case LE:
                s += "<=";
                break;
            case GT:
                s += ">";
                break;
            case GE:
                s += ">=";
                break;
            case AND:
                s += " and ";
                break;
            case OR:
                s += " or ";
                break;
            case XOR:
                s += " xor ";
                break;
            }
            return s + std::string(*m_right_statement) + "]";
        }

        virtual std::shared_ptr<node> derivate() final
        {
            return std::make_shared<constante>(0.);
        }

        std::shared_ptr<node> optimize() final
        {
            m_left_statement = m_left_statement->optimize();
            m_right_statement = m_right_statement->optimize();
            auto& l = *m_left_statement;
            auto& r = *m_right_statement;
            bool l_stmt_is_cte = typeid(l).hash_code() == typeid(constante).hash_code();
            bool r_stmt_is_cte = typeid(r).hash_code() == typeid(constante).hash_code();
            // Si left stmt et right stmt sont deux constantes -> on evalue
            if ( l_stmt_is_cte && r_stmt_is_cte ) return std::make_shared<constante>((*this)(0));
            return std::make_shared<logical_operator>(m_op, m_left_statement, m_right_statement);
        }
    private:
        std::shared_ptr<ast::node> m_left_statement, m_right_statement;     
        OP m_op;
    };

}

#endif