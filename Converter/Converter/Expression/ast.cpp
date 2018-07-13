#include "Expression/ast.hpp"
#include "Expression/binary_operator.hpp"
#include "Expression/constante.hpp"
#include "Expression/function.hpp"
#include "Expression/logical_operator.hpp"
#include "Expression/node.hpp"
#include "Expression/parser.hpp"
#include "Expression/simd_vector_wrapper.hpp"
#include "Expression/unary_operator.hpp"
#include <algorithm>

// CONTINUER ICI POUR NETTOYAGE
namespace Expression {
    class ast::tree {
      public:
        tree(symbol_table_t &table, parser &p) { m_root = build_tree(table, *p.getTree()); }
        tree(std::shared_ptr<ast::node> root) : m_root(root) {}
        ast::node &root() { return *m_root; }

        void optimize() { m_root = m_root->optimize(); }

        std::string to_string() const { return std::string(*m_root); }

      private:
        std::shared_ptr<ast::node> build_tree(symbol_table_t &table, parser::tree &t);
        std::shared_ptr<ast::node> m_root;
        friend ast;
    };

    std::shared_ptr<ast::node> ast::tree::build_tree(symbol_table_t &table, parser::tree &t) {
        parser::token tok = *t;
        switch (tok.second) {
        case lexer::FUNCTION: {
            std::transform(tok.first.begin(), tok.first.end(), tok.first.begin(), ::tolower);
            try {
                auto f               = function::list_of_functions.at(tok.first);
                auto right_statement = build_tree(table, *t[parser::tree::right_child]);
                return f->clone(right_statement);
            } catch(std::out_of_range& obj)
            {
                std::string errTxt = std::string("Unknown function ") + tok.first + " provided... Existing";
                throw std::out_of_range(errTxt);
            }
            return nullptr;
        } break;
        case lexer::NUMBER: {
            return std::make_shared<constante>(std::stod(tok.first));
        } break;
        case lexer::OPERATOR: // Binary operator for that :
        case lexer::MINMAXOP:
        {
            auto left_statement  = build_tree(table, *t[parser::tree::left_child]);
            auto right_statement = build_tree(table, *t[parser::tree::right_child]);
            // Analyse de l'opérateur en question :
            bool                 is_logic_op = false;
            binary_operator::OP  aop;
            logical_operator::OP lop;
            if (tok.first == "=")   aop = binary_operator::ASSIGN;
            if (tok.first == "+")   aop = binary_operator::ADD;
            if (tok.first == "-")   aop = binary_operator::SUB;
            if (tok.first == "*")   aop = binary_operator::MULT;
            if (tok.first == "/")   aop = binary_operator::DIV;
            if (tok.first == "**")  aop = binary_operator::POW;
            if (tok.first == "min") aop = binary_operator::MIN;
            if (tok.first == "max") aop = binary_operator::MAX;
            if (tok.first == "logical_and") {
                lop = logical_operator::AND;
                is_logic_op = true;
            }
            if (tok.first == "logical_or") {
                lop = logical_operator::OR;
                is_logic_op = true;
            }
            if (tok.first == "logical_xor") {
                lop = logical_operator::XOR;
                is_logic_op = true;
            }

            if (tok.first == "==") {
                lop         = logical_operator::EQU;
                is_logic_op = true;
            }
            if (tok.first == "!=") {
                lop         = logical_operator::NE;
                is_logic_op = true;
            }
            if (tok.first == "<") {
                lop         = logical_operator::LT;
                is_logic_op = true;
            }
            if (tok.first == "<=") {
                lop         = logical_operator::LE;
                is_logic_op = true;
            }
            if (tok.first == ">") {
                lop         = logical_operator::GT;
                is_logic_op = true;
            }
            if (tok.first == ">=") {
                lop         = logical_operator::GE;
                is_logic_op = true;
            }
            if (not is_logic_op)
                return std::make_shared<binary_operator>(aop, left_statement, right_statement);
            else
                return std::make_shared<logical_operator>(lop, left_statement, right_statement);
        } break;
        case lexer::UNARY_OPERATOR: {
            auto               right_statement = build_tree(table, *t[parser::tree::right_child]);
            unary_operator::OP uop;
            if (tok.first == "+") uop = unary_operator::ADD;
            if (tok.first == "-") uop = unary_operator::SUB;
            return std::make_shared<unary_operator>(uop, right_statement);
        } break;
        case lexer::VARIABLE:
            return std::make_shared<variable>(tok.first);
            break;
        default:
            // Ameliorer le message d'erreur pour mieux aider l'utilisateur
            throw std::logic_error("Unexcepted statement ! Error in your expression...");
        }
        return nullptr;
    }
    // =====================================================================================================
    ast::ast(const std::string &expr) : m_pt_tree() {
        // On passe le lexer et le parser
        lexer lx;
        lx.analyze(expr);
        parser parse(lx);
        auto & st = symbol_table::get();
        // On construit la table des symboles :
        for (auto &var_name : parse.get_symbol_table())
            st[var_name] = symbol_table::data_t();
        // Puis l'arbre d'expression abstraite :
        m_pt_tree = std::make_shared<tree>(st, parse);
        m_pt_tree->optimize();
    }

    ast::ast(std::shared_ptr<tree> tr) : m_pt_tree(tr) {}

    std::vector<double> ast::operator()(const std::unordered_map<std::string, vector_view<double>> &params) {
        auto &st = symbol_table::get();
        for (auto &v : st)
            v.second = symbol_table::data_t();
        std::size_t sz = 1;
        for (auto v : params) {
            auto &d = st.at(v.first);
            d.size  = v.second.size();
            if (d.size > 1) {
                if (sz == 1) sz = d.size;
                if (sz != d.size) throw std::invalid_argument("Uncompatible dimensions between vectors");
            }
            d.data = v.second.data();
        }
        auto &              root = m_pt_tree->root();
        std::vector<double> result(sz);
            if (sz > simd_vector_wrapper::max_size) {
#pragma omp parallel for 
                for (std::size_t i = 0; i < sz - simd_vector_wrapper::max_size; i += simd_vector_wrapper::max_size) {
                    simd_vector_wrapper(result.data() + i) = root.eval_simd(i);
                }
            }
#pragma omp parallel for
        for (std::size_t i = (sz/simd_vector_wrapper::max_size)*simd_vector_wrapper::max_size; i < sz; ++i) {
                result[i] = root(i);
            }
        return result;
    }

    void ast::eval(const std::unordered_map<std::string, vector_view<double>> &params, vector_view<double> &result) {
        auto &st = symbol_table::get();
        for (auto &v : st)
            v.second = symbol_table::data_t();
        std::size_t sz = 0;
        for (auto v : params) {
            auto &d = st.at(v.first);
            d.size  = v.second.size();
            if (d.size > 1) {
                if (sz == 0) sz = d.size;
                if (sz != d.size) throw std::invalid_argument("Uncompatible dimensions between vectors");
            }
            d.data = v.second.data();
        }
        auto &root = m_pt_tree->root();
        if (sz > simd_vector_wrapper::max_size) {
            if (sz > 2*simd_vector_wrapper::max_size)
#               pragma omp parallel for if (sz > 8*simd_vector_wrapper::max_size)
                for (std::size_t i = 0; i < sz/simd_vector_wrapper::max_size; i ++) {
                    simd_vector_wrapper(result.data() + i*simd_vector_wrapper::max_size) = root.eval_simd(i*simd_vector_wrapper::max_size);
                }
        }
#       pragma omp parallel for
        for (std::size_t i = (sz/simd_vector_wrapper::max_size)*simd_vector_wrapper::max_size; i < sz; ++i) {
            result[i] = root(i);
        }
    }

    void ast::eval(const std::unordered_map<std::string, vector_view<double>> &params) {
        auto &st = symbol_table::get();
        for (auto &v : st)
            v.second = symbol_table::data_t();
        std::size_t sz = 0;
        for (auto v : params) {
            auto &d = st.at(v.first);
            d.size  = v.second.size();
            if (d.size > 1) {
                if (sz == 0) sz = d.size;
                if (sz != d.size) throw std::invalid_argument("Uncompatible dimensions between vectors");
            }
            d.data = v.second.data();
        }
        auto &root = m_pt_tree->root();
            if (sz > simd_vector_wrapper::max_size) {
#pragma omp parallel for if (sz > 8*simd_vector_wrapper::max_size)
                for (std::size_t i = 0; i < sz/simd_vector_wrapper::max_size; i++) {
                    root.eval_simd(i*simd_vector_wrapper::max_size);
                }
            }
#pragma omp parallel for
        for (std::size_t i = (sz/simd_vector_wrapper::max_size)*simd_vector_wrapper::max_size; i < sz; ++i) {
                root(i);
            }
    }

    ast ast::derivate() {
        ast a(std::make_shared<tree>(m_pt_tree->root().derivate()));
        a.m_pt_tree->optimize();
        return a;
    }

    std::string ast::to_string() const { return m_pt_tree->to_string(); }

} // namespace Expression
