#ifndef _CONVERTER_EXPRESSION_VARIABLE_HPP_
#define _CONVERTER_EXPRESSION_VARIABLE_HPP_
#include <string>
#include "Expression/node.hpp"
#include "Expression/symbol_table.hpp"
#include "Expression/simd_vector_wrapper.hpp"

namespace Expression
{
    class variable : public ast::node
    {
    public:
        variable( const std::string& name ) : m_name(name), m_data(symbol_table::get()[name]) {}

        double operator() ( std::size_t i ) const final
        {
            return m_data.size==1 ? m_data.data[0] : m_data.data[i];
        }

        simd_vector_wrapper eval_simd( std::size_t i ) const final
        {
            if (m_data.size == 1)
                return { m_data.data[0] };
            return { m_data.data+i };
        }        

        double& operator[] ( std::size_t i ) const 
        {
            return m_data.size==1 ? m_data.data[0] : m_data.data[i];
        }

        explicit operator std::string() const final
        {
            return m_name;
        }

        std::shared_ptr<node> derivate( ) final
        {
            std::string d_name = "d_"+m_name;
            symbol_table::get().insert(d_name, symbol_table::data_t{});
            return std::make_shared<variable>(d_name);
        }

        std::shared_ptr<node> optimize() final
        {
            return std::make_shared<variable>(m_name);
        }

    private:
        std::string m_name;
        symbol_table::data_t& m_data;
    };
}

#endif