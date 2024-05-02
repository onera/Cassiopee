#ifndef _CONVERTER_EXPRESSION_HPP_
#define _CONVERTER_EXPRESSION_HPP_
#include <unordered_map>
#include <string>
#include <memory>
#include <stdexcept>

namespace Expression
{
	/**
	 * @brief Symbol table of variables for expressions.
	 * @details The symbol's table manage the variables used in expressions.
     * It's a singleton's pattern and so, all expressions shared the same
     * symbol's table. Modifying a data of a variable for an expression modify the same data
     * for all other expressions, so avoid to use expressions in a  multithreading context
     * ( in a distributed computation, all will be fine ). 
     * The table of symbol doesn't owned the data of the variables.
     * Usage :
     * -------
     *     auto& sym_tab = Expression::symbol_table::get();
     *     sym_tab["x"]    = symbol_table::data_t{100, pointer_on_array};
     *     sym_tab["scal"] = symbol_table::data_t{1, &scalar};// ok, it's a scalar
     *     double* values = sym_tab.at("x").data(); // OK, x exists
     *     double* values2= sym_tab.at("y").data(); // out_of_range exception, y doesn't exists !
	 * 
	 */
	class symbol_table
	{
	public:
	    struct data_t {
    	    std::size_t size;
        	double *    data;
        	data_t( ) : size(0), data(nullptr)
        	{}
    	};
		using dictionnary = std::unordered_map<std::string, data_t>;
		using iterator    = dictionnary::iterator;
		using const_iterator = dictionnary::const_iterator;

        /**
         * @brief Create an empty table of symbol
         */
        symbol_table() = default;
    	symbol_table( const symbol_table& ) = delete;
    	symbol_table( symbol_table&& )      = delete;
        ~symbol_table() = default;

        /**
         * @brief Return the number of variables
         * @details Return the number of variables, scalar and vectors !
         * @return The size of the symbol's table
         */
    	std::size_t size() const { return m_data.size(); }
        /**
         * @brief Return a variable
         * @details Return a variable descriptor for the variable with name var.
         *          If the variable doesn't exist, throw a invalid_argument exception
         *          The variable is only readable.
         * 
         * @param var The name of the variable
         * @return The data description of the variable ( size and adress )
         */
    	data_t data(const std::string &var) const
    	{
        	auto search = m_data.find(var);
        	if (search == m_data.end( )) {
            	throw std::invalid_argument("These parameter doesn't exist");
        	}
        	return search->second;
    	}

        /**
         * @brief Return a variable
         * @details Return a variable descriptor for the variable with name var
         *          If the variable doesn't exist, throw a invalid_argument exception
         *          The variable is in writable and readable.
         * 
         * @param var The name of the variable
         * @return The data description of the variable ( size and adress )
         */
    	data_t &data(const std::string &var)
    	{
        	auto search = m_data.find(var);
        	if (search == m_data.end( )) {
            	throw std::invalid_argument("These parameter doesn't exist");
        	}
        	return search->second;
    	}

        /**
         * @brief Return the variable with name var
         * @details Return the variable with name var.
         *          If this variable doesn't exist, create a new variable
         *          in the symbol's table
         * 
         * @param var [description]
         */
    	data_t& operator[] ( const std::string &var )
    	{
    		return m_data[var];
    	}

        /**
         * @brief Return the data for the variable with name var
         * @details Return the data for the variable with name var.
         *          If the variable doesn't exist, return an out_of_range exception
         *          The data is readable and writable
         * 
         * @param var The name of the variable
         * @return The data description of the variable ( size and address )
         */
    	data_t& at ( const std::string& var )
    	{
    		return m_data.at(var);
    	}

        /**
         * @brief Return the data for the variable with name var
         * @details Return the data for the variable with name var.
         *          If the variable doesn't exist, return an out_of_range exception
         *          The data is only readable
         * 
         * @param var The name of the variable
         * @return The data description of the variable ( size and address )
         */
    	const data_t& at ( const std::string& var ) const 
    	{
    		return m_data.at(var);
    	}

        /**
         * @brief Return an iterator on the beginning of the symbol's table
         * @details Return an iterator on the beginning of the symbol's table)
         * @return An iterator
         */
    	iterator begin() { return m_data.begin(); }
        /**
         * @brief Return a const_iterator on the beginning of the symbol's table
         * @details Return a const_iterator on the beginning of the symbol's table)
         * @return A const_iterator
         */
    	const_iterator cbegin() { return m_data.cbegin(); }
        /**
         * @brief Return a const_iterator on the beginning of the symbol's table
         * @details Return a const_iterator on the beginning of the symbol's table)
         * @return A const_iterator
         */
    	const_iterator begin() const { return m_data.cbegin(); }

        /**
         * @brief Return an iterator at the end of the symbol's table
         * @details Return an iterator at the end of the symbol's table)
         * @return An iterator
         */
    	iterator end() { return m_data.end(); }
        /**
         * @brief Return a const_iterator at the end of the symbol's table
         * @details Return a const_iterator at the end of the symbol's table)
         * @return A const_iterator
         */
    	const_iterator cend() { return m_data.cend(); }
        /**
         * @brief Return a const_iterator at the end of the symbol's table
         * @details Return a const_iterator at the end of the symbol's table)
         * @return A const_iterator
         */
    	const_iterator end() const { return m_data.cend(); }

    	void insert( const std::string& name, const data_t& d )
    	{
    		m_data.insert({name,d});
    	}

        /**
         * @brief Return the symbol's table
         * @details Return the unique symbol's table shared by all expressions
         * @return The singleton symbol's table
         */
    	static symbol_table& get();

	private:
		dictionnary m_data;
	};
}

#endif
