#include <algorithm>
#include <cassert>
#include <cctype>
#include <iostream>
#if defined(USE_C_REGEX)
#include <new>
#include <regex.h>
#else
#include <regex>
#endif
#include <sstream>
#include <stdexcept>

#include "Expression/lexer.hpp"

namespace Expression {
    void ast::lexer::analyze(const std::string &expr) {
        m_tokens.reserve(expr.size()) ;
        std::size_t index = 0;
#if defined(USE_C_REGEX)
        int     error;
        regex_t func_regex, var_regex, number_regex;
        size_t nmatch_func_regex, nmatch_var_regex, nmatch_number_regex;
        error = regcomp(&func_regex, "(\\w+ *\\()", REG_EXTENDED);
        if (error) {
            char * error_msg;
            size_t err_sz      = regerror(error, &func_regex, NULL, 0);
            error_msg          = new char[err_sz];
            err_sz = regerror(error, &func_regex, error_msg, err_sz);
            std::string msg    = std::string("Failed to compile regular expression : ") + std::string(error_msg);
            delete[] error_msg;
            throw std::runtime_error(msg);
        }
        nmatch_func_regex      = func_regex.re_nsub;
        error = regcomp(&var_regex, "(\\{[_a-zA-Z0-9]+\\})", REG_EXTENDED);
        if (error) {
            char * error_msg;
            size_t err_sz      = regerror(error, &func_regex, NULL, 0);
            error_msg          = new char[err_sz];
            err_sz = regerror(error, &func_regex, error_msg, err_sz);
            std::string msg    = std::string("Failed to compile regular expression : ") + std::string(error_msg);
            delete[] error_msg;
            throw std::runtime_error(msg);
        }
        nmatch_var_regex       = var_regex.re_nsub;
        error = regcomp(&number_regex, R"RAW(^[-+]?[0-9]*\.?[0-9]*([eE][-+]?[0-9]+)?)RAW",REG_EXTENDED);
            //[:digit:]+(?:[eE][+-]?[:digit:].)RAW", REG_EXTENDED);
        //error = regcomp(&number_regex, R"RAW(([\d.]+(?:[eE][+-]?\d+)?))RAW", REG_EXTENDED);
        if (error) {
            char * error_msg;
            size_t err_sz      = regerror(error, &func_regex, NULL, 0);
            error_msg          = new char[err_sz];
            err_sz = regerror(error, &func_regex, error_msg, err_sz);
            std::string msg    = std::string("Failed to compile regular expression : ") + std::string(error_msg);
            delete[] error_msg;
            throw std::runtime_error(msg);
        }
        nmatch_number_regex    = number_regex.re_nsub;
        regmatch_t * sm_func   = (regmatch_t *)malloc(sizeof(regmatch_t) * nmatch_func_regex);
        regmatch_t * sm_var    = (regmatch_t *)malloc(sizeof(regmatch_t) * nmatch_var_regex);
        regmatch_t * sm_number = (regmatch_t *)malloc(sizeof(regmatch_t) * nmatch_number_regex);
#else
        std::regex  func_regex(R"RAW(\w+ *\()RAW");
        std::regex  var_regex(R"RAW(\{[_a-zA-Z1-90]+\})RAW");
        std::regex  number_regex(R"RAW([\d.]+(?:[eE][+-]?\d+)?)RAW");
        std::smatch sm_var, sm_func, sm_number;
#endif
        while (index < expr.size()) {
            switch (expr[index]) {
            case '(': {
                m_tokens.push_back({"(", LEFT_BRACKET});
                index += 1;
            } break;
            case ')': {
                m_tokens.push_back({")", RIGHT_BRACKET});
                index += 1;
            } break;
            case '+': {
                m_tokens.push_back({"+", OPERATOR});
                index += 1;
            } break;
            case '-': {
                m_tokens.push_back({"-", OPERATOR});
                index += 1;
            } break;
            case '*': {
                if ((index < expr.size() - 1) && (expr[index + 1] == '*')) { // Cas de l'opérateur puissance :
                    m_tokens.push_back({"**", OPERATOR});
                    index += 2;
                } else {
                    m_tokens.push_back({"*", OPERATOR});
                    index += 1;
                }
            } break;
            case ',': {
                m_tokens.push_back({",", SEPARATOR});
                index += 1;
            } break;
            case '/': {
                m_tokens.push_back({"/", OPERATOR});
                index += 1;
            } break;
            case '<': {
                if ((index < expr.size() - 1) && (expr[index + 1] == '=')) {
                    m_tokens.push_back({"<=", OPERATOR});
                    index += 2;
                } else {
                    m_tokens.push_back({"<", OPERATOR});
                    index += 1;
                }
            } break;
            case '>': {
                if ((index < expr.size() - 1) && (expr[index + 1] == '=')) {
                    m_tokens.push_back({">=", OPERATOR});
                    index += 2;
                } else {
                    m_tokens.push_back({">", OPERATOR});
                    index += 1;
                }
            } break;
            case '=': {
                if ((index < expr.size() - 1) && (expr[index + 1] == '=')) {
                    m_tokens.push_back({"==", OPERATOR});
                    index += 2;
                } else {
                    m_tokens.push_back({"=", OPERATOR});
                    index += 1;
                }
            } break;
            case '!': {
                if ((index >= expr.size() - 1) || (expr[index + 1] != '=')) {
                    throw std::logic_error("Wrong operator : ! doesn't exist. Want to use != ? ");
                }
                m_tokens.push_back({"!=", OPERATOR});
            }
            case '{': // Cas d'une variable
            {
                std::string substr(expr, index);
#if defined(USE_C_REGEX)
                int match;
                match = regexec(&var_regex, substr.c_str(), nmatch_var_regex, sm_var, 0);
                if (match != 0) {
                    if (match == REG_NOMATCH) {
                        std::ostringstream ostr;
                        ostr << expr << " : Ill-formed variable at position " << index << " !";
                        throw std::logic_error(ostr.str());
                    } else {
                        char * text;
                        size_t size;
                        size = regerror(match, &var_regex, NULL, 0);
                        text = new char[size];
                        if (text) {
                            regerror(match, &var_regex, text, size);
                            std::string msg(text);
                            delete[] text;
                            throw std::logic_error(msg);
                        } else {
                            throw std::bad_alloc();
                        }
                    }
                }
                int    start = sm_var[0].rm_so;
                int    end   = sm_var[0].rm_eo;
                size_t size  = size_t(end - start);
                m_tokens.push_back({std::string(expr, index + 1, size - 2), VARIABLE});
                index += size;
#else
                bool b = std::regex_search(substr, sm_var, var_regex);
                if (not b) {
                    std::ostringstream ostr;
                    ostr << expr << " : Ill-formed variable at position " << index << " !";
                    throw std::logic_error(ostr.str());
                }
                m_tokens.push_back({std::string(expr, index + 1, sm_var.length() - 2), VARIABLE});
                index += sm_var.length();
#endif
            } break;
            case ' ':
            case '\\':
            case '\n': {
                index += 1;
            } break;
            default: // Cas des fonctions et des nombres et de pi ! :
            {
                std::string expr_v = std::string(expr, index);
                std::transform(expr_v.begin(), expr_v.end(), expr_v.begin(), ::tolower);
                if (std::isalpha(expr[index])) { // On a affaire à une fonction ou à pi !
                    if (((expr_v[0] == 'p') and (expr_v[1] == 'i')) and
                        ((index > expr.size() - 2) or (not std::isalpha(expr[index + 2])))) {
                        m_tokens.push_back({"pi", PI_VALUE});
                        index += 2;
                    } else if (expr_v.find("minimum") == 0) // Min case
                    {
                        m_tokens.push_back({"min", MINMAXOP});
                        index += 7;
                    } else if (expr_v.find("maximum") == 0) // Max case
                    {
                        m_tokens.push_back({"max", MINMAXOP});
                        index += 7;
                    } else if (expr_v.find("logical_and") == 0) // logical_and case :
                    {
                        m_tokens.push_back({"logical_and", MINMAXOP});
                        index += 11;
                    } else if (expr_v.find("logical_or") == 0) // logical_and case :
                    {
                        m_tokens.push_back({"logical_or", MINMAXOP});
                        index += 10;
                    } else if (expr_v.find("logical_xor") == 0) // logical_and case :
                    {
                        m_tokens.push_back({"logical_xor", MINMAXOP});
                        index += 11;
                    } else {
                        std::string substr(expr, index);
#if defined(USE_C_REGEX)
                        int match;
                        match = regexec(&func_regex, substr.c_str(), nmatch_func_regex, sm_func, 0);
                        if (match != 0) {
                            if (match == REG_NOMATCH) {
                                std::ostringstream ostr;
                                ostr << expr << " : Ill-formed function at position " << index << " !";
                                throw std::logic_error(ostr.str());
                            } else {
                                char * text;
                                size_t size;
                                size = regerror(match, &func_regex, NULL, 0);
                                text = new char[size];
                                if (text) {
                                    regerror(match, &func_regex, text, size);
                                    std::string msg(text);
                                    delete[] text;
                                    throw std::logic_error(msg);
                                } else {
                                    throw std::bad_alloc();
                                }
                            }
                        }
                        int         start     = sm_func[0].rm_so;
                        int         end       = sm_func[0].rm_eo;
                        size_t      size      = size_t(end - start);
                        std::string func_name = std::string(expr, index, size - 1);
                        func_name.erase(std::remove_if(func_name.begin(), func_name.end(), ::isspace), func_name.end());
                        m_tokens.push_back({func_name, FUNCTION});
                        index += size - 1;
#else
                        bool b = std::regex_search(substr, sm_func, func_regex);
                        if (not b) {
                            std::ostringstream ostr;
                            ostr << expr << " : Ill-formed function at position " << index << " !";
                            throw std::logic_error(ostr.str());
                        }
                        std::string func_name = std::string(expr, index, sm_func.length() - 1);
                        func_name.erase(std::remove_if(func_name.begin(), func_name.end(), ::isspace), func_name.end());
                        m_tokens.push_back({func_name, FUNCTION});
                        index += sm_func.length() - 1;
#endif
                    }
                } else if (std::isdigit(expr[index])) { // On a affaire à un nombre !
                    std::string substr(expr, index);
#if defined(USE_C_REGEX)
                        int match;
                        match = regexec(&number_regex, substr.c_str(), nmatch_number_regex, sm_number, 0);
                        if (match != 0) {
                            if (match == REG_NOMATCH) {
                                std::ostringstream ostr;
                                ostr << expr << " : Ill-formed number at position " << index << " !";
                                throw std::logic_error(ostr.str());
                            } else {
                                char * text;
                                size_t size;
                                size = regerror(match, &number_regex, NULL, 0);
                                text = new char[size];
                                if (text) {
                                    regerror(match, &number_regex, text, size);
                                    std::string msg(text);
                                    delete[] text;
                                    throw std::logic_error(msg);
                                }
                            }
                        }
                        int         start     = sm_number[0].rm_so;
                        int         end       = sm_number[0].rm_eo;
                        size_t      size      = size_t(end - start);
                        m_tokens.push_back({std::string(expr, index, size), NUMBER});
                        index += size;
#else
                    bool b = std::regex_search(substr, sm_number, number_regex);
                    if (not b) {
                        std::ostringstream ostr;
                        ostr << expr << " : Ill-formed number at position " << index << " !";
                        throw std::logic_error(ostr.str());
                    }
                    m_tokens.push_back({std::string(expr, index, sm_number.length()), NUMBER});
                    index += sm_number.length();
#endif
                } else {
                    std::ostringstream ostr;
                    ostr << expr << " : Ill-formed expression at position " << index << " !";
                    throw std::logic_error(ostr.str());
                }
            }
            }
        }
#if defined(USE_C_REGEX)
        regfree(&var_regex);
        regfree(&number_regex);
        regfree(&func_regex);
        free(sm_number);
        free(sm_func);
        free(sm_var);
#endif
    }
} // namespace Expression
