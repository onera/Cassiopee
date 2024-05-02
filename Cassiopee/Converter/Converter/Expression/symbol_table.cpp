#include "Expression/symbol_table.hpp"

namespace {
	static std::shared_ptr<Expression::symbol_table> m_singleton;
}

namespace Expression {
symbol_table& symbol_table::get()
{
	if ( not m_singleton )
		m_singleton = std::make_shared<symbol_table>();
	return *m_singleton;
}
}
