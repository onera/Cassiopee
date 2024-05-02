#include <cassert>
#include <iostream>
#include "shared_ptr.hpp"

class CpteurRef {
public:
    static int nbRef;
    CpteurRef( ) { nbRef += 1; }
    CpteurRef( const CpteurRef& ref ) { nbRef += 1; }
    ~CpteurRef( ) { nbRef -= 1; }
};
int           CpteurRef::nbRef = 0;
std::ostream& operator<<( std::ostream& out, const CpteurRef& cpt ) {
    out << CpteurRef::nbRef;
    return out;
}

int main( ) {
    CMP::shared_ptr<int> pt_i( new int( 3 ) );
    CMP::shared_ptr<int> pt_j = pt_i;
    assert( pt_i == pt_j );
    assert( *pt_i == 3 );
    *pt_i = 2;
    assert( *pt_j == 2 );
    CMP::shared_ptr<int> pt_k( new int( 4 ) );
    pt_j.swap( pt_k );
    assert( *pt_j == 4 );
    assert( *pt_k == 2 );

    CMP::shared_ptr<CpteurRef> pt_ref  = new CpteurRef;
    CMP::shared_ptr<CpteurRef> pt_ref2 = pt_ref;
    assert( CpteurRef::nbRef == 1 );
    assert( pt_ref.use_count( ) == 2 );
    CMP::shared_ptr<CpteurRef> pt_ref3 = new CpteurRef;
    pt_ref2.swap( pt_ref3 );
    assert( CpteurRef::nbRef == 2 );
    assert( pt_ref.use_count( ) == 2 );
    assert( pt_ref2.use_count( ) == 1 );
    return EXIT_SUCCESS;
}
