// test vector_view
#include <iostream>
#include <vector>
#include "vector_view.hpp"

int main( ) {
    std::vector<double> arr( 100 );
    for ( std::size_t        i = 0; i < arr.size( ); ++i ) arr[i] = 0.1 * ( i + 1 );
    CMP::vector_view<double> view_arr( arr );
    for ( CMP::vector_view<double>::iterator it = view_arr.begin( ); it != view_arr.end( ); ++it )
        std::cout << *it << " ";
    std::cout << std::endl;
    
    int indices[] = { 0,1,2, 3,1,2, 5,3,1, 3,4,1, 2,5,1 };
    CMP::vector_view<int> view_ind( indices, 15 );
    for ( std::size_t i = 0; i < view_ind.size(); ++ i ) std::cout << view_ind[i] << " ";
    std::cout << std::endl;
    
    return 0;
}
