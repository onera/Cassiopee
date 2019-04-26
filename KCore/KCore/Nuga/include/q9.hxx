/*
 
 
 
              NUGA 
 
 
 
 */

#ifndef NUGA_Q9_HXX
#define NUGA_Q9_HXX

namespace NUGA
{
  class Q9
  {
  public:
    template <typename crd_t>
    static void splitQ4(crd_t& crd, const E_Int* nodes, E_Int* q41, E_Int* q42, E_Int* q43, E_Int* q44)
    {
      q41[0] = nodes[0]; q41[1] = nodes[4]; q41[2] = nodes[8]; q41[3] = nodes[7];
      q42[0] = nodes[4]; q42[1] = nodes[1]; q42[2] = nodes[5]; q42[3] = nodes[8];
      q43[0] = nodes[8]; q43[1] = nodes[5]; q43[2] = nodes[2]; q43[3] = nodes[6];
      q44[0] = nodes[7]; q44[1] = nodes[8]; q44[2] = nodes[6]; q44[3] = nodes[3];
    }
    template <typename crd_t>
    static void splitQ4T(crd_t& crd, const E_Int* nodes, E_Int* q41, E_Int* q42, E_Int* q43, E_Int* q44) // définition PGs fils pour tétra (convention)
    {
      q41[0] = nodes[0]; q41[1] = nodes[3]; q41[2] = nodes[5]; 
      q42[0] = nodes[3]; q42[1] = nodes[1]; q42[2] = nodes[4]; 
      q43[0] = nodes[5]; q43[1] = nodes[4]; q43[2] = nodes[2]; 
      q44[0] = nodes[3]; q44[1] = nodes[4]; q44[2] = nodes[5]; 
    }
  };
}

#endif
