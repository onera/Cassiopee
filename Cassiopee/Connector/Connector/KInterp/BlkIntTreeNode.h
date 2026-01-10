/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _CONNECTOR_BLKINTTREENODE_H_
#define _CONNECTOR_BLKINTTREENODE_H_
// ============================================================================
/*
   The node of the ADT (binary tree of cells).
*/
//=============================================================================
namespace K_KINTERP
{
class BlkIntTreeNode
{
  public:
    BlkIntTreeNode();
    BlkIntTreeNode(E_Int ind, E_Float xmax, E_Float ymax, E_Float zmax,
                   E_Float xmin, E_Float ymin, E_Float zmin);
    ~BlkIntTreeNode();
    
    void setLeftAt(BlkIntTreeNode*);
    void setRightAt(BlkIntTreeNode*);
    BlkIntTreeNode* getLeft() const;
    BlkIntTreeNode* getRight() const;
    BlkIntTreeNode** getLeftPt();
    BlkIntTreeNode** getRightPt();
    void getCellBB(E_Float& xmax, E_Float& ymax, E_Float& zmax,
                   E_Float& xmin, E_Float& ymin, E_Float& zmin) const;
    E_Int getIndex() const;

  public: // for speed
    E_Int _ind;
    E_Float _BB[6];
    K_KINTERP::BlkIntTreeNode* _left;
    K_KINTERP::BlkIntTreeNode* _right;
};
}
// ============================================================================
// INLINE
inline void K_KINTERP::BlkIntTreeNode::setLeftAt(K_KINTERP::BlkIntTreeNode* left)
{
   _left = left;
}

//=============================================================================
inline void K_KINTERP::BlkIntTreeNode::setRightAt(K_KINTERP::BlkIntTreeNode* right)
{
   _right = right;
}

//=============================================================================
inline K_KINTERP::BlkIntTreeNode* K_KINTERP::BlkIntTreeNode::getLeft() const
{
   return _left;
}

//=============================================================================
inline K_KINTERP::BlkIntTreeNode* K_KINTERP::BlkIntTreeNode::getRight() const
{
   return _right;
}

//=============================================================================
inline K_KINTERP::BlkIntTreeNode** K_KINTERP::BlkIntTreeNode::getLeftPt()
{
   return &_left;
}

//=============================================================================
inline K_KINTERP::BlkIntTreeNode** K_KINTERP::BlkIntTreeNode::getRightPt()
{
   return &_right;
}

//=============================================================================
inline void K_KINTERP::BlkIntTreeNode::getCellBB(E_Float& xmax, E_Float& ymax,
                                                   E_Float& zmax, E_Float& xmin,
                                                   E_Float& ymin, E_Float& zmin) const
{
  xmax = _BB[0];
  ymax = _BB[1];
  zmax = _BB[2];
  xmin = _BB[3];
  ymin = _BB[4];
  zmin = _BB[5];
}

//=============================================================================
inline E_Int K_KINTERP::BlkIntTreeNode::getIndex() const
{
  return _ind;
}
#endif
// ============================================================================


