/*    
    Copyright 2013-2019 Onera.

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
#ifndef _KCORE_OCTREE_H_
#define _KCORE_OCTREE_H_
#include "Def/DefTypes.h"
namespace K_SEARCH
{
//=============================================================================
/* Classe definissant un noeud d octree */
//=============================================================================
class OctreeNode
{
  public:

    ///+ 1- Constructor/Destructor
    OctreeNode(E_Float xmin, E_Float ymin, E_Float zmin, E_Float dh, E_Int level, 
               E_Int ind1=-1, E_Int ind2=-1, E_Int ind3=-1, E_Int ind4=-1,
               E_Int ind5=-1, E_Int ind6=-1, E_Int ind7=-1, E_Int ind8=-1);
    
    ~OctreeNode();

    E_Float getXmin();
    E_Float getYmin();
    E_Float getZmin();
    E_Float getDh();
    E_Int getInd1();
    E_Int getInd2();
    E_Int getInd3();
    E_Int getInd4();
    E_Int getInd5();
    E_Int getInd6();
    E_Int getInd7();
    E_Int getInd8();
    void setNext1(OctreeNode* fils);
    void setNext2(OctreeNode* fils);
    void setNext3(OctreeNode* fils);
    void setNext4(OctreeNode* fils);
    void setNext5(OctreeNode* fils);
    void setNext6(OctreeNode* fils);
    void setNext7(OctreeNode* fils);
    void setNext8(OctreeNode* fils);
    void setNext9(OctreeNode* fils);
    void setNext10(OctreeNode* fils);
    void setNext11(OctreeNode* fils);
    void setNext12(OctreeNode* fils);
    void setNext13(OctreeNode* fils);
    void setNext14(OctreeNode* fils);
    void setNext15(OctreeNode* fils);
    void setNext16(OctreeNode* fils);
    void setNext17(OctreeNode* fils);
    void setNext18(OctreeNode* fils);
    void setNext19(OctreeNode* fils);
    void setNext20(OctreeNode* fils);
    void setNext21(OctreeNode* fils);
    void setNext22(OctreeNode* fils);
    void setNext23(OctreeNode* fils);
    void setNext24(OctreeNode* fils);
    void setNext25(OctreeNode* fils);
    void setNext26(OctreeNode* fils);
    void setNext27(OctreeNode* fils);

    void setVoisin1(OctreeNode* voisin);
    void setVoisin2(OctreeNode* voisin);
    void setVoisin3(OctreeNode* voisin);
    void setVoisin4(OctreeNode* voisin);
    void setVoisin5(OctreeNode* voisin);
    void setVoisin6(OctreeNode* voisin);

    OctreeNode* getNext1();
    OctreeNode* getNext2();
    OctreeNode* getNext3();
    OctreeNode* getNext4();
    OctreeNode* getNext5();
    OctreeNode* getNext6();
    OctreeNode* getNext7();
    OctreeNode* getNext8();
    OctreeNode* getNext9();
    OctreeNode* getNext10();
    OctreeNode* getNext11();
    OctreeNode* getNext12();
    OctreeNode* getNext13();
    OctreeNode* getNext14();
    OctreeNode* getNext15();
    OctreeNode* getNext16();
    OctreeNode* getNext17();
    OctreeNode* getNext18();
    OctreeNode* getNext19();
    OctreeNode* getNext20();
    OctreeNode* getNext21();
    OctreeNode* getNext22();
    OctreeNode* getNext23();
    OctreeNode* getNext24();
    OctreeNode* getNext25();
    OctreeNode* getNext26();
    OctreeNode* getNext27();

    OctreeNode* getVoisin1();
    OctreeNode* getVoisin2();
    OctreeNode* getVoisin3();
    OctreeNode* getVoisin4();
    OctreeNode* getVoisin5();
    OctreeNode* getVoisin6();

    E_Int getLevel();

 private:
    E_Float _xmin; // bounding box de la cellule
    E_Float _ymin;
    E_Float _zmin;
    E_Float _dh;//spacing
    E_Int _level; // niveau : 0 le plus grossier, croissant vers les plus fins
    E_Int _ind1;
    E_Int _ind2;
    E_Int _ind3;
    E_Int _ind4;
    E_Int _ind5;
    E_Int _ind6;
    E_Int _ind7;
    E_Int _ind8;
    OctreeNode* _next1; // noeuds fils
    OctreeNode* _next2;
    OctreeNode* _next3;
    OctreeNode* _next4;
    OctreeNode* _next5;
    OctreeNode* _next6;
    OctreeNode* _next7;
    OctreeNode* _next8;
    OctreeNode* _next9; 
    OctreeNode* _next10;
    OctreeNode* _next11;
    OctreeNode* _next12;
    OctreeNode* _next13;
    OctreeNode* _next14;
    OctreeNode* _next15;
    OctreeNode* _next16;
    OctreeNode* _next17;
    OctreeNode* _next18;
    OctreeNode* _next19;
    OctreeNode* _next20;
    OctreeNode* _next21;
    OctreeNode* _next22;
    OctreeNode* _next23;
    OctreeNode* _next24;
    OctreeNode* _next25;
    OctreeNode* _next26;
    OctreeNode* _next27;

    OctreeNode* _voisin1; // noeuds voisin
    OctreeNode* _voisin2;
    OctreeNode* _voisin3;
    OctreeNode* _voisin4;
    OctreeNode* _voisin5;
    OctreeNode* _voisin6;

}; //end class 

  inline E_Int OctreeNode::getLevel(){return _level;}
  inline E_Float OctreeNode::getXmin(){return _xmin;}
  inline E_Float OctreeNode::getYmin(){return _ymin;}
  inline E_Float OctreeNode::getZmin(){return _zmin;}
  inline E_Float OctreeNode::getDh(){return _dh;}
  inline E_Int OctreeNode::getInd1(){return _ind1;}
  inline E_Int OctreeNode::getInd2(){return _ind2;}
  inline E_Int OctreeNode::getInd3(){return _ind3;}
  inline E_Int OctreeNode::getInd4(){return _ind4;}
  inline E_Int OctreeNode::getInd5(){return _ind5;}
  inline E_Int OctreeNode::getInd6(){return _ind6;}
  inline E_Int OctreeNode::getInd7(){return _ind7;}
  inline E_Int OctreeNode::getInd8(){return _ind8;}
  inline void OctreeNode::setNext1(OctreeNode* fils){_next1 = fils;}
  inline void OctreeNode::setNext2(OctreeNode* fils){_next2 = fils;}
  inline void OctreeNode::setNext3(OctreeNode* fils){_next3 = fils;}
  inline void OctreeNode::setNext4(OctreeNode* fils){_next4 = fils;}
  inline void OctreeNode::setNext5(OctreeNode* fils){_next5 = fils;}
  inline void OctreeNode::setNext6(OctreeNode* fils){_next6 = fils;}
  inline void OctreeNode::setNext7(OctreeNode* fils){_next7 = fils;}
  inline void OctreeNode::setNext8(OctreeNode* fils){_next8 = fils;}
  inline void OctreeNode::setNext9(OctreeNode* fils){_next9 = fils;}
  inline void OctreeNode::setNext10(OctreeNode* fils){_next10 = fils;}
  inline void OctreeNode::setNext11(OctreeNode* fils){_next11 = fils;}
  inline void OctreeNode::setNext12(OctreeNode* fils){_next12 = fils;}
  inline void OctreeNode::setNext13(OctreeNode* fils){_next13 = fils;}
  inline void OctreeNode::setNext14(OctreeNode* fils){_next14 = fils;}
  inline void OctreeNode::setNext15(OctreeNode* fils){_next15 = fils;}
  inline void OctreeNode::setNext16(OctreeNode* fils){_next16 = fils;}
  inline void OctreeNode::setNext17(OctreeNode* fils){_next17 = fils;}
  inline void OctreeNode::setNext18(OctreeNode* fils){_next18 = fils;}
  inline void OctreeNode::setNext19(OctreeNode* fils){_next19 = fils;}
  inline void OctreeNode::setNext20(OctreeNode* fils){_next20 = fils;}
  inline void OctreeNode::setNext21(OctreeNode* fils){_next21 = fils;}
  inline void OctreeNode::setNext22(OctreeNode* fils){_next22 = fils;}
  inline void OctreeNode::setNext23(OctreeNode* fils){_next23 = fils;}
  inline void OctreeNode::setNext24(OctreeNode* fils){_next24 = fils;}
  inline void OctreeNode::setNext25(OctreeNode* fils){_next25 = fils;}
  inline void OctreeNode::setNext26(OctreeNode* fils){_next26 = fils;}
  inline void OctreeNode::setNext27(OctreeNode* fils){_next27 = fils;}

  inline void OctreeNode::setVoisin1(OctreeNode* voisin){_voisin1 = voisin;}
  inline void OctreeNode::setVoisin2(OctreeNode* voisin){_voisin2 = voisin;}
  inline void OctreeNode::setVoisin3(OctreeNode* voisin){_voisin3 = voisin;}
  inline void OctreeNode::setVoisin4(OctreeNode* voisin){_voisin4 = voisin;}
  inline void OctreeNode::setVoisin5(OctreeNode* voisin){_voisin5 = voisin;}
  inline void OctreeNode::setVoisin6(OctreeNode* voisin){_voisin6 = voisin;}

  inline OctreeNode* OctreeNode::getNext1(){return _next1;}
  inline OctreeNode* OctreeNode::getNext2(){return _next2;}
  inline OctreeNode* OctreeNode::getNext3(){return _next3;}
  inline OctreeNode* OctreeNode::getNext4(){return _next4;}
  inline OctreeNode* OctreeNode::getNext5(){return _next5;}
  inline OctreeNode* OctreeNode::getNext6(){return _next6;}
  inline OctreeNode* OctreeNode::getNext7(){return _next7;}
  inline OctreeNode* OctreeNode::getNext8(){return _next8;}
  inline OctreeNode* OctreeNode::getNext9(){return _next9;}
  inline OctreeNode* OctreeNode::getNext10(){return _next10;}
  inline OctreeNode* OctreeNode::getNext11(){return _next11;}
  inline OctreeNode* OctreeNode::getNext12(){return _next12;}
  inline OctreeNode* OctreeNode::getNext13(){return _next13;}
  inline OctreeNode* OctreeNode::getNext14(){return _next14;}
  inline OctreeNode* OctreeNode::getNext15(){return _next15;}
  inline OctreeNode* OctreeNode::getNext16(){return _next16;}
  inline OctreeNode* OctreeNode::getNext17(){return _next17;}
  inline OctreeNode* OctreeNode::getNext18(){return _next18;}
  inline OctreeNode* OctreeNode::getNext19(){return _next19;}
  inline OctreeNode* OctreeNode::getNext20(){return _next20;}
  inline OctreeNode* OctreeNode::getNext21(){return _next21;}
  inline OctreeNode* OctreeNode::getNext22(){return _next22;}
  inline OctreeNode* OctreeNode::getNext23(){return _next23;}
  inline OctreeNode* OctreeNode::getNext24(){return _next24;}
  inline OctreeNode* OctreeNode::getNext25(){return _next25;}
  inline OctreeNode* OctreeNode::getNext26(){return _next26;}
  inline OctreeNode* OctreeNode::getNext27(){return _next27;}

  inline OctreeNode* OctreeNode::getVoisin1(){return _voisin1;}
  inline OctreeNode* OctreeNode::getVoisin2(){return _voisin2;}
  inline OctreeNode* OctreeNode::getVoisin3(){return _voisin3;}
  inline OctreeNode* OctreeNode::getVoisin4(){return _voisin4;}
  inline OctreeNode* OctreeNode::getVoisin5(){return _voisin5;}
  inline OctreeNode* OctreeNode::getVoisin6(){return _voisin6;}

} // end namespace
#endif
//========================= KCore/Search/OctreeNode.h =========================
