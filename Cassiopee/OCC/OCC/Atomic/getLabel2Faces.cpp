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

#include "occ.h"
#include "TopoDS.hxx"
#include "TopoDS_Shape.hxx"
#include "TopoDS_Edge.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "BRep_Builder.hxx"

#include "XCAFDoc_ShapeTool.hxx"
#include "XCAFDoc_DocumentTool.hxx"
#include "XCAFDoc_ShapeMapTool.hxx"
#include "TDocStd_Document.hxx"
#include "TDataStd_Name.hxx"

//=====================================================================
// return label2Faces map
//=====================================================================
void K_OCC::getLabel2Faces(TDocStd_Document& doc, std::map< E_Int, std::vector<E_Int> >& label2Faces)
{
  Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc.Main());  
  TDF_LabelSequence labels;
  shapeTool->GetShapes(labels);

  E_Int istart = 1; E_Int iend = 1;
  for (Standard_Integer i = 1; i <= labels.Length(); i++)
  {
    TDF_Label label = labels.Value(i);
    TopoDS_Shape shape = shapeTool->GetShape(label);
    
    TopTools_IndexedMapOfShape faces = TopTools_IndexedMapOfShape();
    TopExp::MapShapes(shape, TopAbs_FACE, faces);

    iend = istart + faces.Extent()-1;

    std::vector<E_Int>& p = label2Faces[i];
    p.resize(faces.Extent());
    for (E_Int j = istart; j <= iend; j++) p[j-istart] = j-istart+1;
  }
  /*
  for (Standard_Integer i = 1; i <= labels.Length(); i++)
  {
    printf("compound %d\n", i);
    std::vector<E_Int>& p = label2Faces[i];
    for (size_t j = 0; j < p.size(); j++) printf("%d ", p[j]);
    printf("\n");
  }*/
}
