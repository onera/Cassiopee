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

#include "TDF_Label.hxx"
#include "XCAFDoc_DocumentTool.hxx"
#include "TDocStd_Document.hxx"
#include "XCAFDoc_ShapeTool.hxx"
#include "XCAFDoc_ShapeMapTool.hxx"

//=====================================================================
// copy OCAF to top shape
//=====================================================================
TopoDS_Shape* K_OCC::copyOCAF2TopShape(TDocStd_Document& doc)
{
  // Get free shapes
  Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc.Main());
  TDF_LabelSequence labels;
  shapeTool->GetFreeShapes(labels);
  
  // Rebuild a single compound from topShape + shape
  BRep_Builder builder;
  TopoDS_Compound* compound = new TopoDS_Compound();
  builder.MakeCompound(*compound);
  
  for (Standard_Integer i = 1; i <= labels.Length(); i++)
  {
    TDF_Label label = labels.Value(i);
    TopoDS_Shape shape = shapeTool->GetShape(label);

    TopTools_IndexedMapOfShape sf = TopTools_IndexedMapOfShape();
    TopExp::MapShapes(shape, TopAbs_FACE, sf);
    TopTools_IndexedMapOfShape se = TopTools_IndexedMapOfShape();
    TopExp::MapShapes(shape, TopAbs_EDGE, se);

    for (E_Int i = 1; i <= sf.Extent(); i++)
    {
      TopoDS_Face F = TopoDS::Face(sf(i));
      builder.Add(*compound, F);
    }
    for (E_Int i = 1; i <= se.Extent(); i++)
    {
      TopoDS_Edge E = TopoDS::Edge(se(i));
      builder.Add(*compound, E);
    }
  }
  return compound;
}
