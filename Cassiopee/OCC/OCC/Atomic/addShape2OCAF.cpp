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
#include "TDF_Tool.hxx"

//=====================================================================
// Add a shape (shape, compound) to OCAF
//=====================================================================
void K_OCC::addShape2OCAF(TopoDS_Shape& shape, char* labelName, TDocStd_Document& doc)
{
  Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc.Main());  

  // if label already exists, AddShape make a compound of previous and new
  TDF_Label label = shapeTool->AddShape(shape);
  TDataStd_Name::Set(label, labelName);
  return;

  /*
  // test if label exists
  E_Bool exists = false;

  TDF_LabelSequence labels;
  shapeTool->GetShapes(labels);

  TDF_Label label;
  for (Standard_Integer i = 1; i <= labels.Length(); i++)
  {
    label = labels.Value(i); // one label of collection

    TCollection_AsciiString es;
    TDF_Tool::Entry(label, es); // retourne l'entry de la label (3:0:3 indiquant sa position dans la hierarchie)

    Handle(TDataStd_Name) name = new TDataStd_Name(); 
    if (label.FindAttribute(TDataStd_Name::GetID(), name)) // retourne tous les attributs de type string
    { 
      TCollection_ExtendedString labelName2 = name->Get();
      if (TCollection_ExtendedString(labelName) == labelName2) 
      { exists = true; break; } 
    }
  }

  if (exists)
  {
    TopoDS_Shape prevShape = shapeTool->GetShape(label);
    
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);

    {
      TopTools_IndexedMapOfShape sf = TopTools_IndexedMapOfShape();
      TopExp::MapShapes(prevShape, TopAbs_FACE, sf);
      TopTools_IndexedMapOfShape se = TopTools_IndexedMapOfShape();
      TopExp::MapShapes(prevShape, TopAbs_EDGE, se);

      for (E_Int i = 1; i <= sf.Extent(); i++)
      {
        TopoDS_Face F = TopoDS::Face(sf(i));
        builder.Add(compound, F);
      }
      for (E_Int i = 1; i <= se.Extent(); i++)
      {
        TopoDS_Edge E = TopoDS::Edge(se(i));
        builder.Add(compound, E);
      }
    }

    {
      TopTools_IndexedMapOfShape sf = TopTools_IndexedMapOfShape();
      TopExp::MapShapes(shape, TopAbs_FACE, sf);
      TopTools_IndexedMapOfShape se = TopTools_IndexedMapOfShape();
      TopExp::MapShapes(shape, TopAbs_EDGE, se);

      for (E_Int i = 1; i <= sf.Extent(); i++)
      {
        TopoDS_Face F = TopoDS::Face(sf(i));
        builder.Add(compound, F);
      }
      for (E_Int i = 1; i <= se.Extent(); i++)
      {
        TopoDS_Edge E = TopoDS::Edge(se(i));
        builder.Add(compound, E);
      }
    }
    TDF_Label label = shapeTool->AddShape(compound);
    TDataStd_Name::Set(label, labelName);
  }
  else
  {
    // set shape to doc as new label
    TDF_Label label = shapeTool->AddShape(shape);
    TDataStd_Name::Set(label, labelName);
  }
  */
}
