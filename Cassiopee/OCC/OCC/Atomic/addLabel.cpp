/*    
    Copyright 2013-2024 Onera.

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
// add label to CAD
#include "occ.h"

#include <XCAFDoc_ShapeTool.hxx>
#include <TDF_Label.hxx>
#include <TopoDS_Shape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs.hxx>

// Assuming you have a document and a shape
Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(myDocument->Main());
TopoDS_Shape myShape = ...; // Your shape

// Iterate over faces
for (TopExp_Explorer exp(myShape, TopAbs_FACE); exp.More(); exp.Next()) 
{
  TopoDS_Shape face = exp.Current();
  TDF_Label label;
  if (shapeTool->FindShape(face, label)) 
  {
    // Do something with the label
  }
}
