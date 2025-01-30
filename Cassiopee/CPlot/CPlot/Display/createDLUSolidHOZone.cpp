/*    
    Copyright 2013-2025 Onera.

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
#include "../DataDL.h"
#include "../ZoneImplDL.h"

//=============================================================================
/*
  Display une zone en solid ou en material.
*/
//=============================================================================
void DataDL::createGPUUSolidHOZone(UnstructZone *zonep, E_Int zone, E_Int zonet)
{
    E_Int ind, ind_elt;
    ZoneImplDL *zImpl = static_cast<ZoneImplDL *>( zonep->ptr_impl );
    GLenum error = glGetError();
    if (error != GL_NO_ERROR)
    {
        std::cerr << __PRETTY_FUNCTION__ << ": get error 0 n°0x" << std::hex << error << std::dec << std::flush << std::endl;
    }
    zImpl->_DLsolid = glGenLists( 1 );
    E_Int* connect = zonep->connect[0];
    E_Int stride = zonep->nec[0];
    E_Int eltSize = zonep->eltSize[0];
    E_Int ne = zonep->nec[0];

    double* x = zonep->x;
    double* y = zonep->y;
    double* z = zonep->z;
    
    glNewList(zImpl->_DLsolid, GL_COMPILE);
    glPatchParameteri(GL_PATCH_VERTICES, GLint(eltSize));
    glBegin(GL_PATCHES);

    for (E_Int ielts = 0; ielts < ne; ++ielts) 
    {
        ind_elt = ielts;
        for ( E_Int inode = 0; inode < eltSize; inode++ ) 
        {
            ind = connect[ ind_elt + inode * stride ] - 1;
            glVertex3f( (float)x[ ind ], (float)y[ ind ], (float)z[ ind ] );
        }
    }
    glEnd();
    glEndList();
}
