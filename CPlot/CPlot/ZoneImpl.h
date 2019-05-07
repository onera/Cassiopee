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
#ifndef _CPLOT_ZONEIMPL_H_
#define _CPLOT_ZONEIMPL_H_
#include "CPlotState.h"
#include "Shaders/shaders_id.h"

struct ZoneImpl {
    ZoneImpl( ) : _GPUResUse( 1 ) {}
    virtual ~ZoneImpl( ) {}
    int          _GPUResUse;  // =1 si cette zone utilise les ressources internes du GPU ( Display Lists ou VBO )
    virtual void freeGPURes( CPlotState* state, bool freeIso = true ) = 0;
    virtual void destroyIsoField( ) = 0;
};

#endif
