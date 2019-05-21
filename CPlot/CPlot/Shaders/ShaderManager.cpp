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
#include "ShaderManager.h"
#include "../Data.h"
#include "FragmentShader.h"
#include "GeomShader.h"
#include "TesselationControlShader.hpp"
#include "TesselationEvaluationShader.hpp"
#include "VertexShader.h"
#include <stdexcept>
#include <string.h>
#include <iostream>
using namespace CPlot;

//=============================================================================
ShaderManager::ShaderManager()
    : _shaderList(),
      _currentActiveShader( 0 ),
      m_previous_shader( nullptr )
{
}
//=============================================================================
ShaderManager::~ShaderManager()
{
    deactivate();
    for ( std::vector<Shader *>::iterator itShad = _shaderList.begin();
          itShad != _shaderList.end(); itShad++ ) {
        delete ( *itShad );
    }
}

//=============================================================================
/* Cree un shader a partir de 1 fichier geomShader */
//==============================================================================
Shader *ShaderManager::addFromFile( const char *geomFile, const char *vertexFile, const char *fragmentFile )
{
    Shader *shad = new Shader;
    if ( vertexFile   != nullptr ) shad->add( std::make_shared<VertexShader>  () );
    if ( geomFile     != nullptr ) shad->add( std::make_shared<GeomShader>    () );
    if ( fragmentFile != nullptr ) shad->add( std::make_shared<FragmentShader>() );
    try
    {
        if ( (vertexFile != nullptr) && not shad->vertex_shader->compile( ShaderObject::load( vertexFile ) ) )
        {
            std::cerr << "Failed compiling vertex shader " << vertexFile << " : " << std::endl;
            std::cerr << shad->vertex_shader->getCompilerLog().c_str() << std::endl;
            delete shad;
            return nullptr;
        }
        if ( ( geomFile != nullptr ) && ( not shad->geometry_shader->compile( ShaderObject::load( geomFile ) ) ) )
        {
            std::cerr << "Failed compiling geometry shader " << geomFile << " : " << std::endl;
            std::cerr << shad->geometry_shader->getCompilerLog().c_str() << std::endl;
            delete shad;
            return nullptr;            
        }
        if ( ( fragmentFile != nullptr ) && ( not shad->fragment_shader->compile( ShaderObject::load( fragmentFile ) ) ) )
        {
            std::cerr << "Failed compiling fragment shader " << fragmentFile << " : " << std::endl;
            std::cerr << shad->fragment_shader->getCompilerLog().c_str() << std::endl;
            delete shad;
            return nullptr;
        }
    } catch ( std::runtime_error &err ) {
        std::cerr << "Failed to load shader : " << err.what() << std::endl;
        return nullptr;
    }
    // Rajout shader sans tesselation
    _shaderList.push_back( shad );
    // Puis rajout des shaders avec tesselation :
    for ( int idTess = 1; idTess <= tesselationManager.numberOfShaders(); ++idTess )
    {
        auto tes_shaders = tesselationManager[ idTess ];
        Shader *tes_pt_shader = new Shader;
        if ( shad->vertex_shader != nullptr )
            tes_pt_shader->add( shad->vertex_shader );
        if ( shad->geometry_shader != nullptr )
            tes_pt_shader->add( shad->geometry_shader );
        tes_pt_shader->add( tes_shaders.first );
        tes_pt_shader->add( tes_shaders.second );
        if ( shad->fragment_shader != nullptr )
            tes_pt_shader->add( shad->fragment_shader );
        _shaderList.push_back( tes_pt_shader );
    }

    return shad;
}
//=============================================================================
/* Cree un shader a partir de 2 fichiers vertexShader, fragmentShader */
//==============================================================================
Shader *ShaderManager::addFromFile( const char *vertexFile,
                                    const char *fragmentFile )
{
    Shader *shad = new Shader;
    if ( vertexFile   != nullptr ) shad->add( std::make_shared<VertexShader>  () );
    if ( fragmentFile != nullptr ) shad->add( std::make_shared<FragmentShader>() );
    try
    {
        if ( (vertexFile != nullptr) && not shad->vertex_shader->compile( ShaderObject::load( vertexFile ) ) )
        {
            std::cerr << "Failed compiling vertex shader " << vertexFile << " : " << std::endl;
            std::cerr << shad->vertex_shader->getCompilerLog().c_str() << std::endl;
            delete shad;
            return nullptr;
        }
        if ( ( fragmentFile != nullptr ) && ( not shad->fragment_shader->compile( ShaderObject::load( fragmentFile ) ) ) )
        {
            std::cerr << "Failed compiling fragment shader " << fragmentFile << " : " << std::endl;
            std::cerr << shad->fragment_shader->getCompilerLog().c_str() << std::endl;
            delete shad;
            return nullptr;
        }
    } catch ( std::runtime_error &err ) {
        std::cerr << "Failed to load shader : " << err.what() << std::endl;
        return nullptr;
    }
    _shaderList.push_back( shad );
    // Puis rajout des shaders avec tesselation :
    for ( int idTess = 1; idTess <= tesselationManager.numberOfShaders(); ++idTess )
    {
        auto tes_shaders = tesselationManager[ idTess ];
        Shader *tes_pt_shader = new Shader;
        if ( shad->vertex_shader != nullptr )
            tes_pt_shader->add( shad->vertex_shader );
        tes_pt_shader->add( tes_shaders.first );
        tes_pt_shader->add( tes_shaders.second );
        if ( shad->fragment_shader != nullptr )
            tes_pt_shader->add( shad->fragment_shader );
        _shaderList.push_back( tes_pt_shader );
    }
    return shad;
}
//=============================================================================
unsigned short ShaderManager::getId( Shader *shad ) const
{
    unsigned short id = 0;
    std::vector<Shader *>::const_iterator it = _shaderList.begin();
    while ( it != _shaderList.end() ) {
        if ( ( *it ) == shad ) break;
        id++;
    }
    if ( id == _shaderList.size() ) return 0;
    return id;
}
// ============================================================================
bool ShaderManager::eraseShader( Shader *obj )
{
    std::vector<Shader *>::iterator it = _shaderList.begin();
    while ( it != _shaderList.end() ) {
        if ( ( *it ) == obj ) {
            _shaderList.erase( it );
            delete obj;
            return true;
        }
    }
    return false;
}
//=============================================================================
void ShaderManager::set_tesselation( unsigned short idTess )
{
    if ( idTess == 0 ) {
        unset_tesselation();
        return;
    }
    // On regarde tout d'abord si on avait deja une tesselation active :
    if ( tesselationManager.currentShader() == idTess ) return;
    tesselationManager.activate( idTess );
}
// ----------------------------------------------------------------------------
void ShaderManager::unset_tesselation()
{
    tesselationManager.deactivate();
}
// ----------------------------------------------------------------------------
void ShaderManager::activate( unsigned short id )
{
        if ( m_previous_shader != nullptr ) {
            m_previous_shader->end();
            m_previous_shader = nullptr;
        }
        //if ( _currentActiveShader > 0 ) _shaderList[ _currentActiveShader - 1 ]->end();
        _currentActiveShader = 0;
        if ( id == 0 ) return;
        if ( id >= _shaderList.size() ) return;
        _currentActiveShader = id;
        if ( !_shaderList[ id  ]->start() )
            throw std::runtime_error( "Fail to start shader !" );
        else
            m_previous_shader = _shaderList[ id ];
}
//=============================================================================
void ShaderManager::activate( Shader *shad )
{
    if ( m_previous_shader != nullptr ) m_previous_shader->end();
    //if ( _currentActiveShader > 0 ) _shaderList[ _currentActiveShader - 1 ]->end();
    _currentActiveShader = 0;
    unsigned short id = 0;
    std::vector<Shader *>::iterator it = _shaderList.begin();
    while ( it != _shaderList.end() ) {
        id++;
        if ( ( *it ) == shad ) {
            activate( id );
            _currentActiveShader = id;
            m_previous_shader = shad;
            return;
        }
    }
}
//=============================================================================
void ShaderManager::deactivate()
{
    if ( m_previous_shader != nullptr ) m_previous_shader->end();
    activate( (unsigned short)0 );
    _currentActiveShader = 0;
}

//=============================================================================
int ShaderManager::init()
{
    glewInit();
    if ( GLEW_ARB_vertex_shader && GLEW_ARB_fragment_shader )
        return 1;  // success
    else
        return 0;  // fail
}

//=============================================================================
int ShaderManager::load()
{
    tesselationManager.load();

    char vert[ 256 * 8 ];
    char frag[ 256 * 8 ];
    char geom[ 256 * 8 ];
    Data *d = Data::getInstance();
    char *path = d->ptrState->shaderPath;

    // - 0 - Pas de shader materiel, mais shader tesselation
    strcpy( vert, path );
    strcat( vert, "simple.vert" );
    addFromFile( vert, nullptr );

    // - 1 - Phong unidirectionnel
    strcpy( vert, path );
    strcat( vert, "phong.vert" );
    strcpy( frag, path );
    strcat( frag, "phong.frag" );
    addFromFile( vert, frag );

    // - 2 - Phong bidirectionnel
    strcpy( vert, path );
    strcat( vert, "phong2s.vert" );
    strcpy( frag, path );
    strcat( frag, "phong2s.frag" );
    addFromFile( vert, frag );

    // - 3 - Glass (environ mapping - utilise la texture environ + frame buffer)
    strcpy( vert, path );
    strcat( vert, "glass.vert" );
    strcpy( frag, path );
    strcat( frag, "glass.frag" );
    addFromFile( vert, frag );

    // - 4 - Chrome (environ mapping - utilise la texture environ)
    strcpy( vert, path );
    strcat( vert, "envmap.vert" );
    strcpy( frag, path );
    strcat( frag, "envmap.frag" );
    addFromFile( vert, frag );

    // - 5 - Metal (anisotropic)
    strcpy( vert, path );
    strcat( vert, "anisotropic.vert" );
    strcpy( frag, path );
    strcat( frag, "anisotropic.frag" );
    addFromFile( vert, frag );

    // - 6 - Wood (utilise la texture noise3D)
    strcpy( vert, path );
    strcat( vert, "wood.vert" );
    strcpy( frag, path );
    strcat( frag, "wood.frag" );
    addFromFile( vert, frag );

    // - 7 - Marble (utilise la texture noise3D)
    strcpy( vert, path );
    strcat( vert, "marble.vert" );
    strcpy( frag, path );
    strcat( frag, "marble.frag" );
    addFromFile( vert, frag );

    // - 8 - Smoke (pas OK)
    strcpy( vert, path );
    strcat( vert, "raymarching.vert" );
    strcpy( frag, path );
    strcat( frag, "raymarching.frag" );
    addFromFile( vert, frag );

    // - 9 - XRay
    strcpy( vert, path );
    strcat( vert, "xray.vert" );
    strcpy( frag, path );
    strcat( frag, "xray.frag" );
    addFromFile( vert, frag );

    // - 10 - Iso banded-colormap shader
    strcpy( vert, path );
    strcat( vert, "iso.vert" );
    strcpy( frag, path );
    strcat( frag, "iso.frag" );
    addFromFile( vert, frag );

    // - 11 - Granite (utilise la texture noise3D)
    strcpy( vert, path );
    strcat( vert, "granite.vert" );
    strcpy( frag, path );
    strcat( frag, "granite.frag" );
    addFromFile( vert, frag );

    // - 12 - Sphere billboarding
    strcpy( vert, path );
    strcat( vert, "spheres.vert" );
    strcpy( frag, path );
    strcat( frag, "spheres.frag" );
    addFromFile( vert, frag );

    // - 13 - Anaglyph shader (monochrome)
    strcpy( vert, path );
    strcat( vert, "anaglyph.vert" );
    strcpy( frag, path );
    strcat( frag, "anaglyph.frag" );
    addFromFile( vert, frag );

    // - 14 - Anaglyph shader (color)
    strcpy( vert, path );
    strcat( vert, "anaglyphColor.vert" );
    strcpy( frag, path );
    strcat( frag, "anaglyphColor.frag" );
    addFromFile( vert, frag );

    // - 15 - Iso continuous-colormap shader
    strcpy( vert, path );
    strcat( vert, "iso2.vert" );
    strcpy( frag, path );
    strcat( frag, "iso2.frag" );
    addFromFile( vert, frag );

    // - 16 - Brick shader
    strcpy( vert, path );
    strcat( vert, "brick.vert" );
    strcpy( frag, path );
    strcat( frag, "brick.frag" );
    addFromFile( vert, frag );

    // - 17 - Cloud shader
    strcpy( vert, path );
    strcat( vert, "cloud.vert" );
    strcpy( frag, path );
    strcat( frag, "cloud.frag" );
    addFromFile( vert, frag );

    // - 18 - iso + granite
    strcpy( vert, path );
    strcat( vert, "isoGranite.vert" );
    strcpy( frag, path );
    strcat( frag, "isoGranite.frag" );
    addFromFile( vert, frag );

    // - 19 - shadow mapping
    strcpy( vert, path );
    strcat( vert, "phong2s.vert" );
    strcpy( frag, path );
    strcat( frag, "phong2s.frag" );
    addFromFile( vert, frag );

    // - 20 - DOF!
    strcpy( vert, path );
    strcat( vert, "dof.vert" );
    strcpy( frag, path );
    strcat( frag, "dof.frag" );
    addFromFile( vert, frag );

    // - 21 - Silhouette + gooch shader
    strcpy( vert, path );
    strcat( vert, "gooch.vert" );
    strcpy( frag, path );
    strcat( frag, "gooch.frag" );
    addFromFile( vert, frag );

    // - 22 - Flat shader
    strcpy( vert, path );
    strcat( vert, "flat.vert" );
    strcpy( frag, path );
    strcat( frag, "flat.frag" );
    addFromFile( vert, frag );

    // - 23 - Billboard shader
    strcpy( vert, path );
    strcat( vert, "billboard.vert" );
    strcpy( frag, path );
    strcat( frag, "billboard.frag" );
    addFromFile( vert, frag );

    // - 24 - iso + flat
    strcpy( vert, path );
    strcat( vert, "isoFlat.vert" );
    strcpy( frag, path );
    strcat( frag, "isoFlat.frag" );
    addFromFile( vert, frag );

    // - 25 - iso + chrome
    strcpy( vert, path );
    strcat( vert, "isoEnvmap.vert" );
    strcpy( frag, path );
    strcat( frag, "isoEnvmap.frag" );
    addFromFile( vert, frag );

    // - 26 - iso + glass
    strcpy( vert, path );
    strcat( vert, "isoGlass.vert" );
    strcpy( frag, path );
    strcat( frag, "isoGlass.frag" );
    addFromFile( vert, frag );

    // - 27 - vector rgb direct
    strcpy( vert, path );
    strcat( vert, "rgb.vert" );
    strcpy( frag, path );
    strcat( frag, "rgb.frag" );
    addFromFile( vert, frag );

    // - 28 - iso + brick
    strcpy( vert, path );
    strcat( vert, "isoBrick.vert" );
    strcpy( frag, path );
    strcat( frag, "isoBrick.frag" );
    addFromFile( vert, frag );

    // - 29 - iso colored lines
    strcpy( vert, path );
    strcat( vert, "iso3.vert" );
    strcpy( frag, path );
    strcat( frag, "iso3.frag" );
    addFromFile( vert, frag );

    // - 30 - iso + xray
    strcpy( vert, path );
    strcat( vert, "isoXRay.vert" );
    strcpy( frag, path );
    strcat( frag, "isoXRay.frag" );
    addFromFile( vert, frag );

    // - 31 - iso + gooch
    strcpy( vert, path );
    strcat( vert, "isoGooch.vert" );
    strcpy( frag, path );
    strcat( frag, "isoGooch.frag" );
    addFromFile( vert, frag );

    // - 32 - iso + metal (anistropic)
    strcpy( vert, path );
    strcat( vert, "isoAnisotropic.vert" );
    strcpy( frag, path );
    strcat( frag, "isoAnisotropic.frag" );
    addFromFile( vert, frag );

    // - 33 - Vector line shader ( geom + frag + vert )
    strcpy( geom, path );
    strcat( geom, "streamline.geom" );
    strcpy( vert, path );
    strcat( vert, "streamline.vert" );
    strcpy( frag, path );
    strcat( frag, "streamline.frag" );
    addFromFile( geom, vert, frag );

    // - 34 - Vector triangle shader ( geom + frag + vert (34) )
    strcpy( geom, path );
    strcat( geom, "streamarrow.geom" );
    strcpy( vert, path );
    strcat( vert, "streamarrow.vert" );
    strcpy( frag, path );
    strcat( frag, "streamarrow.frag" );
    addFromFile( geom, vert, frag );

    // - 35 - Textured material shader
    strcpy( vert, path );
    strcat( vert, "texmat.vert" );
    strcpy( frag, path );
    strcat( frag, "texmat.frag" );
    addFromFile( vert, frag );

    return 1;
}
