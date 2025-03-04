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
#include "Shader.h"
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace CPlot;

//=============================================================================
Shader::Shader()
    : vertex_shader(), 
      tesselation_control_shader(), tesselation_evaluation_shader(), 
      geometry_shader(), fragment_shader(), 
      _isActive( false ), _isLinked( false ), _isModified( false ), _programId( 0 )
{
    _programId = glCreateProgramObjectARB();
    CHECK_GL_ERROR();
}
//==============================================================================
void Shader::unlink()
{
    if ( _isLinked ) {
        if ( vertex_shader != nullptr ) {
            glDetachShader( _programId, vertex_shader->getShaderId() );
            CHECK_GL_ERROR();
        }
        if ( tesselation_control_shader != nullptr ) {
            glDetachShader( _programId, tesselation_control_shader->getShaderId() );
            CHECK_GL_ERROR();
        }
        if ( tesselation_control_shader != nullptr ) {
            glDetachShader( _programId, tesselation_control_shader->getShaderId() );
            CHECK_GL_ERROR();
        }
        if ( geometry_shader != nullptr ) {
            glDetachShader( _programId, geometry_shader->getShaderId() );
            CHECK_GL_ERROR();
        }
        if ( fragment_shader != nullptr ) {
            glDetachShader( _programId, fragment_shader->getShaderId() );
            CHECK_GL_ERROR();
        }
    }
    _isLinked = false;
    _isModified = false;
}
//==============================================================================
Shader::~Shader()
{
    unlink();
    glDeleteObjectARB( _programId );
    CHECK_GL_ERROR();
}
// ============================================================================
void Shader::remove( Shader::stage st )
{
    _isModified = true;
    switch ( st ) {
        case VERTEX:
            vertex_shader.reset();
            break;
        case TESSELATION_CONTROL:
            tesselation_control_shader.reset();
            break;
        case TESSELATION_EVALUATION:
            tesselation_evaluation_shader.reset();
            break;
        case GEOMETRY:
            geometry_shader.reset();
            break;
        case FRAGMENT:
            fragment_shader.reset();
            break;
        default:
            std::cerr << "Warning : Wrong kind of shader stage..." << std::endl;
            _isModified = false;
    }
}
//==============================================================================
bool Shader::contain( const ShaderObject &proglet ) const
{
    return ( vertex_shader->getShaderId() == proglet.getShaderId() ) or
        ( tesselation_control_shader->getShaderId() == proglet.getShaderId() ) or
        ( tesselation_evaluation_shader->getShaderId() == proglet.getShaderId() ) or
        ( geometry_shader->getShaderId() == proglet.getShaderId() ) or
        ( fragment_shader->getShaderId() == proglet.getShaderId() );
}
//==============================================================================
int Shader::nb_proglets() const
{
    return ( vertex_shader != nullptr ? 1 : 0 ) +
        ( tesselation_control_shader != nullptr ? 1 : 0 ) +
        ( tesselation_evaluation_shader != nullptr ? 1 : 0 ) +
        ( geometry_shader != nullptr ? 1 : 0 ) +
        ( fragment_shader != nullptr ? 1 : 0 );
}
//==============================================================================
bool Shader::link()
{
    if ( _isLinked ) return true;
    if ( vertex_shader != nullptr ) {
        glAttachShader( _programId, vertex_shader->getShaderId() );
        CHECK_GL_ERROR();
    }
    if ( tesselation_control_shader != nullptr ) {
        glAttachShader( _programId, tesselation_control_shader->getShaderId() );
        CHECK_GL_ERROR();
    }
    if ( tesselation_evaluation_shader != nullptr ) {
        glAttachShader( _programId, tesselation_evaluation_shader->getShaderId() );
        CHECK_GL_ERROR();
    }
    if ( geometry_shader != nullptr ) {
        glAttachShader( _programId, geometry_shader->getShaderId() );
        CHECK_GL_ERROR();
    }
    if ( fragment_shader != nullptr ) {
        glAttachShader( _programId, fragment_shader->getShaderId() );
        CHECK_GL_ERROR();
    }
    glLinkProgram( _programId );
    CHECK_GL_ERROR();
    GLint linked;
    glGetProgramiv( _programId, GL_LINK_STATUS, (int *)&linked );
    if ( not linked ) {
        GLint maxLength = 0;
        glGetProgramiv( _programId, GL_INFO_LOG_LENGTH, &maxLength );

        // The maxLength includes the NULL character
        std::vector<GLchar> infoLog( maxLength );
        glGetProgramInfoLog( _programId, maxLength, &maxLength, &infoLog[ 0 ] );
        for ( auto c : infoLog )
            std::cerr << c;
        std::cerr << std::endl;
    }
    return (linked != 0);
}
//==============================================================================
bool Shader::start()
{
    if ( _isModified && _isLinked ) unlink();
    bool success = true;
    if ( !_isLinked )
        success = this->link();
    if ( success ) {
# if defined(GL_TRACE)
        std::cout << "\tok, link succeed " << std::endl;
# endif
        glUseProgramObjectARB( _programId );
        _isLinked = true;
        _isModified = false;
    }
    CHECK_GL_ERROR();
# if defined(GL_TRACE)
    std::cout << "End " << __PRETTY_FUNCTION__ << " with status "
              << std::boolalpha << success << std::endl;
# endif
    return success;;
}
//==============================================================================
std::string Shader::getLinkerLog( void ) const
{
    if ( _programId == 0 ) return std::string( "Not a valid object !" );
    GLint status;
    glGetObjectParameterivARB( _programId, GL_OBJECT_LINK_STATUS_ARB, &status );
    CHECK_GL_ERROR();

    char *linker_log = 0;
    if ( status != 1 ) {
        int infoLgth;
        glGetObjectParameterivARB( _programId, GL_OBJECT_INFO_LOG_LENGTH_ARB, &infoLgth );
        CHECK_GL_ERROR();
        linker_log = new char[ infoLgth ];
        if ( linker_log != 0 ) {
            printf( "ERROR: Could not allocate linker_log buffer\n" );
            return std::string( "Out of memory" );
        }

        int msgLgth;
        glGetInfoLogARB( _programId, infoLgth, &msgLgth, linker_log );
        CHECK_GL_ERROR();
    }

    std::string errMsg( "Unknown linker error" );
    if ( linker_log != 0 )
        errMsg = std::string( linker_log );
    else
        errMsg = std::string( "Linker log is not available !" );
    delete[] linker_log;
    return errMsg;
}
//==============================================================================
void Shader::end()
{
    glUseProgramObjectARB( 0 );
    CHECK_GL_ERROR();
}
// ============================================================================
unsigned int
Shader::getAttributeLocation( const char *attributeName ) const
{
    GLint loc = glGetAttribLocationARB( _programId, attributeName );
    if ( loc == -1 ) throw std::invalid_argument( "Doesn't find attribute variable." );
    CHECK_GL_ERROR();
    return (unsigned int)loc;
}
// ----------------------------------------------------------------------------
unsigned int Shader::getUniformLocation( const char *uniformName ) const
{
    GLint loc = glGetUniformLocationARB( _programId, uniformName );
    if ( loc == -1 ) {
        printf( "Doesn't find uniform variable %s.\n", uniformName );
        throw std::invalid_argument( "Doesn't find uniform variable." );
    }
    CHECK_GL_ERROR();
    return (unsigned int)loc;
}
// ----------------------------------------------------------------------------
unsigned int Shader::getStorageBufferLocation( const char *storageName ) const
{
    GLuint index = glGetProgramResourceIndex( _programId, GL_SHADER_STORAGE_BLOCK,
                                              storageName );

    if ( index == GL_INVALID_INDEX ) {
        printf( "Doesn't find storage buffer variable %s.\n", storageName );
        throw std::invalid_argument( "Doesn't find storage buffer variable." );
    }
    CHECK_GL_ERROR();
    return (unsigned int)index;
}
// ============================================================================
void Shader::setUniform(const char *varName, float v0)
{
    GLint loc = GLint(getUniformLocation(varName));
    glUniform1fARB(loc, v0);
    CHECK_GL_ERROR();
}
// ----------------------------------------------------------------------------
void Shader::setUniform(unsigned int loc, float v0)
{
    glUniform1fARB(GLint(loc), v0);
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform(const char *varName, float v0, float v1)
{
    GLint loc = GLint( getUniformLocation( varName ) );
    glUniform2fARB( loc, v0, v1 );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( unsigned int loc, float v0, float v1 )
{
    glUniform2fARB( GLint( loc ), v0, v1 );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( const char *varName, float v0, float v1, float v2 )
{
    GLint loc = GLint( getUniformLocation( varName ) );
    glUniform3fARB( loc, v0, v1, v2 );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( unsigned int loc, float v0, float v1, float v2 )
{
    glUniform3fARB( GLint( loc ), v0, v1, v2 );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( const char *varName,
                         float v0, float v1, float v2, float v3 )
{
    GLint loc = GLint( getUniformLocation( varName ) );
    glUniform4fARB( loc, v0, v1, v2, v3 );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( unsigned int loc,
                         float v0, float v1, float v2, float v3 )
{
    glUniform4fARB( GLint( loc ), v0, v1, v2, v3 );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform(const char *varName, int v0)
{
    GLint loc = GLint( getUniformLocation( varName ) );
    glUniform1iARB( loc, v0 );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform(unsigned int loc, int v0)
{
    glUniform1iARB( GLint( loc ), v0 );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( const char *varName, int v0, int v1 )
{
    GLint loc = GLint( getUniformLocation( varName ) );
    glUniform2iARB( loc, v0, v1 );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( unsigned int loc, int v0, int v1 )
{
    glUniform2iARB( GLint( loc ), v0, v1 );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( const char *varName, int v0, int v1, int v2 )
{
    GLint loc = GLint( getUniformLocation( varName ) );
    glUniform3iARB( loc, v0, v1, v2 );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( unsigned int loc, int v0, int v1, int v2 )
{
    glUniform3iARB( GLint( loc ), v0, v1, v2 );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( const char *varName,
                         int v0, int v1, int v2, int v3 )
{
    GLint loc = GLint( getUniformLocation( varName ) );
    glUniform4iARB( loc, v0, v1, v2, v3 );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( unsigned int loc, int v0, int v1, int v2, int v3 )
{
    glUniform4iARB( GLint( loc ), v0, v1, v2, v3 );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( const char *varName, size_t sz, size_t count,
                         float *value )
{
    unsigned int loc = getUniformLocation( varName );
    setUniform( loc, sz, count, value );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( unsigned int loc, size_t sz, size_t count, float *value )
{
    switch ( sz ) {
        case 1:
            glUniform1fvARB( GLint( loc ), count, value );
            break;
        case 2:
            glUniform2fvARB( GLint( loc ), count, value );
            break;
        case 3:
            glUniform3fvARB( GLint( loc ), count, value );
            break;
        case 4:
            glUniform4fvARB( GLint( loc ), count, value );
            break;
        default:
            throw std::length_error( "Size of vectors must be 1,2, 3 or 4" );
    }
    _isModified = true;
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( const char *varName, size_t sz, size_t count,
                         int *value )
{
    unsigned int loc = getUniformLocation( varName );
    setUniform( loc, sz, count, value );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniform( unsigned int loc, size_t sz, size_t count, int *value )
{
    switch ( sz ) {
        case 1:
            glUniform1ivARB( GLint( loc ), count, value );
            break;
        case 2:
            glUniform2ivARB( GLint( loc ), count, value );
            break;
        case 3:
            glUniform3ivARB( GLint( loc ), count, value );
            break;
        case 4:
            glUniform4ivARB( GLint( loc ), count, value );
            break;
        default:
            throw std::length_error( "Size of vectors must be 1,2, 3 or 4" );
    }
    _isModified = true;
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniformMatrix( const char *varname, size_t dim, size_t count,
                               bool transpose, float *value )
{
    unsigned int loc = getUniformLocation( varname );
    setUniformMatrix( loc, dim, count, transpose, value );
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
void Shader::setUniformMatrix( unsigned int loc, size_t dim, size_t count,
                               bool transpose, float *value )
{
    switch ( dim ) {
        case 2:
            glUniformMatrix2fvARB( GLint( loc ), count, GLboolean( transpose ), value );
            break;
        case 3:
            glUniformMatrix3fvARB( GLint( loc ), count, GLboolean( transpose ), value );
            break;
        case 4:
            glUniformMatrix4fvARB( GLint( loc ), count, GLboolean( transpose ), value );
            break;
        default:
            throw std::length_error( "Size of matrices must be 2, 3 or 4." );
    }
    _isModified = true;
    CHECK_GL_ERROR();
}
// ============================================================================
void Shader::setAttribute( unsigned int index, float v0 )
{
    glVertexAttrib1fARB( GLint( index ), GLfloat( v0 ) );
    CHECK_GL_ERROR();
}
//==============================================================================
void Shader::setAttribute( unsigned int index, float v0, float v1 )
{
    glVertexAttrib2fARB( GLint( index ), GLfloat( v0 ), GLfloat( v1 ) );
    CHECK_GL_ERROR();
}
//==============================================================================
void Shader::setAttribute( unsigned int index, float v0, float v1, float v2 )
{
    glVertexAttrib3fARB( GLint( index ), GLfloat( v0 ), GLfloat( v1 ), GLfloat( v2 ) );
    CHECK_GL_ERROR();
}
//==============================================================================
void Shader::setAttribute( unsigned int index, float v0, float v1, float v2,
                           float v3 )
{
    glVertexAttrib4fARB( GLint( index ), GLfloat( v0 ), GLfloat( v1 ), GLfloat( v2 ),
                         GLfloat( v3 ) );
    CHECK_GL_ERROR();
}
//==============================================================================
void Shader::setAttribute( unsigned int index, unsigned int dim,
                           const float *v )
{
    switch ( dim ) {
        case 1:
            glVertexAttrib1fvARB( GLint( index ), (GLfloat *)v );
            break;
        case 2:
            glVertexAttrib2fvARB( GLint( index ), (GLfloat *)v );
            break;
        case 3:
            glVertexAttrib3fvARB( GLint( index ), (GLfloat *)v );
            break;
        case 4:
            glVertexAttrib4fvARB( GLint( index ), (GLfloat *)v );
            break;
        default:
            throw std::length_error( "Size of vectors must be 1, 2, 3 or 4." );
    }
    CHECK_GL_ERROR();
}
// ------------------------------------------------------------------------
unsigned int
Shader::initStorageBuffer( int nfloats, const float *values )
{
    GLuint ssbo;
    glGenBuffers( 1, &ssbo );
    glBindBuffer( GL_SHADER_STORAGE_BUFFER, ssbo );
    glBufferData( GL_SHADER_STORAGE_BUFFER, nfloats * sizeof( float ), (void *)values, GL_DYNAMIC_COPY );
    glBindBuffer( GL_SHADER_STORAGE_BUFFER, 0 );
    CHECK_GL_ERROR();
    return (unsigned int)ssbo;
}
// ------------------------------------------------------------------------
void Shader::updateStorageBuffer( unsigned int ssbo, int nfloats, const float *values )
{
    glBindBuffer( GL_SHADER_STORAGE_BUFFER, ssbo );
    GLvoid *p = glMapBuffer( GL_SHADER_STORAGE_BUFFER, GL_WRITE_ONLY );
    memcpy( p, (void *)values, nfloats * sizeof( float ) );
    glUnmapBuffer( GL_SHADER_STORAGE_BUFFER );
    CHECK_GL_ERROR();
}
