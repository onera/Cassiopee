#ifndef _CPLOT_SHADERS_SHADERBUFFEROBJECT_HPP_
#define _CPLOT_SHADERS_SHADERBUFFEROBJECT_HPP_
#include "GL/glew.h"

/**
 * @brief Memory management for Shader Buffer object
 * @details This class manages the allocation/deallocation of shaders buffer objects
 * which are memory arrays inside the GPGPU. This buffer can be used as data arrays for
 * shaders ( geometry in particular )
 * 
 */
class ShaderBufferObject
{
public:
    ShaderBufferObject() : m_identity{0}
    {}

    template<typename K>
    ShaderBufferObject( unsigned long nsize, const K* data, GLenum usage = GL_STATIC_DRAW )
    {
        glGenBuffers(1, &m_identity);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_identity);
        glBufferData(GL_SHADER_STORAGE_BUFFER, nsize*sizeof(K), (void*)data, usage);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    }

    ShaderBufferObject( const ShaderBufferObject& ) = delete;
    ShaderBufferObject( ShaderBufferObject&& ssbo ) = default;
    ~ShaderBufferObject() { glDeleteBuffers(1, &m_identity); }

    ShaderBufferObject& operator = ( const ShaderBufferObject& ) = delete;
    ShaderBufferObject& operator = ( ShaderBufferObject&& ) = default;

    bool is_allocated() const {
        GLboolean is_buffer = glIsBuffer(this->m_identity);
        return is_buffer == GL_TRUE;
    }

    static unsigned max_number_of_ssbo_per_shader() {
        GLint max_ssbo; 
        glGetIntegerv(GL_MAX_COMPUTE_SHADER_STORAGE_BLOCKS, &max_ssbo);
        return unsigned(max_ssbo);
    }

    static unsigned max_number_of_ssbo_for_all_shaders() {
        GLint max_ssbo; 
        glGetIntegerv(GL_MAX_COMBINED_SHADER_STORAGE_BLOCKS, &max_ssbo);
        return unsigned(max_ssbo);        
    }

    static unsigned max_number_of_ssbo_for_geometric_shader() {
        GLint max_ssbo; 
        glGetIntegerv(GL_MAX_GEOMETRY_SHADER_STORAGE_BLOCKS, &max_ssbo);
        return unsigned(max_ssbo);
    }

    static unsigned long max_storage_size() {
        GLint max_ssbo; 
        glGetIntegerv(GL_MAX_SHADER_STORAGE_BLOCK_SIZE, &max_ssbo);
        return unsigned(max_ssbo);
    }

private:
    GLuint m_identity;
};
#endif
