# Use the Aurora compiler.
#
# This code sets the following variables:
#
#  NEC_FORTRAN_COMPILER
#
# See also FindAURORA.cmake
#=============================================================================

function( compile_nec _name generated_file)
  # Default to assuming all files are Fortran.
  set( nec_object )

  file(GLOB_RECURSE fortran_sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/*.f90 ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/*.for)
  file(GLOB_RECURSE c_sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/*.c)
  file(GLOB_RECURSE cpp_sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/utilities/*.nec.cpp)

  list(FILTER fortran_sources EXCLUDE REGEX ".*\.in\.for$")
  list(FILTER fortran_sources EXCLUDE REGEX ".*\.in\.f90$")

  foreach( fortran_file ${fortran_sources} )
    # message( "fortran_file : " ${fortran_file} )

    # get_filename_component(python_file_directory "${python_file}" DIRECTORY)
    get_filename_component(fortran_file_name_we "${fortran_file}" NAME_WE)
    get_filename_component(fortran_type "${fortran_file}" LAST_EXT)

    file(RELATIVE_PATH _fortran_rel ${CMAKE_CURRENT_SOURCE_DIR} ${fortran_file})

    get_filename_component(_fortran_dir_rel "${_fortran_rel}" DIRECTORY)

    # message( "fortran_type : " ${fortran_type} )
    # message( "_fortran_rel : " ${_fortran_rel} )
    # message( "_fortran_dir_rel : " ${_fortran_dir_rel} )

    set( _pre_generated_file "${CMAKE_CURRENT_BINARY_DIR}/${_fortran_dir_rel}/${fortran_file_name_we}.fpp${fortran_type}" )
    set_source_files_properties( ${_pre_generated_file} PROPERTIES GENERATED TRUE )


    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_fortran_dir_rel})
    # message("Create directory : " ${CMAKE_CURRENT_BINARY_DIR}/${_fortran_dir_rel})

    # message( "_generated_file : " ${_generated_file} )
    # message( "_pre_generated_file : " ${_pre_generated_file} )
    set( comment "Pre-Processing Fortran source for ${fortran_file_name_we}..." )
    add_custom_command( OUTPUT  ${_pre_generated_file}
                        COMMAND cpp
                        ARGS    ${fortran_file} -DE_DOUBLEREAL -D__ve__ -I/${PROJECT_SOURCE_DIR} -o ${_pre_generated_file}
                        DEPENDS ${fortran_file}
                        COMMENT ${comment_preprocessor}
    )
    # Si 2 fichier on le mÃªme nom --> Marche pas
    # Manage same name of file in different module
    string(REPLACE "/" "_" flat_name_nec_pre ${_pre_generated_file} )
    # add_custom_target(nec_pre_${fortran_file_name_we} ALL DEPENDS "${_pre_generated_file}")
    add_custom_target(${flat_name_nec_pre} ALL DEPENDS "${_pre_generated_file}")


    set( _generated_file "${CMAKE_CURRENT_BINARY_DIR}/${_fortran_dir_rel}/${fortran_file_name_we}_nec.o" )
    set_source_files_properties( ${_generated_file} PROPERTIES GENERATED TRUE )
    # set( ${generated_file} ${_generated_file} PARENT_SCOPE )


    # message( "_generated_file : " ${_generated_file} )
    set( comment "Compiling Aurora Fortran source for ${fortran_file_name_we}..." )
    add_custom_command( OUTPUT  ${_generated_file}
                        COMMAND ${NEC_FORTRAN_COMPILER}
                        ARGS    ${_pre_generated_file} -fpp -fdiag-vector=2 -fopenmp -fPIC -traceback -DE_DOUBLEREAL -fextend-source -report-all -c -I/${PROJECT_SOURCE_DIR} -o ${_generated_file}
                        DEPENDS ${_pre_generated_file}
                        COMMENT ${comment}
    )
    string(REPLACE "/" "_" flat_name_nec ${_generated_file} )

    # message( "flat_name_nec : " ${flat_name_nec} )
    add_custom_target(${flat_name_nec} ALL DEPENDS "${_generated_file}")
    # add_custom_target(nec_obj_${fortran_file_name_we} ALL DEPENDS "${_generated_file}")

    list( APPEND nec_object "${_generated_file}" )

  endforeach()

  # Manage C Files
  foreach( c_file ${c_sources} )
    # message( "c_file : " ${c_file} )

    # get_filename_component(python_file_directory "${python_file}" DIRECTORY)
    get_filename_component(c_file_name_we "${c_file}" NAME_WE)

    file(RELATIVE_PATH _c_rel ${CMAKE_CURRENT_SOURCE_DIR} ${c_file})
    get_filename_component(_c_dir_rel "${_c_rel}" DIRECTORY)

    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_c_dir_rel})


    set( _generated_file "${CMAKE_CURRENT_BINARY_DIR}/${_c_dir_rel}/${c_file_name_we}_nec.o" )
    set_source_files_properties( ${_generated_file} PROPERTIES GENERATED TRUE )
    # set( ${generated_file} ${_generated_file} PARENT_SCOPE )


    # message( "_generated_file : " ${_generated_file} )
    set( comment "Compiling Aurora C source for ${c_file_name_we}..." )
    add_custom_command( OUTPUT  ${_generated_file}
                        COMMAND ${NEC_C_COMPILER}
                        ARGS    ${c_file} -fopenmp -fPIC -I/${PROJECT_SOURCE_DIR} -c -o ${_generated_file}
                        DEPENDS ${c_file}
                        COMMENT ${comment}
    )
    string(REPLACE "/" "_" flat_name_nec_pre ${_generated_file} )
    # add_custom_target(nec_obj_${c_file_name_we} ALL DEPENDS "${_generated_file}")
    add_custom_target(${flat_name_nec_pre} ALL DEPENDS "${_generated_file}")

    list( APPEND nec_object "${_generated_file}" )

  endforeach()

  # Manage CPP Files
  foreach( cpp_file ${cpp_sources} )
    message( "cpp_file : " ${cpp_file} )

    # get_filename_component(python_file_directory "${python_file}" DIRECTORY)
    get_filename_component(cpp_file_name_we "${cpp_file}" NAME_WE)

    file(RELATIVE_PATH _cpp_rel ${CMAKE_CURRENT_SOURCE_DIR} ${cpp_file})
    get_filename_component(_cpp_dir_rel "${_cpp_rel}" DIRECTORY)

    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_cpp_dir_rel})


    set( _generated_file "${CMAKE_CURRENT_BINARY_DIR}/${_cpp_dir_rel}/${cpp_file_name_we}.nec_nec.o" )
    set_source_files_properties( ${_generated_file} PROPERTIES GENERATED TRUE )
    # set( ${generated_file} ${_generated_file} PARENT_SCOPE )


    # message( "_generated_file : " ${_generated_file} )
    set( comment "Compiling Aurora C++ source for ${cpp_file_name_we}..." )
    add_custom_command( OUTPUT  ${_generated_file}
                        COMMAND ${NEC_CXX_COMPILER}
                        ARGS    ${cpp_file} -std=c++17 -fopenmp -fPIC -I/${PROJECT_SOURCE_DIR} -c -o ${_generated_file}
                        DEPENDS ${cpp_file}
                        COMMENT ${comment}
    )
    string(REPLACE "/" "_" flat_name_nec_pre ${_generated_file} )
    # add_custom_target(nec_obj_${c_file_name_we} ALL DEPENDS "${_generated_file}")
    add_custom_target(${flat_name_nec_pre} ALL DEPENDS "${_generated_file}")

    list( APPEND nec_object "${_generated_file}" )

  endforeach()


  # message("listall : " ${nec_object})
  # foreach( g  ${nec_object} )
  #   message( "oooa : " ${g})
  # endforeach()

  list(LENGTH nec_object nobject)

  # message( "nobject : " ${nobject})

  # Create Shared library / TODO Static
  if( ${nobject} GREATER 0 )
    set( commentlink "Linking all object into Shared Library source for ${_name}..." )
    add_custom_command( OUTPUT  lib${_name}.so
                        COMMAND ${NEC_C_COMPILER}
                        ARGS    -fopenmp  -shared -o lib${_name}.so ${nec_object}
                        DEPENDS ${nec_object}
                        COMMENT ${commentlink}
    )
    add_custom_target(${_name} ALL DEPENDS "lib${_name}.so")

    install(FILES       "${CMAKE_CURRENT_BINARY_DIR}/lib${_name}.so"
            DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

    # > Tentative non fructueuse
    # add_executable( ${_name} IMPORTED GLOBAL "lib${_name}.so")
    # install(TARGETS ${_name}
    #         DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

  endif()

  # install(TARGETS "${mod_name}"
  #         LIBRARY DESTINATION ${SITE_PACKAGES_OUTPUT_DIRECTORY}/${rel}/${_name})

endfunction()



function( nec_add_module _name )

  compile_nec( ${_name} generated_file )

endfunction()
