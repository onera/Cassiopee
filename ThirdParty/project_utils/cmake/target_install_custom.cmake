function(target_install_custom project_name path file_ext)
  set(annex_copied_modules_${project_name}_${file_ext}) # creation of the empty list

  # message("project name: \"${project_name}\", file extension:\"${file_ext}\"")

  file(GLOB_RECURSE _annex_files CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/${project_name}/*.${file_ext}")
  foreach (annex_file IN LISTS _annex_files)

    file(RELATIVE_PATH annex_rel_file  ${CMAKE_CURRENT_SOURCE_DIR} ${annex_file})
    set(output_annex_file "${CMAKE_CURRENT_BINARY_DIR}/${annex_rel_file}")

    # message("output_annex_file::" ${output_annex_file})

    add_custom_command(OUTPUT  "${output_annex_file}"
                       DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${annex_rel_file}"
                       COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                       "${CMAKE_CURRENT_SOURCE_DIR}/${annex_rel_file}"
                       "${output_annex_file}"
                       COMMENT "Copying ${annex_rel_file} to the binary directory")

    get_filename_component(annex_file_directory "${annex_rel_file}" DIRECTORY)
    install(FILES       "${annex_rel_file}"
            DESTINATION "${path}/${annex_file_directory}"
            COMPONENT   "annex")
    list(APPEND annex_copied_modules_${project_name}_${file_ext} "${output_annex_file}")
  endforeach ()

  add_custom_target(project_copy_${project_name}_${file_ext} ALL
                    DEPENDS
                    ${annex_copied_modules_${project_name}_${file_ext}})
endfunction()
