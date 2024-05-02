# CMake commands #

While CMake has improved to become more of a dependency-based build tool, not everything works from scratch.

In particular, a `*Config.cmake` file should include a list of all its dependencies (through calls to `find_package`). But the fact is that these dependencies have already been mentionned, either through `find_package` or `add_subdirectory`, in the `CMakeLists.txt` (because they are needed for the build). The commands provided in `target_add_dependency.cmake` aim to replace `find_package`, `add_subdirectory` and custom code for the generation of `*Config.cmake`, by respectively `target_add_thirdparty_dependency`, `target_add_dependency` and `target_install`.

`target_add_thirdparty_dependency` and `target_add_dependency` are just calls to `find_package` and `add_subdirectory` under the hood, but they also record that the dependency is needed for the target, so that once `target_install` is called, it can automatically generate the dependencies again for the `*Config.cmake` file.
