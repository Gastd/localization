localization
============

**Description**

Dependencies
------------
<!-- * Only [pthreads](https://en.wikipedia.org/wiki/POSIX_Threads), that is usually installed within any Linux system. -->
* [GMatrix](https://github.com/lara-unb/gMatrix)
* [GDataLogger](https://github.com/lara-unb/gdatalogger)

Usage
-----

A CMake code is provided in hopes of help the link process.
Link this library against your project using the code below in your CMakeLists.txt file.

```bash
include(ExternalProject)
ExternalProject_Add(localization
    GIT_REPOSITORY https://github.com/Gastd/localization
    GIT_TAG master
    # PREFIX ${CMAKE_CURRENT_BINARY_DIR}/localization
    SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/localization-src"
    BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/localization-src/build"
    # CONFIGURE_COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" "${CMAKE_BINARY_DIR}/localization-src/"
    BUILD_COMMAND "${CMAKE_COMMAND}" --build .
    INSTALL_COMMAND ""
    TEST_COMMAND ""
)
# ExternalProject_Get_Property(localization install_dir)
ExternalProject_Get_Property(localization binary_dir)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/localization-src/cmex/lib)

add_executable( <YOUR_EXECUTABLE> <YOUR_CODE>.cpp )
add_dependencies( <YOUR_EXECUTABLE> localization )
target_link_libraries( <YOUR_EXECUTABLE>
  ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}localization.a
)
```

Documentation
-------------

Documentation (and further explanations) for the code is available in the source code files.
