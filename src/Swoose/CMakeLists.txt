# Load lists of source files
include(Files.cmake)

# Find / Import dependecies
include(ImportUtilsOS)
include(ImportCore)
include(ImportMolassembler)
include(ImportYamlCpp)
include(AddEigen)

# Only include database and mongocxx if required
if(SWOOSE_COMPILE_DATABASE)
  include(ImportDatabase)
  include(ImportMongoCxx)
endif()

# Do not build tests of dependencies
set(_build_tests ${SCINE_BUILD_TESTS})
set(SCINE_BUILD_TESTS OFF)
import_utils_os()
import_core()
import_molassembler()
if(SWOOSE_COMPILE_DATABASE)
  import_database()
  import_mongocxx()
endif()
import_yamlcpp()
set(SCINE_BUILD_TESTS ${_build_tests})

find_package(Boost COMPONENTS program_options REQUIRED)
find_package(OpenMP)

# Shared Swoose library / Swoose module
add_library(Swoose SHARED ${SWOOSE_MODULE_FILES})
set_property(TARGET Swoose PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(Swoose PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

set_target_properties(Swoose PROPERTIES
  PREFIX ""
  OUTPUT_NAME swoose.module
)

target_compile_definitions(Swoose PRIVATE $<$<CONFIG:Debug>:EIGEN_INITIALIZE_MATRICES_BY_NAN>)
if(SWOOSE_COMPILE_DATABASE)
  target_compile_definitions(Swoose PUBLIC SWOOSE_COMPILE_DATABASE)
endif()

target_compile_options(Swoose PRIVATE $<$<BOOL:${OpenMP_CXX_FOUND}>:${OpenMP_CXX_FLAGS}>)

target_link_libraries(Swoose
  PUBLIC
    Scine::Core
    Scine::UtilsOS
    $<$<BOOL:${SWOOSE_COMPILE_DATABASE}>:Scine::Database>
  PRIVATE
    Scine::Molassembler
    yaml-cpp
    $<$<BOOL:${SWOOSE_COMPILE_DATABASE}>:MongoDBCXX>
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
  )

add_eigen(Swoose PUBLIC)

# Add namespaced aliases
add_library(Scine::Swoose ALIAS Swoose)

#-- Define Install
# Headers
install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/Swoose/
  DESTINATION include/Swoose
  FILES_MATCHING PATTERN "*.h"
)

# Targets
install(
  TARGETS Swoose
  EXPORT swooseTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include/Swoose
)

# Component - level cmake package files
scine_install_component_cmake_files(
  COMPONENT  ScineSwoose
  EXPORT_NAME swooseTargets
)

# Set the RPATH to be used when installing.$ORIGIN refers to the installation
# directory of the binary.
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")

# Executables

add_executable(SwooseApp ${SWOOSE_APP_FILES})
add_executable(Scine::SwooseApp ALIAS SwooseApp)

set_target_properties(SwooseApp PROPERTIES OUTPUT_NAME swoose)

target_link_libraries(SwooseApp PRIVATE
  Boost::program_options
  Scine::Swoose
  Scine::UtilsOS
  Scine::Core
  yaml-cpp
)

# Targets
install(TARGETS SwooseApp RUNTIME DESTINATION bin)

# Tests
if(SCINE_BUILD_TESTS)
  include(ImportGTest)
  import_gtest()
  add_subdirectory(Tests)
  add_executable(Swoose_tests
    ${SWOOSE_TEST_FILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/App/AppUtils/Tasks.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/App/AppUtils/Tasks.h
  )
  target_link_libraries(Swoose_tests PRIVATE
    GTest::Main
    GMock::GMock
    $<$<BOOL:${SWOOSE_COMPILE_DATABASE}>:Scine::Database>
    Scine::Swoose
    Boost::program_options
    yaml-cpp
    Scine::Molassembler
  )
  add_test(NAME Swoose COMMAND Swoose_tests)
  # Hard - code location of files into the test by configuring a header file for it.
  set(Tests_dir "${CMAKE_CURRENT_SOURCE_DIR}/Tests/Files/")
  configure_file(${Tests_dir}/tests_file_location.h.in
                 ${Tests_dir}/tests_file_location.h)
endif()

# Python bindings
if(SCINE_BUILD_PYTHON_BINDINGS)
  include(ImportPybind11)
  import_pybind11()
  set(PYBIND11_PYTHON_VERSION ${PYTHONVERSION})

  # Python module
  pybind11_add_module(scine_swoose
    ${SWOOSE_PYTHON_FILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/App/AppUtils/Tasks.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/App/AppUtils/Tasks.h
    ${CMAKE_CURRENT_SOURCE_DIR}/App/AppUtils/TaskManagement.h
  )
  set_target_properties(scine_swoose PROPERTIES
    SUFFIX ".so"
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/scine_swoose
  )
  target_include_directories(scine_swoose PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/App
    ${YAML_CPP_INCLUDE_DIR}
  )
  target_link_libraries(scine_swoose PRIVATE
    Scine::UtilsOS
    Scine::Core
    Scine::Swoose
    yaml-cpp
  )
  install(CODE
    "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} --upgrade --no-deps ${CMAKE_CURRENT_BINARY_DIR}
    RESULT_VARIABLE retcode)
  if(NOT \$\{retcode\} EQUAL 0)
    message(FATAL_ERROR \"Fatal error when installing Python module using PIP.\")
  endif()"
  )

  # Copy the setup file and README.md into the binary dir so you can install
  # using pip or make a distribution with the compiled binary
  file(
    COPY ${CMAKE_CURRENT_SOURCE_DIR}/Python/__init__.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_swoose
  )
  add_custom_command(
    TARGET Swoose POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Scine::Swoose> $<TARGET_FILE:Scine::UtilsOS> $<TARGET_FILE:Scine::Core> ${CMAKE_CURRENT_BINARY_DIR}/scine_swoose
  )

  include(TargetLibName)
  target_lib_filename(Scine::Core _core_name)
  if (BUILD_SHARED_LIBS)
    if(APPLE)
      set_target_properties(Swoose PROPERTIES
        INSTALL_RPATH "@loader_path;@loader_path/../lib"
      )
    elseif(UNIX)
      set_target_properties(Swoose PROPERTIES
        INSTALL_RPATH "\$ORIGIN;\$ORIGIN/../lib"
      )
    endif()
    add_custom_command(
      TARGET Swoose POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Swoose> $<TARGET_FILE:Scine::UtilsOS> ${CMAKE_CURRENT_BINARY_DIR}/scine_swoose
      COMMENT "Copying dependent shared libraries into python package directory"
    )
    target_lib_filename(Scine::UtilsOS _utils_name)
    target_lib_filename(Swoose _swoose_name)
    message(STATUS "deduced library names: ${_core_name}, ${_utils_name}, ${_swoose_name}")
    set(swoose_PY_DEPS ", \"${_core_name}\", \"${_utils_name}\", \"${_swoose_name}\"")
    unset(_utils_name)
    unset(_swoose_name)
  else()
    message(STATUS "deduced library names: ${_core_name}")
    set(swoose_PY_DEPS ", \"${_core_name}\"")
  endif()
  # Add setuptools file
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Python/setup.py
    ${CMAKE_CURRENT_BINARY_DIR}/setup.py
    @ONLY
  )
  unset(_core_name)
  unset(swoose_PY_DEPS)
  file(
    COPY
    ${PROJECT_SOURCE_DIR}/README.rst
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )
endif() # endif python bindings

# if (WIN32)
#  set_target_properties(Swoose PROPERTIES PREFIX "")
#  set_target_properties(Swoose PROPERTIES IMPORT_PREFIX "")
#  set_target_properties(Swoose PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
# endif()
#
## Add link library specific for MinGW
# if (WIN32)
#  target_link_libraries(UtilsOS PRIVATE $ <$ <BOOL : MINGW> : ws2_32>)
# endif()
