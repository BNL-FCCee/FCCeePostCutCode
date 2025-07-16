cmake_minimum_required(VERSION 3.8)

project(PostCutCode)
message("announcing PostCutCode")
set(CMAKE_BUILD_TYPE Debug)

add_compile_options(-std=c++17 -fPIC -Wall -Wextra -Wno-unused-parameter -O2)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS RIO EG Minuit RooFit RooFitCore RooStats HistFactory Core Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread MultiProc TMVA)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Collect source files
file(GLOB sources "Root/*.cxx")

# Generate the dictionary
set(DICT_HEADERS PostCutCode/RVecDict.h)
set(DICT_LINKDEF PostCutCode/LinkDef.h)

ROOT_GENERATE_DICTIONARY(RVecDict
    ${DICT_HEADERS}
    LINKDEF ${DICT_LINKDEF}
)

# Add the shared library and include the dictionary source variable
add_library(PostCutCodeLib SHARED ${sources} ${RVecDict_dictionary})
target_include_directories(PostCutCodeLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(PostCutCodeLib ${ROOT_LIBRARIES})

# Build executables
file(GLOB_RECURSE exe_files "utils/*.cxx")
foreach(_exe_file ${exe_files})
    get_filename_component(_theExec ${_exe_file} NAME_WE)
    add_executable(${_theExec} ${_exe_file})
    target_link_libraries(${_theExec} PostCutCodeLib)
endforeach()

# Setup script
set(SETUP ${CMAKE_CURRENT_BINARY_DIR}/../setup.sh)
file(WRITE ${SETUP} "#!/bin/bash\n")
file(APPEND ${SETUP} "# this is an auto-generated setup script\n")
file(APPEND ${SETUP} "export PATH=\${PATH}:${CMAKE_CURRENT_BINARY_DIR}/\n")

# Optional: export data env
set(ENV{WSCANNER_DATADIR} ${CMAKE_CURRENT_SOURCE_DIR}/data/)

