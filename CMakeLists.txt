cmake_minimum_required(VERSION 3.8)

#cmake_policy(SET CMP0022 OLD)
project(PostCutCodeAnalysis)

set(CMAKE_BUILD_TYPE Debug)

# Recurse into the "Hello" and "Demo" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory(PostCutCode)
add_subdirectory(WSMaker)
add_subdirectory(wscanner)
