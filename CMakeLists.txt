# Tsung-Yi's first CMakeLists for FixedPointLDPC project. 
# CMakeLists files in this project can 
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} 
# and to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6) 
project (FixedPointLDPC) 

# Recurse into the "Hello" and "Demo" subdirectories. This does not actually 
# cause another cmake executable to run. The same process will walk through 
# the project's entire directory structure. 
#add_subdirectory (FPLDPC) 
#add_subdirectory (Demo)
# --- For now let's just put everything in the same directory

add_library (fpldpc Memory.cxx ArrayLDPC.cxx 

include_directories("${PROJECT_BINARY_DIR}"))


# To be continued 
