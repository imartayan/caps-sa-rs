cmake_minimum_required(VERSION 3.22)
project(caps_sa_wrapper
        VERSION 0.1.0
        LANGUAGES CXX
)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(CaPS-SA)

get_directory_property(CAPSSA_COMPILE_DEFS DIRECTORY CaPS-SA COMPILE_DEFINITIONS)
list(GET CAPSSA_COMPILE_DEFS 1 L1_CACHE_LINE_SIZE_VAR)
string(REPLACE "L1_CACHE_LINE_SIZE=" "" L1_CACHE_LINE_SIZE ${L1_CACHE_LINE_SIZE_VAR})
message("L1_CACHE_LINE_SIZE = ${L1_CACHE_LINE_SIZE}")
message("CAPSSA_COMPILE_DEFS = ${CAPSSA_COMPILE_DEFS}")
add_compile_definitions(L1_CACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE})

add_library(caps_sa_wrapper STATIC wrapper.cpp)
target_link_libraries(caps_sa_wrapper PRIVATE core)
install(TARGETS caps_sa_wrapper core DESTINATION .)
