cmake_minimum_required(VERSION 3.22)
project(caps_sa_wrapper
        VERSION 0.1.0
        LANGUAGES CXX
)
set(CMAKE_CXX_STANDARD 17)
set(CAPSSA_QUIET TRUE)

add_subdirectory(CaPS-SA)

get_directory_property(CAPSSA_COMPILE_DEFS DIRECTORY CaPS-SA COMPILE_DEFINITIONS)
foreach (d ${CAPSSA_COMPILE_DEFS})
  message("CaPS-SA set compile definition: '${d}'")
  add_compile_definitions(${d})
endforeach()

add_library(caps_sa_wrapper STATIC wrapper.cpp)
target_link_libraries(caps_sa_wrapper PRIVATE core)
install(TARGETS caps_sa_wrapper core DESTINATION .)
