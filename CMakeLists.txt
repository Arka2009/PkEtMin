cmake_minimum_required(VERSION 3.4.2)
project(ptss-dse)

set(CMAKE_BUILD_TYPE Debug)
set(PROJECT_DESCRIPTION "Design Space Exploration and Static Analyses for Probabilistic Timed Sensitive Systems (As a part of ECOLAB NUS)")
set(PROJECT_URL "https://github.com/Arka2009/ptss-dse")
set(PTSS_DSE_VERSION_MAJOR 1)
set(PTSS_DSE_VERSION_MINOR 0)
set(EXEC ptssdse)

 
# Target
add_executable(${EXEC} src/ptss_dse_hrt.cpp src/ptss_pkmin.cpp)
target_include_directories(${EXEC} PUBLIC ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR})
target_compile_options(${EXEC} PRIVATE -std=c++11 -Wall)
target_link_libraries(${EXEC} PUBLIC -lm)