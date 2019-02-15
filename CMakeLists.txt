cmake_minimum_required(VERSION 3.4.2)
project(ptss-dse)

set(CMAKE_BUILD_TYPE Debug)
set(PROJECT_DESCRIPTION "Design Space Exploration and Static Analyses for Probabilistic Timed Sensitive Systems (As a part of ECOLAB NUS)")
set(PROJECT_URL "https://github.com/Arka2009/ptss-dse")
set(PTSS_DSE_VERSION_MAJOR 1)
set(PTSS_DSE_VERSION_MINOR 0)
set(EXEC1 ptssdse)
set(EXEC2 ptssutil)
set(EXEC3 testnlopt)

 
# Targets
add_executable(${EXEC1} src/ptss_dse_hrt.cpp src/ptss_pkmin.cpp src/ptss_test_pkmin.cpp src/ptss_nlopt.cpp)
target_include_directories(${EXEC1} PUBLIC ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR})
target_compile_options(${EXEC1} PRIVATE -std=c++11 -Wall)
target_link_libraries(${EXEC1} PUBLIC -lnlopt -lm)

add_executable(${EXEC2} src/ptss_dse_hrt.cpp src/ptss_nlopt.cpp src/ptss_test.cpp)
target_include_directories(${EXEC2} PUBLIC ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR})
target_compile_options(${EXEC2} PRIVATE -std=c++11 -Wall)
target_link_libraries(${EXEC2} PUBLIC -lnlopt -lm)

add_executable(${EXEC3} src/ptss_nlopt.cpp src/ptss_nlopt_test.cpp)
target_include_directories(${EXEC3} PUBLIC ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR})
target_compile_options(${EXEC3} PRIVATE -std=c++11 -Wall)
target_link_libraries(${EXEC3} PUBLIC -lnlopt -lm)
