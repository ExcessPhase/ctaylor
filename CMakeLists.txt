cmake_minimum_required(VERSION 3.10)
project(ctaylor VERSION 1.0 LANGUAGES CXX)

# Ensure C++11 standard or higher is used
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the Boost libraries
find_package(Boost REQUIRED COMPONENTS filesystem system)

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost libraries not found!")
endif()

# Include Boost headers
include_directories(${Boost_INCLUDE_DIRS})

# Source files
set(TEST_SOURCES
    test/test.cpp
    test/test0.cpp
    test/test1.cpp
    test/test2.cpp
    test/test3.cpp
    test/test4.cpp
    test/test5.cpp
    test/test6.cpp
    test/test7.cpp
    test/test8.cpp
    test/test9.cpp
    test/test10.cpp
    test/test11.cpp
    test/test12.cpp
    test/test13.cpp
    test/test14.cpp
    test/test15.cpp
    test/test16.cpp
    test/test17.cpp
    test/test18.cpp
    test/test19.cpp
    test/test20.cpp
    test/test21.cpp
    test/test22.cpp
    test/test23.cpp
    test/test24.cpp
    test/test25.cpp
    test/test26.cpp
    test/test27.cpp
    test/test28.cpp
    test/test29.cpp
    test/test30.cpp
)

# Executable targets
add_executable(ctaylor ctaylor.cpp)
add_executable(vbic95Jac VBIC95Jac/VBIC95Jac.cpp LUFAC/lufac.cpp)
add_executable(vbic95Taylor VBIC95/VBIC95.cpp LUFAC/lufac.cpp)
add_executable(black_scholes BLACK_SCHOLES/autodiff_black_scholes.cpp)
add_executable(black_scholes_orig BLACK_SCHOLES/autodiff_black_scholes_orig.cpp)
add_executable(logistic_regression logistic_regression/logistic_regression.cpp)
add_executable(test_exe ${TEST_SOURCES})

# Include directories
target_include_directories(ctaylor PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(vbic95Jac PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(vbic95Taylor PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(black_scholes PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(black_scholes_orig PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(logistic_regression PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(test_exe PRIVATE ${Boost_INCLUDE_DIRS})

# Link Boost libraries
target_link_libraries(ctaylor ${Boost_LIBRARIES})
target_link_libraries(vbic95Jac ${Boost_LIBRARIES})
target_link_libraries(vbic95Taylor ${Boost_LIBRARIES})
target_link_libraries(black_scholes ${Boost_LIBRARIES})
target_link_libraries(black_scholes_orig ${Boost_LIBRARIES})
target_link_libraries(logistic_regression ${Boost_LIBRARIES})
target_link_libraries(test_exe ${Boost_LIBRARIES})

# Optionally, add include directories for your own headers
target_include_directories(ctaylor PRIVATE include)
# Repeat for other targets as needed

# Enable testing
enable_testing()

# Add tests if any (adjust accordingly)
add_test(NAME ctaylor_test COMMAND ctaylor)
