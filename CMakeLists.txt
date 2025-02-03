cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Check for Boost
if(NOT DEFINED ENV{BOOST_ROOT})
    message(FATAL_ERROR "BOOST_ROOT is not set")
else()
    set(BOOST_ROOT $ENV{BOOST_ROOT})
endif()

find_package(Boost REQUIRED)
include_directories(${BOOST_ROOT}/include)

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
add_executable(test_exe ${TEST_SOURCES})  # Renaming the target to avoid issues

# Custom clean target
add_custom_target(clean_all
    COMMAND ${CMAKE_COMMAND} -E rm -f $(OBJECTS)
    COMMAND ${CMAKE_COMMAND} -E rm -f ctaylor.exe vbic95Jac.exe vbic95Taylor.exe black_scholes.exe cjacobian.exe black_scholes_orig.exe logistic_regression.exe test_exe
)

# Include directories and libraries
target_include_directories(ctaylor PRIVATE ${BOOST_ROOT}/include)
target_include_directories(vbic95Jac PRIVATE ${BOOST_ROOT}/include)
target_include_directories(vbic95Taylor PRIVATE ${BOOST_ROOT}/include)
target_include_directories(black_scholes PRIVATE ${BOOST_ROOT}/include)
target_include_directories(black_scholes_orig PRIVATE ${BOOST_ROOT}/include)
target_include_directories(logistic_regression PRIVATE ${BOOST_ROOT}/include)
target_include_directories(test_exe PRIVATE ${BOOST_ROOT}/include)  # Updated target name

target_link_libraries(ctaylor ${Boost_LIBRARIES})
target_link_libraries(vbic95Jac ${Boost_LIBRARIES})
target_link_libraries(vbic95Taylor ${Boost_LIBRARIES})
target_link_libraries(black_scholes ${Boost_LIBRARIES})
target_link_libraries(black_scholes_orig ${Boost_LIBRARIES})
target_link_libraries(logistic_regression ${Boost_LIBRARIES})
target_link_libraries(test_exe ${Boost_LIBRARIES})  # Updated target name
