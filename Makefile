ifndef BOOST_ROOT
$(error BOOST_ROOT is not set)
endif

USER_CXXFLAGS=-march=native
CXXFLAGS = -std=c++14 -DNDEBUG -O3 $(USER_CXXFLAGS) -flto=auto -fno-stack-protector -isystem $(BOOST_ROOT)/include -MMD -MP -fnon-call-exceptions -fasynchronous-unwind-tables
TEST_OBJECTS_LOCAL=test.o \
	test0.o test1.o test2.o test3.o \
	test4.o test5.o test6.o test7.o test8.o test9.o \
	test10.o test11.o test12.o test13.o test14.o test15.o test16.o \
	test17.o test18.o test19.o test20.o test21.o test22.o test23.o \
	test24.o test25.o test26.o test27.o test28.o test29.o test30.o test31.o \
	test32.o test33.o
TEST_OBJECTS=$(addprefix test/, $(TEST_OBJECTS_LOCAL))
OBJECTS = cjacobian.o ctaylor.o VBIC95Jac/VBIC95Jac.o LUFAC/lufac.o VBIC95/VBIC95.o \
	BLACK_SCHOLES/autodiff_black_scholes.o BLACK_SCHOLES/autodiff_black_scholes_orig.o \
	logistic_regression/logistic_regression.o $(TEST_OBJECTS)

DEPS = $(OBJECTS:.o=.d)

all: ctaylor.exe vbic95Jac.exe vbic95Taylor.exe black_scholes.exe cjacobian.exe black_scholes_orig.exe \
	logistic_regression.exe test/regtest.exe

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test/regtest.exe: $(TEST_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

cjacobian.exe: cjacobian.o
	$(CXX) $(CXXFLAGS) -o $@ $^

ctaylor.exe: ctaylor.o
	$(CXX) $(CXXFLAGS) -o $@ $^

vbic95Jac.exe: VBIC95Jac/VBIC95Jac.o LUFAC/lufac.o
	$(CXX) $(CXXFLAGS) -o $@ $^

vbic95Taylor.exe: VBIC95/VBIC95.o LUFAC/lufac.o
	$(CXX) $(CXXFLAGS) -o $@ $^

black_scholes.exe: BLACK_SCHOLES/autodiff_black_scholes.o
	$(CXX) $(CXXFLAGS) -o $@ $^

black_scholes_orig.exe: BLACK_SCHOLES/autodiff_black_scholes_orig.o
	$(CXX) $(CXXFLAGS) -o $@ $^

logistic_regression.exe: logistic_regression/logistic_regression.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	@find . -type f -name "*.d" | xargs rm -f
	@find . -type f -name "*.o" | xargs rm -f
	@find . -type f -name "*.exe" | xargs rm -f

$(BOOST_ROOT)/include:
	@echo \$$\(BOOST_ROOT\)/include does not exist!
	@exit 1

-include $(DEPS)
