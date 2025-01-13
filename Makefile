ifndef BOOST_ROOT
$(error BOOST_ROOT is not set)
endif

USER_CXXFLAGS=-march=native
CXXFLAGS = -std=c++14 -DNDEBUG -O3 $(USER_CXXFLAGS) -flto=auto -fno-stack-protector -isystem $(BOOST_ROOT)/include -MMD -MP
OBJECTS = cjacobian.o ctaylor.o VBIC95Jac/VBIC95Jac.o LUFAC/lufac.o VBIC95/VBIC95.o \
	BLACK_SCHOLES/autodiff_black_scholes.o BLACK_SCHOLES/autodiff_black_scholes_orig.o \
	logistic_regression/logistic_regression.o

DEPS = $(OBJECTS:.o=.d)

all: ctaylor.exe vbic95Jac.exe vbic95Taylor.exe black_scholes.exe cjacobian.exe black_scholes_orig.exe \
	logistic_regression.exe test/test.exe

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test/test.exe: test/test.o
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
