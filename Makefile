ifndef BOOST_ROOT
$(error BOOST_ROOT is not set)
endif
all: ctaylor.exe vbic95Jac.exe vbic95Taylor.exe black_scholes.exe cjacobian.exe black_scholes_orig.exe
#CXX=g++
CXXFLAGS=-std=c++14 -DNDEBUG -O3 -march=native -flto=auto -isystem $(BOOST_ROOT)/include -MMD -MP
OBJECTS=cjacobian.o ctaylor.o VBIC95Jac/VBIC95Jac.o LUFAC/lufac.o VBIC95/VBIC95.o \
	BLACK_SCHOLES/autodiff_black_scholes.o BLACK_SCHOLES/autodiff_black_scholes_orig.o
DEPS=$(OBJECTS:.o.d)
%.o %.d: %.cpp $(BOOST_ROOT)/include
	$(CXX) -c $< -o $@ $(CXXFLAGS)

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

clean:
	@find . -type f -name "*.d"|xargs rm -f
	@find . -type f -name "*.o"|xargs rm -f
	@find . -type f -name "*.exe"|xargs rm -f

$(BOOST_ROOT)/include:
	@echo \$$\(BOOST_ROOT\)/include does not exist!
	@exit 1

-include $(DEPS)
