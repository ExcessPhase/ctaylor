ifndef BOOST_ROOT
$(error BOOST_ROOT is not set)
endif
all: ctaylor.exe vbic95Jac.exe vbic95Taylor.exe black_scholes.exe cjacobian.exe
#CXX=g++
CFLAGS=-std=c++14 -DNDEBUG -O3 -march=native -flto -isystem $(BOOST_ROOT)/include
%.o: %.cpp $(BOOST_ROOT)/include
	$(CXX) -c $< -o $@ $(CFLAGS)

cjacobian.exe: cjacobian.o
	$(CXX) $(CFLAGS) -o $@ $^

ctaylor.exe: ctaylor.o
	$(CXX) $(CFLAGS) -o $@ $^

vbic95Jac.exe: VBIC95Jac/VBIC95Jac.o LUFAC/lufac.o
	$(CXX) $(CFLAGS) -o $@ $^

vbic95Taylor.exe: VBIC95/VBIC95.o LUFAC/lufac.o
	$(CXX) $(CFLAGS) -o $@ $^

black_scholes.exe: BLACK_SCHOLES/autodiff_black_scholes.o
	$(CXX) $(CFLAGS) -o $@ $^

clean:
	@find . -type f -name "*.o"|xargs rm -f
	@find . -type f -name "*.exe"|xargs rm -f

$(BOOST_ROOT)/include:
	@echo \$$\(BOOST_ROOT\)/include does not exist!
	@exit 1

VBIC95/VBIC95.o: LUFAC/lufac.h \
ctaylor.h \
VBIC95/circuitNodes.h \
VBIC95/currentSources.h \
VBIC95/inputs.h \
VBIC95/members.h \
VBIC95/nodes.h \
VBIC95/parameters.h \
VBIC95/temperatureSetup.h

VBIC95/VBIC95Jac.o: LUFAC/lufac.h \
cjacobian.h \
VBIC95/circuitNodes.h \
VBIC95/currentSources.h \
VBIC95/inputs.h \
VBIC95/members.h \
VBIC95/nodes.h \
VBIC95/parameters.h \
VBIC95/temperatureSetup.h

ctaylor.o:ctaylor.h

cjacobian.o:cjacobian.h

LACK_SCHOLES/autodiff_black_scholes.o:ctaylor.h
