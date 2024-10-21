ifndef BOOST_ROOT
$(error BOOST_ROOT is not set)
endif
all: ctaylor.exe vbic95Jac.exe vbic95Taylor.exe black_scholes.exe cjacobian.exe
CC=g++
CFLAGS=-std=c++14 -DNDEBUG -O3 -march=native -flto -isystem $(BOOST_ROOT)/include
%.o: %.cpp
	$(CXX) -c $< -o $@ $(CFLAGS)

cjacobian.exe: cjacobian.o
	$(CC) $(CFLAGS) -o $@ $^

ctaylor.exe: ctaylor.o
	$(CC) $(CFLAGS) -o $@ $^

vbic95Jac.exe: VBIC95Jac/VBIC95Jac.o LUFAC/lufac.o
	$(CC) $(CFLAGS) -o $@ $^

vbic95Taylor.exe: VBIC95/VBIC95.o LUFAC/lufac.o
	$(CC) $(CFLAGS) -o $@ $^

black_scholes.exe: BLACK_SCHOLES/autodiff_black_scholes.o
	$(CC) $(CFLAGS) -o $@ $^
