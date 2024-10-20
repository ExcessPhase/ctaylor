ifndef BOOST_ROOT
$(error BOOST_ROOT is not set)
endif
all: ctaylor.exe vbic95Jac.exe vbic95Taylor.exe black_scholes.exe cjacobian.exe
CC=g++
CFLAGS=-std=c++14 -DNDEBUG -O3 -march=native -flto -isystem $(BOOST_ROOT)/include

cjacobian.exe: cjacobian.cpp
	$(CC) $(CFLAGS) -o cjacobian.exe cjacobian.cpp

ctaylor.exe: ctaylor.cpp
	$(CC) $(CFLAGS) -o ctaylor.exe ctaylor.cpp

vbic95Jac.exe: VBIC95/VBIC95.cpp LUFAC/lufac.cpp
	$(CC) $(CFLAGS) -o vbic95Jac.exe VBIC95/VBIC95.cpp LUFAC/lufac.cpp -D__JACOBIAN__

vbic95Taylor.exe: VBIC95/VBIC95.cpp LUFAC/lufac.cpp
	$(CC) $(CFLAGS) -o vbic95Taylor.exe VBIC95/VBIC95.cpp LUFAC/lufac.cpp

black_scholes.exe: BLACK_SCHOLES/autodiff_black_scholes.cpp
	$(CC) $(CFLAGS) -o black_scholes.exe BLACK_SCHOLES/autodiff_black_scholes.cpp
