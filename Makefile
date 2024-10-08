all: ctaylor.exe vbic95.exe black_scholes.exe
CC=g++
CFLAGS=-std=c++14 -DNDEBUG -O3 -march=native -flto -isystem $(BOOST_ROOT)/include -ffinite-math-only

ctaylor.exe: ctaylor.cpp
	$(CC) $(CFLAGS) -o ctaylor.exe ctaylor.cpp

vbic95.exe: VBIC95/VBIC95.cpp LUFAC/lufac.cpp
	$(CC) $(CFLAGS) -o vbic95.exe VBIC95/VBIC95.cpp LUFAC/lufac.cpp

black_scholes.exe: BLACK_SCHOLES/autodiff_black_scholes.cpp
	$(CC) $(CFLAGS) -o black_scholes.exe BLACK_SCHOLES/autodiff_black_scholes.cpp
