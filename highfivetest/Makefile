CXX = g++

TARGET_ARCH =
CPPFLAGS = -I../../HighFive/include #-I/usr/include/hdf5/openmpi
CXXFLAGS = -std=c++17 -Wall -Werror -DHIGHFIVE_PARALLEL_HDF5=ON -g -Og # -O3 -pg 
LDFLAGS  = #-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi #-pg
LDLIBS   = -lm -lhdf5 # -lmpi

loadtest: loadtest.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) loadtest.cpp $(LDLIBS) -o loadtest

