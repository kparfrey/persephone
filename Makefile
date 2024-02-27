### Laptop ###
CXX = mpic++
CPPFLAGS = -I../HighFive/include
CXXFLAGS = -std=c++17 -Werror -DHIGHFIVE_PARALLEL_HDF5=ON -O2 #-g -Og # -O3 -pg 
LDFLAGS  = #-pg
LDLIBS   = -lmpi -lm -lhdf5

### Workstation ###
# CXX = mpic++
# CPPFLAGS = -I../HighFive/include -I/usr/include/hdf5/openmpi
# CXXFLAGS = -std=c++17 -Wall -Werror -DHIGHFIVE_PARALLEL_HDF5=ON -g -Og # -O3 -pg 
# LDFLAGS  = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi #-pg
# LDLIBS   = -lmpi -lm -lhdf5

### Pleiades - Intel compiler & HPE MPI ###
# CXX = icpc
# CPPFLAGS = -I../HighFive/include 
# CXXFLAGS = -std=c++17 -Wall -Werror -DHIGHFIVE_PARALLEL_HDF5=ON -O3 -xAVX
# LDFLAGS  =
# LDLIBS   = -lstdc++fs -lmpi++ -lmpi -lm -lhdf5 

SRC = main.cpp kernels.cpp process.cpp element_block.cpp params.cpp  \
	  write_mesh.cpp geometry.cpp matrix.cpp write_screen.cpp  \
	  write_file_utils.cpp write_data.cpp lagrange_polynomials.cpp \
	  face_communicator.cpp edge.cpp transfinite_map.cpp \
	  params_cartesian.cpp initial_state_cartesian.cpp \
	  params_torus.cpp initial_state_torus.cpp cerfon_freidberg.cpp \
	  desc.cpp legendre_roots.cpp

OBJ := $(SRC:%.cpp=%.o)
DEP := $(SRC:%.cpp=%.d)


### Linking
persephone: $(OBJ)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

-include $(DEP)


### Compiling
DEPFLAGS = -MT $@ -MMD -MP -MF $*.d

%.o:%.cpp
	$(CXX) $(DEPFLAGS) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


clean:
	rm -f persephone
	rm -f *.o
	rm -f *.d


### Dependency management based on ideas from:
### http://scottmcpeak.com/autodepend/autodepend.html
### http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
