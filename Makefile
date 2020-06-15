CXX = mpic++

TARGET_ARCH =
CPPFLAGS = -I../HighFive/include
CXXFLAGS = -std=c++17 -Wall -Werror -O3
LDFLAGS  = 
LDLIBS   = -lmpi -lm -lhdf5

SRC = main.cpp kernels.cpp process.cpp element_block.cpp params.cpp params_cartesian.cpp \
	  write_mesh.cpp metric.cpp matrix.cpp initial_state_cartesian.cpp write_screen.cpp \
	  write_file_utils.cpp write_data.cpp 
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
