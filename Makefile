CXX = mpic++

TARGET_ARCH =
CPPFLAGS = 
CXXFLAGS = -std=c++17 -Wall -Werror
LDFLAGS  = 
LDLIBS   = -lmpi -lm

SRC = main.cpp kernels.cpp process.cpp element_block.cpp params.cpp params_cartesian.cpp 
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
