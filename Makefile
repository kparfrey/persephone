CC = mpic++

INCLUDES =
CPPFLAGS = -std=c++17 -Wall -Werror
LDFLAGS  = -lmpi -lm
#HDRS     = process.hpp element_block.hpp params.hpp common.hpp active_params.hpp kernels.hpp
#OBJ      = main.o kernels.o

#%.o: %.cpp $(HDRS)
#	$(CC) -c $(CPPFLAGS) $< -o $@


persephone: main.o kernels.o process.o params_cartesian.o
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -o $@


main.o: main.cpp active_params.hpp process.hpp params.hpp common.hpp
	$(CC) -c $(CPPFLAGS) $< -o $@

kernels.o: kernels.cpp kernels.hpp common.hpp
	$(CC) -c $(CPPFLAGS) $< -o $@

process.o: process.cpp process.hpp common.hpp
	$(CC) -c $(CPPFLAGS) $< -o $@

params_cartesian.o: params_cartesian.cpp params_cartesian.hpp params.hpp common.hpp
	$(CC) -c $(CPPFLAGS) $< -o $@


clean:
	rm -f persephone
	rm -f *.o
