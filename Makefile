CC = mpic++

INCLUDES =
CPPFLAGS = -std=c++17 -Wall -Werror
LDFLAGS  = -lmpi -lm


persephone: main.o kernels.o process.o params_cartesian.o element_block.o
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -o $@


main.o: main.cpp active_params.hpp params_cartesian.hpp params.hpp process.hpp common.hpp
	$(CC) -c $(CPPFLAGS) $< -o $@

kernels.o: kernels.cpp kernels.hpp common.hpp
	$(CC) -c $(CPPFLAGS) $< -o $@

process.o: process.cpp process.hpp params.hpp common.hpp
	$(CC) -c $(CPPFLAGS) $< -o $@

params_cartesian.o: params_cartesian.cpp params_cartesian.hpp params.hpp common.hpp
	$(CC) -c $(CPPFLAGS) $< -o $@

element_block.o: element_block.cpp element_block.hpp common.hpp
	$(CC) -c $(CPPFLAGS) $< -o $@


clean:
	rm -f persephone
	rm -f *.o
