CC = mpic++

INCLUDES =
CPPFLAGS = -std=c++17 -Wall -Werror
LDFLAGS  = -lmpi -lm
HDRS     = process.hpp element_block.hpp params.hpp common.hpp active_params.hpp kernels.hpp
OBJ      = main.o kernels.o

%.o: %.cpp $(HDRS)
	$(CC) -c $(CPPFLAGS) $< -o $@


persephone: $(OBJ)
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -o $@

clean:
	rm -f persephone
	rm -f *.o
