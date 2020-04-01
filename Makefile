CC = mpic++

INCLUDES =
CPPFLAGS = -std=c++17 -Wall -Werror
LDFLAGS  = -lmpi -lm
HDRS     = process.hpp
OBJ      = main.o

%.o: %.cpp $(HDRS)
	$(CC) -c $(CPPFLAGS) $< -o $@


persephone: $(OBJ)
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -o $@

clean:
	rm -f persephone
	rm -f *.o
