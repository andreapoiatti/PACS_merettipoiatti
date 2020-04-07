EXE = program
OBJ = main.o
CXXFLAGS += -std=c++11 -Wall

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAS) $^ -o $@

main.o: carrier.h

all: $(EXE)

run: $(EXE)
	./$(EXE)

clean:
	$(RM) *.o

distclean:
	$(RM) $(EXE)
