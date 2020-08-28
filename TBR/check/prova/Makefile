CXX = g++
CXXFLAGS += -Wall -std=c++11 \
   -Werror \
   -Wextra \
   -Wconversion \
   -Wno-deprecated \
   -Winit-self \
   -Wsign-conversion \
   -Wredundant-decls \
   -Wvla -Wshadow -Wctor-dtor-privacy -Wnon-virtual-dtor -Woverloaded-virtual \
   -Winit-self \
   -Wpointer-arith \
   -Wcast-qual \
   -Wcast-align \
   -Wdouble-promotion \
   -Wold-style-cast -Wno-error=old-style-cast \
   -Wsign-promo \
   -Wswitch-enum \
   -Wswitch-default \
   -Wundef

EXE = main

all: $(EXE)

.PHONY: all clean distclean

$(EXE): main.o
	$(CXX) $(CXXFLAGS) $^ -o $@

main.o: dummy.h


clean:
	$(RM) *.o

distclean: clean
	$(RM) $(EXE)
