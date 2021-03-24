IDIR = ./include

CC = g++
CXXFLAGS = -W -Wall -std=c++17 -I$(IDIR)

SRC = ./src
ODIR = ./obj

LIBS = -lm
DEBUG = -g

_DEPS = $(patsubst $(IDIR)/%,%,$(DEPS))
DEPS = $(wildcard $(IDIR)/*.hpp)

_OBJ = $(patsubst $(SRC)/%.cpp,%.o,$(wildcard $(SRC)/*.cpp))
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SRC)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CXXFLAGS) $(DEBUG)

LZ77: $(OBJ)
	$(CC) -o $@ $^ $(CXXFLAGS) $(LIBS) $(DEBUG)

.PHONY: clean

clean:
	rm -f $(wildcard $(ODIR)/*.o)
