# Locations
BIN := ./bin
SRC := ./src
LIB := ./lib

CXX := g++ -g -Wall

FLAGS := -std=c++17 -I$(LIB) $(shell root-config --cflags)
LIBS  := $(shell root-config --libs)

SOURCES := $(SRC)/main.cc \
           $(SRC)/CreateTuples.cc \
           $(SRC)/kinematics.cc

OBJECTS := main.o \
           CreateTuples.o \
           kinematics.o

EXECUTABLE := $(BIN)/main

.PHONY: all
all: $(EXECUTABLE)

%.o: $(SRC)/%.cc
	@echo "Compiling $<..."
	$(CXX) $(FLAGS) -c $< -o $@

$(EXECUTABLE): $(OBJECTS)
	@mkdir -p $(BIN)
	@echo "Linking $(EXECUTABLE)..."
	$(CXX) $(FLAGS) -o $@ $^ $(LIBS)
	@echo "done"

.INTERMEDIATE: $(OBJECTS)

.PHONY: clean
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
