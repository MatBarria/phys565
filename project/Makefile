# Locations
BIN := ./bin
SRC := ./src
LIB := ./lib

# C++ compiler flags.
CXX   := g++ -g -Wall 
#FLAGS := -std=c++17 $(shell root-config --cflags) $(shell correction config --cflags --ldflags) 
#-I/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/boost/1.80.0-4e27c6908892e3529bfbcf4445a96fe6/include

#LIBS  := $(shell root-config --libs) $(shell correction config --cflags --ldflags) 
#-L/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/boost/1.80.0-4e27c6908892e3529bfbcf4445a96fe6/lib
FLAGS := -std=c++17 $(shell root-config --cflags)
LIBS  := $(shell root-config --libs)

# Tuple sources and header
TUPLE_SRC := $(SRC)/CreateTuples.cc
TUPLE_HEADERS := $(LIB)/CreateTuples.h
#p$(LIB)/Run3Constants.h $(LIB)/Constants.h

SOURCES :=  $(TUPLE_SRC)


TUPLE_OBJECTS := CreateTuples.o

OBJECTS :=  $(TUPLE_OBJECTS)

# Executables
TUPLE_EXECUTABLE := $(BIN)/CreateTuples

# Default target: build both executables
.PHONY: all
all: $(TUPLE_EXECUTABLE)


# Compile all the objects
$(OBJECTS): %.o: $(SRC)/%.cc $(LIB)/%.h
#$(OBJECTS): %.o: $(SOURCES) $(LIB)/%.h
	@echo "Compiling $<..."
	$(CXX) $(FLAGS) -c $< -o $@

# Linking for Tuples
$(TUPLE_EXECUTABLE): $(TUPLE_OBJECTS)
	@echo "Linking $(TUPLE_EXECUTABLE)..."
	$(CXX) $(FLAGS) -o $@ $^ $(LIBS)
	@echo "done"
	
# Specifying the object files as intermediates deletes them automatically after the build process.
.INTERMEDIATE:$(TUPLE_OBJECTS)

# Clean target
.PHONY: clean
clean:
	rm -f $(TUPLE_OBJECTS) $(TUPLE_EXECUTABLE)
