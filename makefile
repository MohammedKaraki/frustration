EXECUTABLE = collect
SOURCES = src/collect.cpp
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))

CXX = g++
CXX_GENERAL_FLAGS = -O3 -std=c++1z -Wall -Wpedantic -Wextra
CXX_LIBRARY_FLAGS = 

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) -o $(EXECUTABLE) $(OBJECTS) $(CXX_LIBRARY_FLAGS)

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXX_GENERAL_FLAGS)  -c $< -o $@

.dep: $(SOURCES)
	rm -f ./.dep
	$(CXX) $(CXX_GENERAL_FLAGS) -MM $^ > .dep


-include .dep

.PHONY: clean

clean:
	rm -rf *.o
	rm -rf .dep
