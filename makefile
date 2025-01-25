# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -g -std=c++17 -fdiagnostics-color=auto

# Target executable
TARGET = Qmesh

# Source files
SRCS = Qmesh.cpp QuadMesh.cpp Display.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default rule to build and run the executable
all: $(TARGET)

# Rule to link object files into the target executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

# Rule to compile .cpp files into .o files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to run the executable
run: $(TARGET)
	./$(TARGET)

# Clean rule to remove generated files
clean:
	rm -f $(TARGET) $(OBJS)

