# -----------------------------------------------------------------------------
# Configuration Variables
# -----------------------------------------------------------------------------
CXX = g++                                       # C++ Compiler (assuming g++)
CXXFLAGS = -std=c++17 -Wall -Wextra -g          # Compilation flags: C++17, active Warnings, -g for debugging
TARGET = wave_solver                            # Name of the final executable file
SOURCES = main.cpp Mesh.cpp                     # List of all source files (.cpp)

# -----------------------------------------------------------------------------
# Rules (Targets)
# -----------------------------------------------------------------------------

# Main Rule: compiles the project
all: $(TARGET)

# Rule for building the executable (linking)
$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

# Rule for execution
run: $(TARGET)
	./$(TARGET)

# Rule for cleaning up (removes the executable)
clean:
	rm -f $(TARGET)