# test_counter build rule

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -I../include
TARGET = test_counter
SRC = test_counter.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)
