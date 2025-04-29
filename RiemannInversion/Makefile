CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall
TARGET = inverter
SRC = inverter.cpp
CSV = inverted_*.csv

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET) $(CSV)

run: 
	./inverter 
