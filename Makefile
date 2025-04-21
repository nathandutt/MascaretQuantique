CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall
TARGET = NLS_sim
SRC = NLS_sim.cpp
CSV = evolution.csv

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET) $(CSV)

run: 
	make clean && make all && ./NLS_sim && python3 plotter.py
