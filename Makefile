# Makefile pentru proiect

# Variabile
CXX = g++
CXXFLAGS = -Iinclude $(shell root-config --glibs --cflags --libs)
SRC = src/*.cpp
TARGET = task

# Regula pentru compilare
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(SRC) $(CXXFLAGS) -o $(TARGET)

# Regula pentru rulare
run: $(TARGET)
	./$(TARGET) 10 Europium data/data.root mDelila_raw data/energy.txt 10.0 1500.0 20 0 10000000 output/ Europiu-152

# Curățare
clean:
	rm -f $(TARGET)

