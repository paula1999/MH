HOME = .
BIN = $(HOME)/bin
INCLUDE = $(HOME)/include
SRC = $(HOME)/src
OBJ = $(HOME)/obj
LIB = $(HOME)/lib
DATA = $(HOME)/data

all: $(BIN)/practicaPAR

$(BIN)/practicaPAR: $(OBJ)/random.o $(OBJ)/PAR.o $(OBJ)/main.o
	g++ $^ -o $@

$(OBJ)/main.o: $(SRC)/main.cpp
	g++ -c -O3 $^ -I$(INCLUDE) -o $@

$(OBJ)/random.o: $(SRC)/random.cpp
	g++ -c -O3 $^ -I$(INCLUDE) -o $@

$(OBJ)/PAR.o: $(SRC)/PAR.cpp
	g++ -c -O3 $^ -I$(INCLUDE) -o $@

clean:
	-rm $(OBJ)/*.o
	-rm $(BIN)/*
