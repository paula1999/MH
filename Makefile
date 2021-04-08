HOME = .
BIN = $(HOME)/bin
INCLUDE = $(HOME)/include
SRC = $(HOME)/src
OBJ = $(HOME)/obj
LIB = $(HOME)/lib
DATA = $(HOME)/data

all: $(BIN)/practica1

$(BIN)/practica1: $(OBJ)/random.o $(OBJ)/PAR.o $(OBJ)/main.o
	g++ $^ -o $@

$(OBJ)/main.o: $(SRC)/main.cpp
	g++ -c $^ -I$(INCLUDE) -o $@

$(OBJ)/random.o: $(SRC)/random.cpp
	g++ -c $^ -I$(INCLUDE) -o $@

$(OBJ)/PAR.o: $(SRC)/PAR.cpp
	g++ -c $^ -I$(INCLUDE) -o $@

clean:
	-rm $(OBJ)/*.o
	-rm $(BIN)/*
