CC=g++-13
#CFLAGS= -Ofast  -m64  -ftree-vectorize -mmmx -msse -msse2
G++FLAG= -std=c++17 -Ofast  -m64  -ftree-vectorize 
#-g
#-pthread -fopenmp  
# SRC_DIR := ./
# OBJ_DIR := ./
# SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
# OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
LDFLAGS := 
#-pthread

LINK_OPTS2 =  -lz 
#-fopenmp -lpthread
I_PATH = -I../include
#L_PATH = ../lib/macos-arm

all : make_directories adze 

adze : main.o  gzstream.o 
	$(CC) $(G++FLAG) -o adze main.o gzstream.o  $(LINK_OPTS2)

main.o : main.cpp
	$(CC) $(G++FLAG) -c main.cpp $(I_PATH)

gzstream.o : gzstream.cpp gzstream.h
	$(CC) $(G++FLAG) -c gzstream.cpp $(I_PATH)


.PHONY: make_directories
make_directories:
	mkdir -p ../bin/

clean:
	rm -rf *.o adze
