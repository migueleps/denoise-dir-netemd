# gtrieScanner: quick discovery of network motifs
#
# Pedro Ribeiro & David Aparício - CRACS & INESC-TEC, DCC/FCUP
# ------------------------------------------------
# Makefile
# ------------------------------------------------

# Executable name
EXEC_NAME=GTScanner

# Compiler to use
CC=g++

# Default Options to use
# CFLAGS= -Wall -Wno-write-strings -O0 -g
CFLAGS= -Wno-write-strings -O3 -g -std=c++11
# CFLAGS= -Wno-write-strings -O3 -g
# CFLAGSFINAL= -O3
CLIBS= -pthread

#NAUTYSRC= nauty.c nautil.c naugraph.c

# Source files
SRC =                   \
	CmdLine.cpp	\
	GraphMatrix.cpp	\
	GraphUtils.cpp	\
	Esu.cpp		\
        Fase.cpp	\
	FaseGTrie.cpp	\
	LSLabeling.cpp	\
	Conditions.cpp	\
	Isomorphism.cpp	\
	Error.cpp	\
	GTrie.cpp	\
	GDA.cpp		\
	Timer.cpp	\
	Random.cpp	\
	GraphTree.cpp	\
	nauty/nauty.c	\
	nauty/nautil.c	\
	nauty/naugraph.c	\
	main.cpp

OBJ =  ${SRC:.cpp=.o}

#------------------------------------------------------------

all: ${EXEC_NAME}

${EXEC_NAME}: ${OBJ}
	${CC} ${CFLAGS} -o ${EXEC_NAME} ${OBJ} ${CLIBS} 

%.o: %.cpp
	${CC} ${CFLAGS} -c -o $@ $+ ${CLIBS} 

clean:
	rm ${EXEC_NAME} *.o *~ *# -rf


