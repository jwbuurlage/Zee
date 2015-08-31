CCPP = g++
CCPP_FLAGS = -std=c++14 -Wfatal-errors -Wall -g

OUTPUT_DIR = .
INCLUDE_DIRS = -I../include -I../unpain/include

LIB_DEPS = -lpthread

all: gmres

gmres: gmres_host.cpp streaming_matrix.hpp
	${CCPP} ${CCPP_FLAGS} ${INCLUDE_DIRS} -o ${OUTPUT_DIR}/$@ $< ${LIB_DEPS}
