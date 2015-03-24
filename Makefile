CCPP = g++
CCPP_FLAGS = -std=c++11

OUTPUT_DIR = bin
INCLUDE_DIRS = include

all: test_spmv

test_spmv: test/test_spmv.cpp 
	${CCPP} ${CCPP_FLAGS} -I${INCLUDE_DIRS} -o ${OUTPUT_DIR}/$@ $<

docs:
	doxygen docs/Doxyfile
