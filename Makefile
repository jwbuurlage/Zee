CCPP = clang++
CCPP_FLAGS = -std=c++11 -Wfatal-errors -Wall

OUTPUT_DIR = bin
INCLUDE_DIRS = -Iinclude

all: dirs test_spmv

dirs:
	mkdir -p ${OUTPUT_DIR}

test_spmv: test/test_spmv.cpp 
	${CCPP} ${CCPP_FLAGS} ${INCLUDE_DIRS} -o ${OUTPUT_DIR}/$@ $<

docs:
	doxygen docs/Doxyfile
