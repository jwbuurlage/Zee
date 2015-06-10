CCPP = clang++
CCPP_FLAGS = -std=c++14 -Wfatal-errors -Wall

OUTPUT_DIR = bin
INCLUDE_DIRS = -Iinclude

LIB_DEPS = -lpthread

.PHONY: lint test docs

all: dirs spmv cycpart part

dirs:
	mkdir -p ${OUTPUT_DIR}

spmv: examples/spmv.cpp 
	${CCPP} ${CCPP_FLAGS} ${INCLUDE_DIRS} -o ${OUTPUT_DIR}/$@ $< ${LIB_DEPS}

cycpart: examples/cyclic_partitioner.cpp
	${CCPP} ${CCPP_FLAGS} ${INCLUDE_DIRS} -o ${OUTPUT_DIR}/$@ $< ${LIB_DEPS}

part: examples/partitioner.cpp
	${CCPP} ${CCPP_FLAGS} ${INCLUDE_DIRS} -o ${OUTPUT_DIR}/$@ $< ${LIB_DEPS}

lint:
	./script/cpplint.py --filter=-whitespace,-build/c++11 --extensions=hpp include/*.hpp include/*/*.hpp

test:
	./script/test.py

docs:
	doxygen docs/Doxyfile
