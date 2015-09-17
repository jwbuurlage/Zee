CCPP = g++
CCPP_FLAGS = -std=c++14 -Wfatal-errors -Wall -g

OUTPUT_DIR = bin
INCLUDE_DIRS = -Iinclude

LIB_DEPS = -lpthread

.PHONY: lint test docs

all: dirs part ir gmres lial

dirs:
	mkdir -p ${OUTPUT_DIR}

part: examples/partitioner.cpp
	${CCPP} ${CCPP_FLAGS} ${INCLUDE_DIRS} -o ${OUTPUT_DIR}/$@ $< ${LIB_DEPS}

ir: examples/ir.cpp
	${CCPP} ${CCPP_FLAGS} ${INCLUDE_DIRS} -o ${OUTPUT_DIR}/$@ $< ${LIB_DEPS}

gmres: examples/gmres.cpp
	${CCPP} ${CCPP_FLAGS} ${INCLUDE_DIRS} -o ${OUTPUT_DIR}/$@ $< ${LIB_DEPS}

lial: examples/linear_algebra.cpp
	${CCPP} ${CCPP_FLAGS} ${INCLUDE_DIRS} -o ${OUTPUT_DIR}/$@ $< ${LIB_DEPS}

dense: examples/dense.cpp
	${CCPP} ${CCPP_FLAGS} ${INCLUDE_DIRS} -o ${OUTPUT_DIR}/$@ $< ${LIB_DEPS}

lint:
	./script/cpplint.py --filter=-whitespace,-build/c++11 --extensions=hpp include/*.hpp include/*/*.hpp

test:
	./script/test.py

docs:
	doxygen docs/Doxyfile
	@cd docs/sphinx && make html
	firefox docs/sphinx/_build/html/index.html&
	@cd ../..
