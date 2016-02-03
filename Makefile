CCPP = g++
CCPP_FLAGS = -std=c++14 -Wfatal-errors -Wall -g

OUTPUT_DIR = bin
INCLUDE_DIRS = -Iinclude

LIB_DEPS = -lpthread

.PHONY: lint test docs

all: dirs \
	$(OUTPUT_DIR)/partitioner \
	$(OUTPUT_DIR)/pulp_partitioner \
	$(OUTPUT_DIR)/pulp_experiments \
	$(OUTPUT_DIR)/vector_partitioner \
	$(OUTPUT_DIR)/ir \
	$(OUTPUT_DIR)/gmres \
	$(OUTPUT_DIR)/gmres_experiments \
	$(OUTPUT_DIR)/lial

.phony: docs

dirs:
	mkdir -p ${OUTPUT_DIR}

$(OUTPUT_DIR)/%: examples/%.cpp
	${CCPP} ${CCPP_FLAGS} ${INCLUDE_DIRS} -o $@ $< ${LIB_DEPS}

lint:
	./script/cpplint.py --filter=-whitespace,-build/c++11 --extensions=hpp include/*.hpp include/*/*.hpp

test:
	./script/test.py

docs:
	doxygen docs/Doxyfile
	@cd docs/sphinx && make html
	firefox docs/sphinx/_build/html/index.html&
	@cd ../..
