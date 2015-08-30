CCPP = g++
CCPP_FLAGS = -std=c++14 -Wfatal-errors -Wall -g

OUTPUT_DIR = .
INCLUDE_DIRS = -I../include

LIB_DEPS = -lpthread

SOURCES = catch.cpp linear_algebra.cpp

all: dirs tests

dirs:
	mkdir -p ${OUTPUT_DIR}

tests: $(SOURCES)
	@echo 'CC $(SOURCES)'
	@${CCPP} ${CCPP_FLAGS} ${INCLUDE_DIRS} -o ${OUTPUT_DIR}/$@ $(SOURCES) ${LIB_DEPS}
