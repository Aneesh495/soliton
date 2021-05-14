CC ?= clang
CFLAGS ?= -O3 -Wall -Wextra -std=c11 -Iengine/include
LDFLAGS ?= -lm

SRC_ENGINE = \
	engine/src/alloc.c \
	engine/src/sym.c \
	engine/src/core_verbs.c \
	engine/src/struct_verbs.c \
	engine/src/adverbs.c \
	engine/src/table.c \
	engine/src/temporal.c \
	engine/src/parser.c \
	engine/src/eval.c \
	engine/src/qsql.c \
	engine/src/ipc.c

OBJ_ENGINE = $(SRC_ENGINE:.c=.o)

all: bin/soliton bin/test_engine symlink

bin/soliton: $(OBJ_ENGINE) engine/src/main.o | bin
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

bin/test_engine: $(OBJ_ENGINE) engine/tests/test_engine.o | bin
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

symlink: bin/soliton
	@ln -sf soliton bin/q

bin:
	@mkdir -p bin

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: bin/test_engine
	@./bin/test_engine
	@echo "Running Q verification harness..."
	@python3 scripts/run_tests.py

bench: bin/soliton
	@python3 tools/benchmark.py

demo: bin/soliton
	@bash scripts/run_demo.sh

verify: test
	@python3 tools/verify_clean.py

clean:
	rm -rf bin $(OBJ_ENGINE) engine/src/main.o engine/tests/test_engine.o

.PHONY: all test bench demo verify clean symlink
