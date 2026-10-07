# BMS firmware - host build for unit tests and software-in-the-loop.
# (On the target MCU the same src/ files are compiled by the vendor toolchain.)
CC      := gcc
CFLAGS  := -std=c99 -O2 -Wall -Wextra -pedantic -Wshadow -Wdouble-promotion -Isrc
LDLIBS  := -lm

ifeq ($(OS),Windows_NT)
  EXE := .exe
  RM  := del /Q
else
  EXE :=
  RM  := rm -f
endif

SRC  := src/bms.c src/protection.c src/balancing.c src/soc_ekf.c
HDRS := $(wildcard src/*.h)

all: bms_sil$(EXE) test_bms$(EXE)

bms_sil$(EXE): $(SRC) sil/sil_main.c $(HDRS)
	$(CC) $(CFLAGS) -o $@ $(SRC) sil/sil_main.c $(LDLIBS)

test_bms$(EXE): $(SRC) tests/test_bms.c $(HDRS)
	$(CC) $(CFLAGS) -o $@ $(SRC) tests/test_bms.c $(LDLIBS)

test: test_bms$(EXE)
	./test_bms$(EXE)

clean:
	-$(RM) bms_sil$(EXE) test_bms$(EXE)

.PHONY: all test clean
