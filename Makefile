CC      ?= gcc
CFLAGS  ?= -O3 -march=native -Wall
LDLIBS  ?= -lm -lpthread -ldivsufsort
PREFIX  ?= /usr/local
DESTDIR ?=

dna_seed: dna_seed.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

.PHONY: test
test: dna_seed
	./test_roundtrip.sh

.PHONY: install
install: dna_seed
	install -Dm755 dna_seed $(DESTDIR)$(PREFIX)/bin/dna_seed

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/dna_seed

.PHONY: clean
clean:
	rm -f dna_seed *.dna *.restored
