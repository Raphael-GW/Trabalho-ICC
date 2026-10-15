
LIKWID_PREFIX ?= /usr/local

LIKWID_INC = -I$(LIKWID_PREFIX)/include
LIKWID_LIB = -L$(LIKWID_PREFIX)/lib

# Flags de Compilação
# -g para o Valgrind e -DLIKWID_PERFMON para ativar as marcações do likwid.h
CFLAGS = -O3 -march=native -mavx -fopt-info-vec -g -DLIKWID_PERFMON $(LIKWID_INC)

# Flags de Linkagem (Bibliotecas)
LDFLAGS = $(LIKWID_LIB) -lm -llikwid -llikwid-hwloc -llikwid-lua \
          -Wl,-rpath,$(LIKWID_PREFIX)/lib \
          -Wl,-rpath-link,$(LIKWID_PREFIX)/lib

nomePrograma = broyden
FLAGS_VALGRIND = --leak-check=full --show-leak-kinds=all --track-origins=yes

all: $(nomePrograma)

$(nomePrograma): main.o broyden.o utils.o
	gcc -o $(nomePrograma) main.o broyden.o utils.o $(LDFLAGS)

main.o: main.c
	gcc -c main.c $(CFLAGS)

broyden.o: broyden.h broyden.c
	gcc -c broyden.c $(CFLAGS)

utils.o: utils.h utils.c
	gcc -c utils.c $(CFLAGS)

valgrind: $(nomePrograma)
	valgrind $(FLAGS_VALGRIND) ./$(nomePrograma)

clean:
	rm -f *.o *.gch $(nomePrograma)