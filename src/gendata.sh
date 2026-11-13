#!/bin/bash

CMD_DIR=`dirname $0`

# ---------------------------------------------------------------------------
# lê o 'broyden':
#   - recebe uma FLAG (-n ou -o) como argumento de linha de comando;
#   - lê, via stdin (scanf), nesta ordem: N, x0, tol, max_iter;
#   - usa marcadores LIKWID FIXOS ("Newton" e "Newton_Otim"), sem o N embutido.
# Por isso X0, TOL e MAX_ITER viraram parâmetros do script, e precisamos
# gravar o N no log manualmente (ver "===> N: ... <==" abaixo), já que o
# genplot.py não consegue mais extrair o N do nome do marcador.
# ---------------------------------------------------------------------------

PROG=${1:-broyden}
CPU=${2:-3}
X0=${3:--1.0}
TOL=${4:-0}
MAXIT=${5:-25}

DATA_DIR="resultados/${PROG}"
TEMPOS="${DATA_DIR}/Tempos.csv"

METRICA="L3 FLOPS_DP L2CACHE"

TAMANHOS_NOTIM="32 64 128 256 512 1000 2000 4000"
TAMANHOS_OTIM="8000 9000 10000 20000"
# Modos de execução do programa: Newton "puro" (-n) e Newton Otimizado (-o).
# Cada um ativa um marcador diferente dentro do main.c.
MODOS="-n -o"

mkdir -p ${DATA_DIR}

likwid-setFrequencies -g performance -c ${CPU}

# O Makefile_1 só tem o target 'all' (compila $(nomePrograma)); não existem
# mais os targets 'avx'/'debug' do matmult.


for m in ${METRICA}
do
    LIKWID_LOG="${DATA_DIR}/${m}.log"
    rm -f ${TEMPOS}
    rm -f ${LIKWID_LOG}

    for n in $TAMANHOS_NOTIM
    do
	# Grava no log qual N corresponde aos blocos de dados que virão a
	# seguir. O genplot.py usa essa linha para saber o eixo X, já que
	# o broyden não embute o N no nome do marcador.
	echo "===> N: ${n} <==" >> ${LIKWID_LOG}

	for modo in ${MODOS}
	do
	    LIKWID_OUT="${DATA_DIR}/${m}_${n}_${modo#-}.txt"

	    echo "--->>  $m: ./${PROG} ${modo}  N=${n}" >/dev/tty

	    # O broyden lê, via stdin (scanf), nesta ordem: N, x0, tol, max_iter
	    printf '%d\n%s\n%s\n%d\n' "${n}" "${X0}" "${TOL}" "${MAXIT}" | \
		likwid-perfctr -O -C ${CPU} -g ${m} -o ${LIKWID_OUT} -m ./${PROG} ${modo} >>${TEMPOS}

	    cat ${LIKWID_OUT} >> ${LIKWID_LOG}
	    rm -f ${LIKWID_OUT}
	done
    done
    for n in $TAMANHOS_OTIM
    do
	# Grava no log qual N corresponde aos blocos de dados que virão a
	# seguir. O genplot.py usa essa linha para saber o eixo X, já que
	# o broyden não embute o N no nome do marcador.
	echo "===> N: ${n} <==" >> ${LIKWID_LOG}

	modo="-n"
	
	LIKWID_OUT="${DATA_DIR}/${m}_${n}_${modo#-}.txt"

	echo "--->>  $m: ./${PROG} ${modo}  N=${n}" >/dev/tty

	# O broyden lê, via stdin (scanf), nesta ordem: N, x0, tol, max_iter
	printf '%d\n%s\n%s\n%d\n' "${n}" "${X0}" "${TOL}" "${MAXIT}" | \
	likwid-perfctr -O -C ${CPU} -g ${m} -o ${LIKWID_OUT} -m ./${PROG} ${modo} >>${TEMPOS}

	cat ${LIKWID_OUT} >> ${LIKWID_LOG}
	rm -f ${LIKWID_OUT}
	
    done

    # A partir dos logs, gera (para cada métrica) um CSV com uma coluna
    # por marcador (Newton, Newton_Otim), indexado por N.
    ${CMD_DIR}/genplot.py < ${LIKWID_LOG} > ${LIKWID_LOG%%.log}.csv
done

likwid-setFrequencies -g powersave -c ${CPU}

echo ""
echo ""
