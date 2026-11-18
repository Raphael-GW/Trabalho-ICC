#!/bin/bash
# ---------------------------------------------------------------------------
# Coleta os dados de desempenho do 'broyden' e gera CSVs e graficos.
#
# Uso: ./gendata.sh [PROG] [CPU] [X0] [TOL] [MAXIT]
#
# O broyden recebe a FLAG (-n = Newton normal, -o = Newton otimizado) na linha
# de comando e le via stdin, nesta ordem: N, x0, tol, max_iter.
#
# Parametros fixados pelo enunciado: TOL=0 e MAXIT=25 para todos os testes.
#
# Saidas (em resultados/<PROG>/):
#   Tempos.csv         tempos em ms medidos com timestamp() (execucao SEM likwid)
#   <GRUPO>.log        saida do likwid-perfctr de todos os N (um bloco por N)
#   *.csv / *.png      gerados pelo genplot.py (12 graficos)
# ---------------------------------------------------------------------------

CMD_DIR=$(dirname "$0")

PROG=${1:-broyden}
CPU=${2:-3}
X0=${3:--1.0}
TOL=${4:-0}
MAXIT=${5:-25}

DATA_DIR="resultados/${PROG}"
TEMPOS="${DATA_DIR}/Tempos.csv"

# Grupos do LIKWID (se o grupo MEM existir na maquina, troque L3 por MEM)
METRICAS="L3 FLOPS_DP FLOPS_AVX L2CACHE"

# Ate 2000 roda -n e -o. A partir de 4000 so a versao otimizada (-o) e executada.
TAMANHOS_AMBOS="32 64 128 256 512 1000 2000"
TAMANHOS_OTIM="4000 8000 9000 10000 20000"

#if ({1} = "-NDEBUG") {
make -DNDEBUG

mkdir -p "${DATA_DIR}"

likwid-setFrequencies -g performance -c "${CPU}"

# Executa o programa com a entrada padrao do enunciado: N, x0, tol, max_iter
entrada() {
    printf '%d\n%s\n%s\n%d\n' "$1" "${X0}" "${TOL}" "${MAXIT}"
}

# ---------------------------------------------------------------------------
# 1) Tempos (sem likwid, para nao medir a sobrecarga da instrumentacao).
#    Formato: N,versao,total_ms,jacobiana_ms,sl_ms
# ---------------------------------------------------------------------------
echo "N,versao,total_ms,jacobiana_ms,sl_ms" > "${TEMPOS}"

roda_tempo() {   # $1 = N, $2 = flag (-n|-o)
    local n=$1 modo=$2 versao saida
    [ "${modo}" = "-n" ] && versao="Normal" || versao="Otimizado"
    echo "--->>  TEMPO: ./${PROG} ${modo}  N=${n}" >&2

    saida=$(entrada "${n}" | ./${PROG} ${modo})
    # linhas: "# Tempo [Otimizado] Total|Jacobiana|SL: <valor>"
    local t_total t_jac t_sl
    t_total=$(echo "${saida}" | grep -E "^# Tempo.* Total:"     | awk '{print $NF}')
    t_jac=$(  echo "${saida}" | grep -E "^# Tempo.* Jacobiana:" | awk '{print $NF}')
    t_sl=$(   echo "${saida}" | grep -E "^# Tempo.* SL:"        | awk '{print $NF}')
    echo "${n},${versao},${t_total},${t_jac},${t_sl}" >> "${TEMPOS}"
}

for n in ${TAMANHOS_AMBOS}; do
    for modo in -n -o; do roda_tempo "${n}" "${modo}"; done
done
for n in ${TAMANHOS_OTIM}; do roda_tempo "${n}" -o; done

# ---------------------------------------------------------------------------
# 2) Contadores LIKWID. Os marcadores do broyden nao carregam o N no nome, por
#    isso gravamos "===> N: <n> <==" no log antes de cada bloco de dados.
# ---------------------------------------------------------------------------
roda_likwid() {  # $1 = grupo, $2 = N, $3 = flag, $4 = arquivo de log
    local m=$1 n=$2 modo=$3 log=$4
    local out="${DATA_DIR}/${m}_${n}_${modo#-}.txt"
    echo "--->>  ${m}: ./${PROG} ${modo}  N=${n}" >&2

    entrada "${n}" | \
        likwid-perfctr -O -C "${CPU}" -g "${m}" -o "${out}" -m ./${PROG} ${modo} > /dev/null

    cat "${out}" >> "${log}"
    rm -f "${out}"
}

for m in ${METRICAS}; do
    LIKWID_LOG="${DATA_DIR}/${m}.log"
    rm -f "${LIKWID_LOG}"

    for n in ${TAMANHOS_AMBOS}; do
        echo "===> N: ${n} <==" >> "${LIKWID_LOG}"
        for modo in -n -o; do roda_likwid "${m}" "${n}" "${modo}" "${LIKWID_LOG}"; done
    done

    for n in ${TAMANHOS_OTIM}; do
        echo "===> N: ${n} <==" >> "${LIKWID_LOG}"
        roda_likwid "${m}" "${n}" -o "${LIKWID_LOG}"
    done
done

likwid-setFrequencies -g powersave -c "${CPU}"

# ---------------------------------------------------------------------------
# 3) CSVs e graficos
# ---------------------------------------------------------------------------
python3 "${CMD_DIR}/genplot.py" "${DATA_DIR}"

echo ""
echo "Concluido. Graficos em ${DATA_DIR}/"
