#!/usr/bin/env python3
# coding=utf8
"""
Le os logs do likwid-perfctr (-O, CSV) e o Tempos.csv gerados pelo gendata.sh
e produz, em <dir>:
  - um CSV por teste (tempo, banda, l2cache, flops), com cabecalho;
  - 12 graficos de linhas PNG: 4 testes x 3 trechos (total, jacobiana, SL).

Uso: python3 genplot.py [diretorio]      (padrao: resultados/broyden)
"""
import csv
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = sys.argv[1] if len(sys.argv) > 1 else "resultados/broyden"

# regiao LIKWID -> (trecho, versao)
REGIOES = {
    "Newton":               ("total", "Normal"),
    "Newton_Otim":          ("total", "Otimizado"),
    "Monta_Jacobiana":      ("jac",   "Normal"),
    "Monta_Jacobiana_Otim": ("jac",   "Otimizado"),
    "Resolucao_SL":         ("sl",    "Normal"),
    "Resolucao_SL_Otim":    ("sl",    "Otimizado"),
}
TRECHOS = [("total", "Metodo de Newton (total)"),
           ("jac",   "Calculo das Jacobianas"),
           ("sl",    "Resolucao dos sistemas lineares")]
VERSOES = ["Normal", "Otimizado"]
COR = {"Normal": "tab:red", "Otimizado": "tab:blue"}

RE_N = re.compile(r"===>\s*N:\s*(\d+)\s*<==")
RE_REGIAO = re.compile(r"TABLE,Region (\S+),Group \d+ Metric,")


def le_log(arquivo, padroes):
    """Retorna {(trecho, versao): {N: valor}} para a metrica `padroes` (regex
    aplicadas no inicio da linha; a 1a que casar dentro da regiao e usada)."""
    dados = {}
    n = None
    regiao = None
    achou = False
    with open(arquivo) as f:
        for linha in f:
            m = RE_N.search(linha)
            if m:
                n, regiao = int(m.group(1)), None
                continue
            m = RE_REGIAO.match(linha)
            if m:
                regiao, achou = m.group(1), False
                continue
            if linha.startswith(("STRUCT", "TABLE")):
                regiao = None
                continue
            if regiao in REGIOES and not achou and n is not None:
                if any(re.match(p, linha) for p in padroes):
                    achou = True
                    try:
                        valor = float(linha.split(",")[1])
                    except (IndexError, ValueError):
                        continue
                    dados.setdefault(REGIOES[regiao], {})[n] = valor
    return dados


def le_tempos(arquivo):
    """Retorna {(trecho, versao): {N: ms}} a partir do Tempos.csv."""
    dados = {}
    colunas = {"total": "total_ms", "jac": "jacobiana_ms", "sl": "sl_ms"}
    with open(arquivo) as f:
        for r in csv.DictReader(f):
            for trecho, col in colunas.items():
                try:
                    dados.setdefault((trecho, r["versao"]), {})[int(r["N"])] = float(r[col])
                except (ValueError, TypeError):
                    pass
    return dados


def primeiro_existente(*nomes):
    for nome in nomes:
        caminho = os.path.join(DIR, nome)
        if os.path.exists(caminho):
            return caminho
    return None


def grava_csv(nome, series):
    """series: {rotulo: {N: valor}} -> CSV com cabecalho 'N,<rotulos>'."""
    ns = sorted({n for s in series.values() for n in s})
    rotulos = sorted(series)
    with open(os.path.join(DIR, nome), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["N"] + rotulos)
        for n in ns:
            w.writerow([n] + [series[r].get(n, "") for r in rotulos])


def grafico(nome, titulo, ylabel, series, logy=False):
    """series: {(versao, rotulo_extra, estilo): {N: valor}}"""
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for (versao, extra, estilo), pontos in series.items():
        if not pontos:
            continue
        ns = sorted(pontos)
        ax.plot(ns, [pontos[n] for n in ns], estilo, color=COR[versao],
                marker="o", markersize=4,
                label=versao + (" - " + extra if extra else ""))
    ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel("N (tamanho do SNLB)")
    ax.set_ylabel(ylabel)
    ax.set_title(titulo)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(DIR, nome), dpi=150)
    plt.close(fig)


def main():
    # cada teste: (id, titulo, ylabel, logy, [(arquivo(s), padroes, estilo, rotulo_extra)])
    banda_log = primeiro_existente("MEM.log", "L3.log")
    banda_pad = [r"Memory bandwidth \[MBytes/s\]"] if banda_log and banda_log.endswith("MEM.log") \
        else [r"L3 bandwidth \[MBytes/s\]"]
    banda_nome = "Banda de memoria (MEM)" if banda_log and banda_log.endswith("MEM.log") \
        else "Banda de memoria (L3)"

    testes = [
        ("tempo", "Tempo de execucao", "Tempo (ms)", True,
         [(os.path.join(DIR, "Tempos.csv"), None, "-", "")]),
        ("banda", banda_nome, "Banda [MBytes/s]", False,
         [(banda_log, banda_pad, "-", "")]),
        ("l2cache", "Cache miss ratio L2", "L2 miss ratio", False,
         [(primeiro_existente("L2CACHE.log"),
           [r"L2 miss ratio", r"data cache miss ratio"], "-", "")]),
        ("flops", "Operacoes aritmeticas", "MFLOP/s", False,
         [(primeiro_existente("FLOPS_DP.log"), [r"DP \[MFLOP/s\]"], "-", "DP"),
          (primeiro_existente("FLOPS_AVX.log"),
           [r"Packed DP .*MFLOP/s", r"AVX DP .*MFLOP/s"], "--", "AVX")]),
    ]

    for tid, titulo, ylabel, logy, fontes in testes:
        csv_series = {}
        lidas = []
        for arquivo, padroes, estilo, extra in fontes:
            if not arquivo or not os.path.exists(arquivo):
                sys.stderr.write("AVISO: arquivo ausente para '%s' (%s)\n" % (tid, extra or arquivo))
                continue
            dados = le_tempos(arquivo) if padroes is None else le_log(arquivo, padroes)
            lidas.append((extra, estilo, dados))
            for (trecho, versao), pontos in dados.items():
                rot = "%s_%s%s" % (trecho, versao, "_" + extra if extra else "")
                csv_series[rot] = pontos
        if not csv_series:
            continue
        grava_csv(tid + ".csv", csv_series)

        for trecho, nome_trecho in TRECHOS:
            series = {}
            for extra, estilo, dados in lidas:
                for versao in VERSOES:
                    series[(versao, extra, estilo)] = dados.get((trecho, versao), {})
            grafico("%s_%s.png" % (tid, trecho), "%s - %s" % (titulo, nome_trecho),
                    ylabel, series, logy)
            print("gerado %s_%s.png" % (tid, trecho))


if __name__ == "__main__":
    main()
