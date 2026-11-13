#! /usr/bin/env python3
# coding=utf8

# Descomente as 2 linhas abaixo se quiser que Python procure por modulos
# em outros diretórios que não os do sistema e o diretório corrente.
## import sys, os
## sys.path.extend([ '~/lib/python', '.', '..' ])
# Ou em seu ambiente shell, crie e exporte a variavel de ambiente PYTHONPATH.
# Coloque a linha abaixo em seu arquivo '~/.bashrc' ou '~/.profile:
## export PYTHONPATH="${HOME}/lib/python:.:.."
#

import re, sys

#  "CACHE"    : "data cache misses",
#  "CACHE"    : "data cache miss ratio",
#  "CACHE"    : "data cache miss rate",
#  "L2CACHE"  : "L2 miss ratio",
#  "L2CACHE"  : "L2 miss rate",
#  "L3"       : "L3 bandwidth \[*MBytes/s\]*",
#  "MEM"      : "Memory bandwidth \[*MBytes/s\]*",
#  "MEM"      : "Memory read bandwidth",
#  "MEM"      : "Memory write bandwidth",
#  "ENERGY"   : "Power \[*W\]*",
#  "TIME"     : "Runtime \(RDTSC\) \[*s\]*",
#  "FLOPS_AVX": "AVX DP \[*MFLOP/s\]*",

campos = { "L2CACHE"  : "L2 miss ratio",
           "FLOPS_DP" : "DP .*MFLOP/s.*",
           "FLOPS_AVX": "Packed DP .*MFLOP/s.*",
           "L3"       : r"L3 bandwidth \[*MBytes/s]*"
         }

# 'STRUCT,Info' inicia Região
# 'TABLE,Region.*Metric,<metrica>' inicia metrica
# '<string>' inicia linha '<string>.*,<valor>,.*', onde
# <string> é obtido a partir de 'campos[]', tendo <metrica> como chave

# Linha gravada pelo gendata.sh ANTES de cada bloco de medidas, indicando
# o tamanho (N) do problema que está sendo medido a seguir. É necessária
# porque, diferente do matmult (onde o N vinha embutido no próprio nome
# do marcador, ex. "MatVet_64"), o broyden usa marcadores fixos
# ("Newton", "Newton_Otim") que não carregam o N.
RE_TAMANHO = re.compile(r'===>\s*N:\s*(\d+)\s*<==')

# Tamanho (N) "corrente": atualizado sempre que uma linha "===> N: ... <=="
# é encontrada no fluxo de entrada; usado como 'ordem' quando o nome do
# marcador não traz o N embutido.
n_atual = None


def lerDados():
    global n_atual

    linha = sys.stdin.readline()
    while linha and (re.match("STRUCT,Info", linha) is None):
        m = RE_TAMANHO.search(linha)
        if m:
            n_atual = int(m.group(1))
        linha = sys.stdin.readline()

    for linha in sys.stdin:
        m = RE_TAMANHO.search(linha)
        if m:
            n_atual = int(m.group(1))
        if re.match("TABLE,Region.*Metric,", linha) is not None:
            break

    if linha:
        linha = linha.split(',')
        metrica = linha[3]
        regiao = linha[1].split('_')

        # Se o último pedaço do nome do marcador for numérico, o N está
        # embutido nele (comportamento original, usado pelo matmult:
        # ex. "MatVet_64", "MatMat_Otim_128").
        if regiao[-1].isnumeric():
            ordem = int(regiao[-1])
            marker = re.sub("Region ", "", '_'.join(regiao[:-1]))
        else:
            # Caso contrário (ex. "Newton", "Newton_Otim", no broyden),
            # usa o N informado externamente pela linha "===> N: ... <=="
            if n_atual is None:
                sys.stderr.write(
                    "AVISO: tamanho (N) desconhecido para o marcador "
                    "'{0}' (nenhuma linha '===> N: ... <==' encontrada "
                    "antes dele)\n".format(linha[1]))
            ordem = n_atual
            marker = re.sub("Region ", "", '_'.join(regiao))

        for linha in sys.stdin:
            if re.match(campos[metrica], linha) is not None:
                break

        if linha:
            linha = linha.split(',')
            valor = float(linha[1])

            return [metrica, ordem, marker, valor]

    return ''

# fim lerDados()


# =====================================================
# Início Programa principal
# =====================================================
#
# Lê saída de likwid-perfctr -O e gera saída para gnuplot/matplotlib:
# um "gráfico" (tabela) para cada métrica, contendo os dados de todos
# os markers, indexados por N.
#

# saida = {<metrica> : {<ordem> : [ [<marker>, <valor>] ]} }
saida = {}

# item possui campos: metrica, ordem(N), marker, valor.
# 'ordem' é inteiro, 'valor' é double. 'marker' e 'metrica' são strings

item = lerDados()
while item:
    if item[0] in saida:
        if item[1] in saida[item[0]]:
            saida[item[0]][item[1]].append(item[2:])
        else:
            saida[item[0]][item[1]] = [item[2:]]
    else:
        saida[item[0]] = {item[1]: [item[2:]]}

    item = lerDados()

# gera tabela para uso em gnuplot ou pyplot/matplotlib
gnuplot = ''

for i in saida:
    ordens = []
    for j in saida[i]:
        ordens.append(j)
    ordens.sort()

    for j in ordens:
        gnuplot += '{0}'.format(j)
        for k in saida[i][j]:
            gnuplot += ',{0}'.format(k[1])

        gnuplot += '\n'
    # FIM for j in ordens

    # grava tabela de pontos
    print(gnuplot, end='')
# FIM for i in saida
