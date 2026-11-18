#!/bin/bash

make

./broyden -n < ../entradas/entrada_N32.txt > saida_normal
./broyden -o < ../entradas/entrada_N32.txt > saida_otim

