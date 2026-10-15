# Resolução do Sistema de Broyden pelo Método de Newton

Este guia descreve o procedimento técnico para resolver o **Sistema de Broyden** utilizando o **Método de Newton**, conforme os princípios de computação científica.

## 1. O Problema
O Sistema de Broyden é um exemplo clássico de sistema de equações não-lineares onde a matriz Jacobiana possui uma estrutura tridiagonal. Para um sistema de $n$ dimensões, a $i$-ésima função $f_i$ é geralmente definida como:

$$f_i(x_1, \dots, x_n) = (3 - 2x_i)x_i - x_{i-1} - 2x_{i+1} + 1 = 0$$

## 2. A Matriz Jacobiana ($J$)
O Método de Newton requer o cálculo das derivadas parciais para formar a matriz Jacobiana $J(X)$. Para o sistema acima, a matriz é esparsa (tridiagonal):

* **Diagonal Principal ($i=j$):** $rac{\partial f_i}{\partial x_i} = 3 - 4x_i$
* **Subdiagonal ($j=i-1$):** $rac{\partial f_i}{\partial x_{i-1}} = -1$
* **Superdiagonal ($j=i+1$):** $rac{\partial f_i}{\partial x_{i+1}} = -2$

## 3. Fluxo de Execução Passo a Passo

Para encontrar a solução aproximada partindo de um chute inicial $X^{(0)}$, repetimos as seguintes etapas na iteração $k$:

### Passo 1: Avaliação do Sistema ($F$)
Calcula-se o vetor de funções $F(X^{(k)})$ com os valores atuais das variáveis.

### Passo 2: Construção da Jacobiana ($J$)
Atualiza-se a matriz $J(X^{(k)})$ calculando os termos da diagonal principal e mantendo os termos constantes das diagonais adjacentes.

### Passo 3: Resolução do Sistema Linear
Resolve-se o sistema linear para encontrar o passo de correção $\Delta^{(k)}$:
$$J(X^{(k)}) \Delta^{(k)} = -F(X^{(k)})$$
*Dica: Como J é tridiagonal, utiliza-se o **Algoritmo de Thomas** para obter complexidade O(n).*

### Passo 4: Atualização do Vetor de Solução
Atualiza-se o palpite da raiz:
$$X^{(k+1)} = X^{(k)} + \Delta^{(k)}$$

### Passo 5: Critérios de Parada
O loop encerra se:
1.  A norma do resíduo for suficientemente pequena: $\|F(X^{(k)})\| < \epsilon_1$
2.  A correção aplicada for insignificante: $\|\Delta^{(k)}\| < \epsilon_2$

## 4. Por que este método?
* **Rapidez:** O método de Newton possui convergência quadrática nas proximidades da raiz.
* **Baixo Custo:** A estrutura tridiagonal do Sistema de Broyden permite que cada iteração de Newton seja resolvida em tempo linear, tornando-o extremamente eficiente para grandes valores de $n$.

---
*Documento gerado com base no material de Introdução à Computação Científica.*
