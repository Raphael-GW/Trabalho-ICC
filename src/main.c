#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <likwid.h>

#include "broyden.h"
#include "utils.h"

int main() {
    int n, max_iter; // define  o tamanho do sistema e o número máximo de iterações
    double x0, tol; // define o chute inicial e a tolerância para convergência

    printf ("Digite a dimensao do sistema (N): ");
    scanf ("%d", &n);
    printf ("Digite o chute inicial (x0): ");
    scanf ("%lf", &x0);
    printf ("Digite a tolerancia (TOL): ");
    scanf ("%lf", &tol);
    printf ("Digite o numero maximo de iteracoes: ");
    scanf ("%d", &max_iter);

    long double *x =  malloc(n * sizeof(long double)); // vetor com as variáveis x
    long double *f =  malloc(n * sizeof(long double)); // vetor com cada resultado de f(x)
    long double *s =  malloc(n * sizeof(long double)); // delta usado para atualizar x (i+1)
    
    /*
    long double **j =  malloc(n * sizeof(long double*)); // matriz jacobiana

    for (int i = 0; i < n; i++){
        j[i] = malloc(n * sizeof(long double));
    }
    */
    
    //Jacobiana otimizada
    long double *a =  malloc(n * sizeof(long double)); // subdiagonal
    long double *b =  malloc(n * sizeof(long double)); // diagonal principal
    long double *c =  malloc(n * sizeof(long double)); // superdiagonal
    

    /* chute inicial para todos */
    for (int i = 0; i < n; i++) x[i] = x0;
    
    LIKWID_MARKER_INIT; // Inicializa o sistema de marcação do LIKWID

    /*
    avaliaF(x, f, n); //monta o vetor F(X)

    rtime_t tempo_montaJ = 0.0, tempo_resolucao = 0.0;

    LIKWID_MARKER_START("Newton");
    rtime_t tempo_newton = timestamp();
    newton(j, f, s, x, tol, max_iter, n, &tempo_montaJ, &tempo_resolucao);
    tempo_newton = timestamp() - tempo_newton;
    LIKWID_MARKER_STOP("Newton");
    */

    
    // Roda Agora com as Otimizacoes
    rtime_t tempo_montaJ_ot = 0.0, tempo_resolucao_ot = 0.0;
    for (int i = 0; i < n; i++) x[i] = x0;
    
    avaliaF_Otim (x, f, n); //monta o vetor F(X)

    LIKWID_MARKER_START("Newton_Otim");
    rtime_t tempo_newton_ot = timestamp();
    newton_Otim (a, b, c, f, s, x, tol, max_iter, n, &tempo_montaJ_ot, &tempo_resolucao_ot);
    tempo_newton_ot = timestamp() - tempo_newton_ot;
    LIKWID_MARKER_STOP("Newton_Otim");
    
    /*
    printf ("##########\n");
    printf ("# Tempo Total: %f\n", tempo_newton);
    printf ("# Tempo Jacobiana: %f\n", tempo_montaJ);
    printf ("# Tempo SL: %f\n", tempo_resolucao);
    printf ("##########\n");
    */
    
    printf ("##########\n");
    printf ("# Tempo Otimizado Total: %f\n", tempo_newton_ot);
    printf ("# Tempo Otimizado Jacobiana: %f\n", tempo_montaJ_ot);
    printf ("# Tempo Otimizado SL: %f\n", tempo_resolucao_ot);
    printf ("##########\n");
    

    free(x); 
    free(f); 
    free(s);
    /*
    for (int i = 0; i < n; ++i){
        free(j[i]);
    }
    free(j);
    */

    free(a);
    free(b);
    free(c);
    
    LIKWID_MARKER_CLOSE; // Finaliza o sistema de marcação do LIKWID
    return 0;
}
