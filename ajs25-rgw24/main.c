#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <likwid.h>
#include <string.h>

#include "broyden.h"
#include "utils.h"

int main(int argc, char* argv[]) {
    int n, max_iter; // define  o tamanho do sistema e o número máximo de iterações
    double x0, tol; // define o chute inicial e a tolerância para convergência

	if (argc > 1 && strcmp(argv[1], "-o") != 0 && strcmp(argv[1], "-n") != 0) {
		fprintf(stderr, "Flag errada \n");
		return 1 ;
	}

    printf ("Digite a dimensao do sistema (N): ");
    scanf ("%d", &n);
    printf ("Digite o chute inicial (x0): ");
    scanf ("%lf", &x0);
    printf ("Digite a tolerancia (TOL): ");
    scanf ("%lf", &tol);
    printf ("Digite o numero maximo de iteracoes: ");
    scanf ("%d", &max_iter);

    double *x =  malloc(n * sizeof(double)); // vetor com as variáveis x
    double *f =  malloc(n * sizeof(double)); // vetor com cada resultado de f(x)
    double *s =  malloc(n * sizeof(double)); // delta usado para atualizar x (i+1)
    

    /* chute inicial para todos */
	for (int i = 0; i < n; i++) x[i] = x0;
		
	LIKWID_MARKER_INIT; // Inicializa o sistema de marcação do LIKWID
	if (argc > 1){
		if ( strcmp(argv[1], "-n") == 0) {

			double **j =  malloc(n * sizeof(double*)); // matriz jacobiana

			for (int i = 0; i < n; i++){
				j[i] = malloc(n * sizeof(double));
			}
			
			avaliaF(x, f, n); //monta o vetor F(X)

			rtime_t tempo_montaJ = 0.0, tempo_resolucao = 0.0;

			LIKWID_MARKER_START("Newton");
			rtime_t tempo_newton = timestamp();
			newton(j, f, s, x, tol, max_iter, n, &tempo_montaJ, &tempo_resolucao);
			tempo_newton = timestamp() - tempo_newton;
			LIKWID_MARKER_STOP("Newton");

			printf ("\n##########\n");
			printf ("# Tempo Total: %f\n", tempo_newton);
			printf ("# Tempo Jacobiana: %f\n", tempo_montaJ);
			printf ("# Tempo SL: %f\n", tempo_resolucao);
			printf ("##########\n");

			for (int i = 0; i < n; ++i){
        		free(j[i]);
    		}
    		free(j);
		
		}
		else {

			//Jacobiana otimizada
			double *a =  malloc(n * sizeof(double)); // subdiagonal
			double *b =  malloc(n * sizeof(double)); // diagonal principal
			double *c =  malloc(n * sizeof(double)); // superdiagonal

			// Roda Agora com as Otimizacoes
			rtime_t tempo_montaJ_ot = 0.0, tempo_resolucao_ot = 0.0;
			
			avaliaF_Otim (x, f, n); //monta o vetor F(X)

			LIKWID_MARKER_START("Newton_Otim");
			rtime_t tempo_newton_ot = timestamp();
			newton_Otim (a, b, c, f, s, x, tol, max_iter, n, &tempo_montaJ_ot, &tempo_resolucao_ot);
			tempo_newton_ot = timestamp() - tempo_newton_ot;
			LIKWID_MARKER_STOP("Newton_Otim");


			printf ("##########\n");
			printf ("# Tempo Otimizado Total: %f\n", tempo_newton_ot);
			printf ("# Tempo Otimizado Jacobiana: %f\n", tempo_montaJ_ot);
			printf ("# Tempo Otimizado SL: %f\n", tempo_resolucao_ot);
			printf ("##########\n");

			free(a);
			free(b);
			free(c);
			
		}
	}
	// executa o otimizado se nao tem argumentos
	else{
		// Roda Agora com as Otimizacoes
	        	double *a =  malloc(n * sizeof(double)); // subdiagonal
			double *b =  malloc(n * sizeof(double)); // diagonal principal
			double *c =  malloc(n * sizeof(double)); // superdiagonal

			rtime_t tempo_montaJ_ot = 0.0, tempo_resolucao_ot = 0.0;
			
			avaliaF_Otim (x, f, n); //monta o vetor F(X)

			LIKWID_MARKER_START("Newton_Otim");
			rtime_t tempo_newton_ot = timestamp();
			newton_Otim (a, b, c, f, s, x, tol, max_iter, n, &tempo_montaJ_ot, &tempo_resolucao_ot);
			tempo_newton_ot = timestamp() - tempo_newton_ot;
			LIKWID_MARKER_STOP("Newton_Otim");


			printf ("##########\n");
			printf ("# Tempo Otimizado Total: %f\n", tempo_newton_ot);
			printf ("# Tempo Otimizado Jacobiana: %f\n", tempo_montaJ_ot);
			printf ("# Tempo Otimizado SL: %f\n", tempo_resolucao_ot);
			printf ("##########\n");

			free(a);
			free(b);
			free(c);
	}

    free(x); 
    free(f); 
    free(s);
    
    LIKWID_MARKER_CLOSE; // Finaliza o sistema de marcação do LIKWID
    return 0;
}
