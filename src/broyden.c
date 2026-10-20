#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <likwid.h>

#include "utils.h"

/* --------- F(x): calcula o vetor F no ponto x --------- */
void avaliaF(long double *x, long double *f, int n) {
    f[0] = -2.0*x[0]*x[0] + 3.0*x[0] - 2.0*x[1] + 1.0;

    for (int i = 1; i < n - 1; i++) {
        f[i] = -2.0*x[i]*x[i] + 3.0*x[i] - x[i-1] - 2.0*x[i+1] + 1.0;
    }
    
    f[n-1] = -2.0*x[n-1]*x[n-1] + 3.0*x[n-1] - x[n-2];
}

void avaliaF_Otim (long double *x, long double *f, int n) {
    f[0] = -2.0L*x[0]*x[0] + 3.0L*x[0] - 2.0L*x[1] + 1.0L;

    int limit = 1 + (n - 2) / 8 * 8; // Limite para o loop unroll, garantindo que não ultrapasse n-1
    int i;
    for (i = 1; i < limit; i += 8) {
        f[i] = -2.0L*x[i]*x[i] + 3.0L*x[i] - x[i-1] - 2.0L*x[i+1] + 1.0L;
        f[i+1] = -2.0L*x[i+1]*x[i+1] + 3.0L*x[i+1] - x[i] - 2.0L*x[i+2] + 1.0L;
        f[i+2] = -2.0L*x[i+2]*x[i+2] + 3.0L*x[i+2] - x[i+1] - 2.0L*x[i+3] + 1.0L;
        f[i+3] = -2.0L*x[i+3]*x[i+3] + 3.0L*x[i+3] - x[i+2] - 2.0L*x[i+4] + 1.0L;
        f[i+4] = -2.0L*x[i+4]*x[i+4] + 3.0L*x[i+4] - x[i+3] - 2.0L*x[i+5] + 1.0L;
        f[i+5] = -2.0L*x[i+5]*x[i+5] + 3.0L*x[i+5] - x[i+4] - 2.0L*x[i+6] + 1.0L;
        f[i+6] = -2.0L*x[i+6]*x[i+6] + 3.0L*x[i+6] - x[i+5] - 2.0L*x[i+7] + 1.0L;
        f[i+7] = -2.0L*x[i+7]*x[i+7] + 3.0L*x[i+7] - x[i+6] - 2.0L*x[i+8] + 1.0L;
    }
    for (; i < n - 1; i++) {
        f[i] = -2.0L*x[i]*x[i] + 3.0L*x[i] - x[i-1] - 2.0L*x[i+1] + 1.0L;
    }
    
    f[n-1] = -2.0L*x[n-1]*x[n-1] + 3.0L*x[n-1] - x[n-2];
}


/*a e c sao constantes (-1 e -2) e NAO sao alterados pelo solver,
  entao sao preenchidos uma unica vez. Apenas b depende de x.
 */
void inicializaJacobiana_Otim(long double *a, long double *c, int n) {
    a[0] = 0.0L;
    for (int i = 1; i < n; i++) a[i] = -1.0L;
    for (int i = 0; i < n - 1; i++) c[i] = -2.0L;
    c[n-1] = 0.0L;
}

/* --------- Monta as 3 diagonais da Jacobiana J(x) ---------
  a: subdiagonal  (indices 1..n-1 usados; a[0] ignorado)
  b: diagonal     (indices 0..n-1)
  c: superdiag.   (indices 0..n-2 usados; c[n-1] ignorado)
*/
void montaJacobiana_Otim(long double *x, long double *a,  long double *b, long double *c, int n) {
// Reconstrói a diagonal principal baseada no x atualizado
    for (int i = 0; i < n; i++) {
        b[i] = -4.0L * x[i] + 3.0L;
    }
    
}

void montaJacobiana(long double *x, long double **j, int n){
    for (int i = 0; i < n; ++i){
        for (int k = 0; k < n; ++k){
            j[i][k] = 0.0; // preenche a matriz com zeros
        }
    }
    
    j[0][0] = -4.0*x[0] + 3.0;
    j[0][1] = -2.0;
    for (int i = 1; i < n-1; ++i){ // monta as 3 diagonais
        j[i][i-1] = -1.0;
        j[i][i] = -4.0*x[i] + 3.0;
        j[i][i+1] = -2.0;
    }
    j[n-1][n-2] = -1.0;
    j[n-1][n-1] = -4.0*x[n-1] + 3.0;
}


// Eliminação de Gauss otimizada
void eliminacaoGauss_Otim (long double *a, long double *b, long double *c, long double *s, long double *f, int n){
    // Triangularização (Algoritmo de Thomas)
    f[0] = -f[0];
    
    int limit = 1 + ((n - 1) / 8) * 8; // Limite para o loop unroll, garantindo que não ultrapasse n-1
    int i;
    for (i = 1; i < limit; i += 8) {
        // Proteção essencial contra divisão por zero ou valores perigosamente próximos de zero
        
        long double m1 = a[i] / b[i-1];
        b[i] -= m1 * c[i-1];
        f[i] -= m1 * f[i-1];

        long double m2 = a[i+1] / b[i];
        b[i+1] -= m2 * c[i];
        f[i+1] -= m2 * f[i];

        long double m3 = a[i+2] / b[i+1];
        b[i+2] -= m3 * c[i+1];
        f[i+2] -= m3 * f[i+1];

        long double m4 = a[i+3] / b[i+2];
        b[i+3] -= m4 * c[i+2];
        f[i+3] -= m4 * f[i+2];

        long double m5 = a[i+4] / b[i+3];
        b[i+4] -= m5 * c[i+3];
        f[i+4] -= m5 * f[i+3];

        long double m6 = a[i+5] / b[i+4];
        b[i+5] -= m6 * c[i+4];
        f[i+5] -= m6 * f[i+4];

        long double m7 = a[i+6] / b[i+5];
        b[i+6] -= m7 * c[i+5];
        f[i+6] -= m7 * f[i+5];

        long double m8 = a[i+7] / b[i+6];
        b[i+7] -= m8 * c[i+6];
        f[i+7] -= m8 * f[i+6];
        
        
    }
    
    //Resto
    for (; i < n; ++i) {
        long double m1 = a[i] / b[i-1];
        b[i] -= m1 * c[i-1];
        f[i] -= m1 * f[i-1];
    }

    // Verifica todos os pivos (a condicao tambem captura NaN/Inf)
    for (int k = 0; k < n; k++) {
        if (!(fabsl(b[k]) > EPS_PIVO)) {
            fprintf(stderr, "Erro Numerico: pivo b[%d] nulo ou invalido.\n", k);
            exit(EXIT_FAILURE);
        }
    }

    // Retro-substituição
    s[n-1] = f[n-1] / b[n-1];
    int i = n - 2; // Declarado fora para manter o estado entre os loops

    // 1. Loop Desenrolado (Sempre que houver pelo menos 8 elementos restantes)
    for (; i >= 7; i -= 8) {
        s[i]   = (f[i]   - c[i]   * s[i+1]) / b[i];
        s[i-1] = (f[i-1] - c[i-1] * s[i])   / b[i-1];
        s[i-2] = (f[i-2] - c[i-2] * s[i-1]) / b[i-2];
        s[i-3] = (f[i-3] - c[i-3] * s[i-2]) / b[i-3];
        s[i-4] = (f[i-4] - c[i-4] * s[i-3]) / b[i-4];
        s[i-5] = (f[i-5] - c[i-5] * s[i-4]) / b[i-5];
        s[i-6] = (f[i-6] - c[i-6] * s[i-5]) / b[i-6];
        s[i-7] = (f[i-7] - c[i-7] * s[i-6]) / b[i-7];
    }

    // 2. Loop de Resto (Continua exatamente de onde o anterior parou até o 0)
    for (; i >= 0; i--) {
        s[i] = (f[i] - c[i] * s[i+1]) / b[i];
    }
}


void eliminacaoGauss (long double **j, long double *s, long double *f, int n){
    for (int i = 0; i < n; ++i) f[i] = -f[i];   // lado direito = -F

    for (int i = 0; i < n; ++i){
        for (int k = i+1; k < n; ++k){
            long double m = j[k][i] / j[i][i];
            j[k][i] = 0;
            for (int l = i+1; l < n; ++l)
                j[k][l] -= j[i][l]*m;
            f[k] -= f[i]*m;
        }
    }

    //retro-substituicao
    s[n-1] = f[n-1] / j[n-1][n-1];
    for (int i = n-2; i >= 0; --i){
        s[i] = f[i];
        for (int l = i+1; l < n; ++l)
            s[i] -= j[i][l]*s[l];
        s[i] /= j[i][i];
    }
}


/* --------- norma euclidiana --------- */
long double norma(long double *v, int n) {
    long double s = 0.0;
    for (int i = 0; i < n; i++) s += v[i] * v[i];
    return sqrtl(s);
}

void newton(long double **j, long double *f, long double *s, long double *x, long double tol, int max_iter, int n, rtime_t *tempoJ, rtime_t *tempo_resolucao) {
  rtime_t tempo;

  LIKWID_MARKER_START("Monta_Jacobiana");
  tempo = timestamp();
  montaJacobiana(x, j, n);
  *tempoJ += timestamp() - tempo;
  LIKWID_MARKER_STOP("Monta_Jacobiana");


  for (int k = 0; k < max_iter; k++) {
        
    long double nF = norma(f, n);

    if (nF < tol) {
      printf("%2d | %.6e |  (convergiu)\n", k, nF);
      printf("\nConvergencia em %d iteracoes.\n", k);
      break;
    }

	
    LIKWID_MARKER_START("Resolucao_SL");
    tempo = timestamp ();
    eliminacaoGauss(j, s, f, n); // resolve J(X) * s = -F(X)
	*tempo_resolucao += timestamp () - tempo;
    LIKWID_MARKER_STOP("Resolucao_SL");
	

    for (int i = 0; i < n; i++) x[i] += s[i]; // atualiza X(i+1)

    for (int i = 0; i < n; i++) {
      printf ("x%d = %.6Lf \n", i+1, x[i]);
    }
    printf("#\n");

    nF = norma(s, n);

    if (nF < tol) {
      printf("%2d | %.6e |  (convergiu)\n", k, nF);
      printf("\nConvergencia em %d iteracoes.\n", k);
      break;
    }

    avaliaF(x, f, n); //monta o vetor F(X)

    tempo = timestamp ();
        
    LIKWID_MARKER_START("Monta_Jacobiana");
    montaJacobiana(x, j, n);
    LIKWID_MARKER_STOP("Monta_Jacobiana");
        
	*tempoJ += timestamp() - tempo;
  }

  
}

/* --------- norma euclidiana otimizado --------- */
long double norma_Otim(long double *v, int n) {
    // Acumuladores independentes quebram a cadeia de dependencia da soma
    long double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
    int lim = n - n % 8;
    int i;
    for (i = 0; i < lim; i += 8){
        s0 += v[i]   * v[i]   + v[i+4] * v[i+4];
        s1 += v[i+1] * v[i+1] + v[i+5] * v[i+5];
        s2 += v[i+2] * v[i+2] + v[i+6] * v[i+6];
        s3 += v[i+3] * v[i+3] + v[i+7] * v[i+7];
    }
    for (; i < n; ++i){
        s0 += v[i] * v[i];
    }
 
    return sqrtl((s0 + s1) + (s2 + s3));
}


//Método de Newton Otimizado
void newton_Otim(long double *a, long double *b, long double *c, 
	long double *f, long double *s, long double *x, long double tol, 
	int max_iter, int n, rtime_t *tempoJ_ot, rtime_t *tempo_resolucao_ot) {
    
    rtime_t tempo;
    // a e c sao constantes: preenchidos uma unica vez
    inicializaJacobiana_Otim(a, c, n);


    LIKWID_MARKER_START("Monta_Jacobiana_Otim");
    tempo = timestamp();
    montaJacobiana_Otim (x, a, b, c, n); 
    *tempoJ_ot += timestamp() - tempo;
    LIKWID_MARKER_STOP("Monta_Jacobiana_Otim");


    for (int k = 0; k < max_iter; k++) {
        long double nF = norma_Otim (f, n);
        if (nF < tol) { 
            printf("%2d | %.6e |  (convergiu)\n", k, nF);
            printf("\nConvergencia em %d iteracoes.\n", k);
            break; 
        }


        LIKWID_MARKER_START("Resolucao_SL_Otim");
        tempo = timestamp ();
        eliminacaoGauss_Otim (a, b, c, s, f, n);
        *tempo_resolucao_ot += timestamp () - tempo;
        LIKWID_MARKER_STOP("Resolucao_SL_Otim");

        for (int i = 0; i < n; i++) x[i] += s[i];

        for (int i = 0; i < n; i++) {
            printf("x%d = %.6Lf \n", i+1, x[i]);
        }
        printf("#\n");

        nF = norma_Otim (s, n);
        if (nF < tol) { 
            printf("%2d | %.6e |  (convergiu)\n", k, nF);
            printf("\nConvergencia em %d iteracoes.\n", k);
            break; 
        }

        avaliaF_Otim (x, f, n);

        tempo = timestamp ();
        LIKWID_MARKER_START("Monta_Jacobiana_Otim");
        montaJacobiana_Otim (x, a, b, c, n); 
        LIKWID_MARKER_STOP("Monta_Jacobiana_Otim");
        *tempoJ_ot += timestamp() - tempo;
    }
}
