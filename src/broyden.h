#ifndef __BROYDEN_H__
#define __BROYDEN_H__

#include "utils.h"

void montaJacobiana(long double *x, long double **j, int n);

void avaliaF(long double *x, long double *f, int n);

//newton otimizado
//void newton(long double *a, long double *b, long double *c, long double *f, long double *s, long double *x, long double tol, int max_iter, int n, rtime_t tempoJ, rtime_t tempo_resolucao);

void newton(long double **j, long double *f, long double *s, long double *x, long double tol, int max_iter, int n, rtime_t *tempoJ, rtime_t *tempo_resolucao);

long double norma(long double *v, int n);


void montaJacobiana_Otim (long double *x, long double *a, long double *b, long double *c, int n);

void avaliaF_Otim (long double *x, long double *f, int n);

void newton_Otim (long double *a, long double *b, long double *c, long double *f, long double *s, 
	long double *x, long double tol, int max_iter, int n, rtime_t *tempoJ_ot, rtime_t *tempo_resolucao_ot);

long double norma_Otim (long double *v, int n);

#endif
