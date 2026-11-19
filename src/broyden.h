#ifndef __BROYDEN_H__
#define __BROYDEN_H__

#include "utils.h"

void montaJacobiana(double *x, double **j, int n);

void avaliaF(double *x, double *f, int n);

void newton(double **j, double *f, double *s, double *x, double tol, int max_iter, int n, rtime_t *tempoJ, rtime_t *tempo_resolucao);

double norma(double *v, int n);


void montaJacobiana_Otim (double *x, double *a, double *b, double *c, int n);

void avaliaF_Otim (double *x, double *f, int n);

void newton_Otim (double *a, double *b, double *c, double *f, double *s, 
	double *x, double tol, int max_iter, int n, rtime_t *tempoJ_ot, rtime_t *tempo_resolucao_ot);

double norma_Otim (double *v, int n);

#endif
