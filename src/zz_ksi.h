#include <stdlib.h>
#include <math.h>

#include "l.h"
#include "zufall.h"


//Metropolis Step mit LogNormal-Naeherung
void ZZ_aus_fc_von_ksi0(double my, double* theta, double* phi, double* psi, double** ksi, double delta, int number_of_agegroups, int number_of_periods, int vielfaches_der_breite, int** y, int** n, int** yes, int** no, int** schalter);
