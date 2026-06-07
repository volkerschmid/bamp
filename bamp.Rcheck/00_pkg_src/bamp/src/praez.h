
#include <stdlib.h>



//datiert Hyperparameter auf (fuer A-P-C)
double hyper(int rw, double* theta, double k_a, double k_b, int n);
double hyper2(double* z, double d_g, double d_h, int n);
double hyper_a(double hyper2, int rw, double* theta, double d_g, double d_h, int n);



//datiert delta auf
double delta_berechnen(double** z, double d_g, double d_h, int n_o_a, int n_o_p);

//datiert delta auf im SAPC-Modell
double delta_berechnen_S(double*** z, double d_g, double d_h, int n_o_a, int n_o_p, int n_o_r);

//datiert Hyperparameter fuer Raum auf
double tau_berechnen(double* alpha, double k_a, double k_b, int** k_alpha, int number_of_regions);
