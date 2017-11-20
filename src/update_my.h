
// ziehen aus der Normalverteilung;

void update_my_1(double &my, double** ksi, double* theta, double* phi, double* psi, int vielfaches_der_breite, int number_of_agegroups, int number_of_periods, double delta);


// ziehen aus vorschlagsdichte;
void update_my_mh(double &my, double** ksi, double* theta, double* phi, double* psi, int vielfaches_der_breite, int number_of_agegroups, int number_of_periods, int lung_summe, int** n, int** y, int& ja_my);
