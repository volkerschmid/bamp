
double taylor1(double x);
double taylor2(double x);
double taylor0(double x);

 

//Erzeugt standardnormalverteilten Zufallsvektor der Laenge noa
void gausssample(double* temp, int noa);

//Blockupdate fuer APC-Modell
void blockupdate(int swit, int age_block, double kappa, double delta, int noa, int nop, double** ksi, double &my, double* theta, double* phi, double* psi, double* ageQ, double* thetatemp, int vdb);

//Blockupdate fuer APC-Modell mit Kovariablen
void blockupdateplus(int swit, int age_block, double kappa, double delta, int noa, int nop, double** ksi, double &my, double* theta, double* phi, double* psi, double* ageQ, double* thetatemp, int vdb, double* data);
void blockupdateplus_linear(int swit, int age_block, double kappa, double delta, int noa, int nop, double** ksi, double &my, double* theta, double* phi, double* psi, double* ageQ, double* thetatemp, int vdb, double* data, double* data2, double alpha, double lambda);

//Blockupdate fuer SAPC-Modell
void blockupdate_S(int swit, int age_block, double kappa, double delta, int noa, int nop, double*** ksi,double &my,  double* theta, double* phi, double* psi, double* alpha, double* ageQ, double* thetatemp, int vdb, int nor);

//Bedinge Ergebnis aus Blockupdate mit Ax=b (L ist ageQ aus blockupdate)
void bedinge(int age_block, int noa, double* theta, double* L, double* thetatemp, int k, double* A, double* b);

//Blockupdate fuer APC-Modell mit unst. Effekten
void blockupdate_a2(int swit, int age_block, double kappa, double kappa2,double delta, int noa, int nop, double** ksi, double &my, double* theta,  double* theta2,double* phi, double* psi,  double* thetatemp, int vdb);

// APC-Modell ohne overdispersion
void blocken(int swit, int age_block, double kappa, double delta, int noa, int nop, double &my, double* theta, double* phi, double* psi, double* ageQ, double* thetatemp, int vdb, int** daten_n, int** daten_y, int &ja);
// APC-Modell ohne overdispersion
void blocken2(int swit, int age_block, double kappa, double delta, int noa, int nop, double &my, double* theta, double* theta2, double* phi, double* psi, double* ageQ, double* thetatemp, int vdb, int** daten_n, int** daten_y, int &ja);
