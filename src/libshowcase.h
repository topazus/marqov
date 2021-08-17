#include <string>
int pyIsing(std::string path, int dim, int len, double beta, double J);
int pyIsingGraph(std::string outpath, std::string graphpath, double beta, double J);
int pyHeisenberg(std::string path, int dim, int len, double beta, double J);
int pyPhi4(std::string path, int dim, int len, double beta, double lambda, double mass);
int pyAshkinTeller(std::string path, int dim, int len, double beta, double J, double K);
int pyBlumeCapel(std::string path, int dim, int len, double beta, double J, double D);
int pyBEG(std::string path, int dim, int len, double beta, double J, double D, double K);
int pyPotts(int q, std::string path, int dim, int len, double beta, double J);

