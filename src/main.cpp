

class Grid {};

/** This 
 */
struct Interaction
{
    const double J;
    double operator(StateVector phi_i, StateVector phi_j);
private:
    double [] Dij;
};

struct OnSite
{
    const double[] h;
    StateVector operator(StateVector phi);
private:
    double [] data;
};

struct MultiSite
{
    const double[] k;
    StateVector operator(StateSpace s);
private:
    double data;
};



class Hamiltonian {
    typedef something StateVector;
    const uint Nalpha;
    const uint Nbeta;
    const uint Ngamma;
    
    Interaction[Nalpha] interactions;
    OnSite[Nbeta] onsite;
    MultiSite[Ngamma] multisite;
};

class Marqov {};

void metropolis();

void wolff();

class Observable {};



int main()
{
}
