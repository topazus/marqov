#ifndef SSH_H
#define SSH_H
#include <array>
#include <tuple>
#include <string>
#include <complex>
#include <functional>
#include "hamiltonianparts.h"
#include "metropolis.h"

namespace ssh
{


	/**
	The fermi function
	@return An occupation propability
	*/
	template <typename T>
	inline T fermi(const T& e) throw()
	{
	    T xhalf = -e/T(2.0);
	    T retval;
	    if(xhalf > std::log(0.01*std::numeric_limits<T>::max()))
	      retval = 1.0;
	    else
	      if(xhalf < std::log(10*std::numeric_limits<T>::min()))
		retval = 0.0;//FIXME: not yet decided whether setting this to zero is better than setting it to epsilon
		else
		  retval = 0.5*std::exp(xhalf)/std::cosh(xhalf);
	    return retval;
	}

	// The Fermi-Dirac function
	double fermi_dirac(double x, double beta=1, double mu=0)
	{
		return 1.0 / (std::exp(beta*(x-mu)) + 1);
	}
		
	
	// Green's function for a regular hypercubic lattice in 1D (i.e. a chain)
	double g1D(std::vector<double> c1, std::vector<double> c2, int sign1, int sign2, int L, double beta, double dtau)
	{
        static double* data = NULL;
        int ntau = std::round(beta/dtau);
        if(data == NULL)
        {
            data = new double[ntau*L];
            for(int j = 0; j < L; ++j)
            {
                double k = (j*2)*M_PI/double(L);
                // 1D disperson relation
            
                const double mu = 0;
                double eps = -2.0*std::cos(k)-mu;
                for(int dt = 0; dt < ntau; ++dt)
                {
                    data[dt * L + j] = 0.5*std::exp((beta/2 - dt*dtau)*eps)/std::cosh(beta*eps/2);
                }
            }
        }
		std::complex<double> retval = 0;
		std::complex<double> jj(0,1);
	
		double step = 2*M_PI/double(L);
	
		// extract spatial and temporal coordinates
		// improve me: if c1=c2, r1-r2=+1/-1

		const double r1 = fmod(c1[0] + sign1*0.5, L); // account for p.b.c
		const double r2 = fmod(c2[0] + sign2*0.5, L);

//		cout << r1 << "  " << r2 << endl;

		double dist = r1-r2;
		double sign = 1;
		if (dist < 0) sign = -1;
		if (fabs(dist) > 0.5*L) 
		{
//			cout << "distance too long" << endl;
			dist = L-fabs(dist);
		}
		dist = dist*sign;

		const double t1 = c1[1];
        int t1i = round(t1);
		const double t2 = c2[1];
        int t2i = round(t2);
        int dti = t1i - t2i;

		// account for anti-symmetric behaviour in delta tau
		double signum = 1;
		if (dti < 0) 
		{
			signum = -1;
			dti = ntau-dti;
		}

		// compute Greens function in Fourier space
		std::complex<double> dexpk = std::exp(-jj*dist);
        std::complex<double> expk = 1.0;
		for (int j = 0; j < L; ++j)
		{
//             double k = (j*2)*M_PI/double(L);
// 			// 1D disperson relation
//             
// 	            const double mu = 0;
// 	            const double eps = -2.0*std::cos(k)-mu;
// 	            const double fac = 0.5*std::exp((beta/2 - dti*dtau)*eps)/std::cosh(beta*eps/2);
//                 std::cout<<fac<<" "<<data[dti*L + j]<<std::endl;
//	            const double occ = ssh::fermi(beta*eps);

			// direct Fourier summation (slow!)
			// todo: use FFT
//			retval += std::exp(-jj*k*(r1-r2))*fac;
//			retval += std::exp(-jj*k*dist)*fac;
            retval += expk.real()*data[dti*L + j];
            expk *= dexpk;//loop-carried dependency. breaks vectorization.
		}
		auto retv = signum*retval.real()/double(L);

		if (retv != retv) cout << "retv: " << retv << endl;

		return retv;
	}


/*

	 j1/i2      b2     j2/i3      b3      j3/i4      b4
  ------x-----------------x------------------x-----------------x---------- // original model: bonds

  ---------------o-----------------o------------------o----------------    // our model: bonds represented as sites


   e.g.  r(b2(j2)) = r(b3(i3))

*/


	std::complex<double> fourpointcorr(std::vector<std::vector<double>> s, std::vector<int> e, int L, double beta, double dtau)
	{
		std::complex<double> retval = 0;

		// the square of the number operator of fermions is the number operator
		if (s[0]==s[1] && s[1]==s[2] && s[2]==s[3])
		{
			if (e[0]==e[1] && e[1]==e[2] && e[2]==e[3])
			{
				retval = ssh::g1D(s[0], s[1], e[0], e[1], L, beta, dtau);
			}
		}

		// Wick decomposition
		else 
		{
			retval = ssh::g1D(s[0],s[1],e[0],e[1], L, beta,dtau)*ssh::g1D(s[2],s[3],e[2],e[3], L, beta,dtau) 
				  - ssh::g1D(s[0],s[3],e[0],e[3], L, beta,dtau)*ssh::g1D(s[2],s[1],e[2],e[1], L, beta,dtau);
		}

		return retval;
	}





	template <class Lattice>
	double suscept(Lattice& grid, int idx1, int idx2, double beta, double dtau)
	{
		auto w1 = grid.getcrds(idx1);
		auto w2 = grid.getcrds(idx2);

//		cout << w1[0] << "  " << w2[0] << endl;

		const int L = grid.len;
	
		// signs
		const int i = -1;
		const int j = +1;

		std::complex<double> retval = 0;

		/*

		< K(b1,t1) K(b2,t2) > = 	  
							  < c†(i1) c(j1) c†(i2) c(j2) > 
							+ < c†(i1) c(j1) c†(j2) c(i2) > 
							+ < c†(j1) c(i1) c†(i2) c(j2) > 
							+ < c†(j1) c(i1) c†(j2) c(i2) > // temporal indices supressed; b1 = (i1,j1), b2 = (i2,j2)
						  =
						  	  <i1 j1> <i2 j2> - <i1 i2> <i2 i1>
						  	+ <i1 j1> <j2 i2> - <i1 i2> <j2 j1>
						  	+ <j1 i1> <i2 j2> - <j1 j2> <i2 i1>
						  	+ <j1 i1> <j2 i2> - <j1 j2> <j2 j1> // operators supressed
		*/

		std::vector<std::vector<double>> sites = {w1,w1,w2,w2};

		retval += ssh::fourpointcorr(sites, {i,j,i,j}, L, beta, dtau);
		retval += ssh::fourpointcorr(sites, {i,j,j,i}, L, beta, dtau);
		retval += ssh::fourpointcorr(sites, {j,i,i,j}, L, beta, dtau);
		retval += ssh::fourpointcorr(sites, {j,i,j,i}, L, beta, dtau);

// 		debug things
//		------------
//		auto c1 = ssh::fourpointcorr(sites, {i,j,i,j}, L, beta, dtau);
//		auto c2 = ssh::fourpointcorr(sites, {i,j,j,i}, L, beta, dtau);
//		auto c3 = ssh::fourpointcorr(sites, {j,i,i,j}, L, beta, dtau);
//		auto c4 = ssh::fourpointcorr(sites, {j,i,j,i}, L, beta, dtau);
//		if (c1 != c1) cout << "c1: " << c1 << endl;
//		if (c2 != c2) cout << "c2: " << c2 << endl;
//		if (c3 != c3) cout << "c3: " << c3 << endl;
//		if (c4 != c4) cout << "c4: " << c4 << endl;


		/*

		< K(b1,t1) > < K(b2,t2 >

		*/
		const std::complex<double> K1 = ssh::g1D(w1,w1,i,j,L,beta,dtau)+ssh::g1D(w1,w1,j,i,L,beta,dtau);
		const std::complex<double> K2 = ssh::g1D(w2,w2,i,j,L,beta,dtau)+ssh::g1D(w2,w2,j,i,L,beta,dtau);
		retval -= K1*K2;

		return std::real(retval);
	}
}


// ------------------------------ OBSERVABLES ---------------------------

class SSHMag
{
	public:
		std::string name, desc;
		template <class StateSpace, class Grid> 
		double measure(const StateSpace& statespace, const Grid& grid) 
		{
			const int N = grid.size();
			double mag = 0.0;

			for (int i=0; i<N; i++) 
			{
				 mag += statespace[i][0];
				}
			return mag/double(N);
		}
		SSHMag() : name("m"), desc("The Magnetization of the SSH Modell") {}
};


class SSHMagSq
{
	public:
		std::string name, desc;
		template <class StateSpace, class Grid> 
		double measure(const StateSpace& statespace, const Grid& grid) 
		{
			const int N = grid.size();
			double mag = 0.0;

			for (int i=0; i<N; i++) 
			{
				 mag += pow(statespace[i][0],2);
				}
			return mag/double(N);
		}

		SSHMagSq() : name("msq"), desc("...") {}
};



#include <math.h> // for debug, provides "isnan"

class SSHTwoPointCorrSpace
{
	public:
	std::string name, desc;
	template <class StateSpace, class Grid>
	std::array<double, 32> measure(const StateSpace& statespace, const Grid& grid)
//	std::vector<double> measure(const StateSpace& statespace, const Grid& grid) // not yet supported
	{

		const int Ls = grid.len;
		const int Lt = grid.lentime;

		const double norml = 1.0 / Ls / Ls;

//		std::vector<double> retval;
//		for (int i=0; i<10; i++) retval.push_back(0);

		std::array<double,32> retval;
		for (int i=0; i<32; i++) retval[i] = 0;

		for (int i=0; i<Lt; i++)
		{
			for (int j=0; j<Ls; j++)
			{
				for (int k=0; k<Ls; k++)
				{
					if (k==j) continue;

					const int idx1 = i*Ls+j;
					const int idx2 = i*Ls+k;

				
					const double x1 = statespace[idx1][0];
					const double x2 = statespace[idx2][0];
					
					//// actual coordinates not necessary for 1D
					// const auto c1 = grid.getcrds(idx1)[0];
					// const auto c2 = grid.getcrds(idx2)[0];
					// auto dist = c1-c2;
					// if (fabs(dist) > 0.5) dist = 1.0 - fabs(dist);

					const int diff = std::abs(j-k) % (Ls/2);

					retval[diff] += norml*x1*x2;
				}
			}
		}

		return retval;
	}
					
	
	SSHTwoPointCorrSpace() : name("corr"), desc("Spatial two-point correlator") {}

};

template <class StateVector>
class SSH_interaction : public Interaction<StateVector> 
{
public:
	SSH_interaction(double m, double dtau)
	{
		this->J = -m/dtau;
	}
	StateVector get (const StateVector& phi) {return phi;};
};




template <class StateVector>
class SSH_onsite : public OnSite<StateVector, double> 
{
	public:
		SSH_onsite(double m, double k, double dtau)
		{
			this->h = m/dtau + k*dtau/2;
		}
		double get (const StateVector& phi) {return dot(phi,phi);}; 
};



template <class StateSpace, class StateVector>
class SSH_multisite
{
	public:
		double k, beta, dtau;

		SSH_multisite(double g, double b, double d) : k(-0.5*g*g), beta(b), dtau(d) {
        };

		// computes the difference already!
		// can be simplified if "suscept" is symmetric

		template <class Lattice>
		double diff (const int rsite,
					const StateVector& svold, 
					const StateVector& svnew, 
					std::vector<int>& nbrs, 
					StateSpace& s,
					Lattice& grid)
		{
			double retval = 0;

			for (int i=0; i<nbrs.size(); i++)
			{

				auto b1 = svnew-svold;
				auto b2 = ssh::suscept(grid, rsite, nbrs[i], beta,dtau);
				auto b3 = s[i];

//				if (b1[0] != b1[0]) cout << "b1: " << b1[0] << endl;
//				if (b2 != b2)       cout << "b2: " << b2    << endl;
//				if (b3[0] != b3[0]) cout << "b3: " << b3[0] << endl;

				auto a1 = dot(b1, mult(b2,b3));
				auto a2 = dot(b3, mult(b2,b1));

//				if (a1 != a1) cout << "a1: " << a1 << endl;
//				if (a2 != a2) cout << "a2: " << a2 << endl;

				retval = retval + a1 + a2;

//				retval +=	dot((svnew-svold), mult(ssh::suscept(grid, rsite, nbrs[i], beta,dtau), s[i]));
//				retval +=	dot(s[i], mult(ssh::suscept(grid, nbrs[i], rsite, beta,dtau), (svnew-svold)));
			}

			retval += dot(svnew, mult(ssh::suscept(grid, rsite, rsite, beta, dtau), svnew));
			retval -= dot(svold, mult(ssh::suscept(grid, rsite, rsite, beta, dtau), svold));

			return retval;
		}
		~SSH_multisite() {}
};







template <class StateVector, class RNG>
class SSH_Initializer
{
	private:
		RNG& rng;
	public:
		SSH_Initializer()   {}
		SSH_Initializer(RNG& rn) : rng(rn) {}

		// specifies how a random new state vector is generated
		// in this case a simple spin flip
		StateVector newsv(const StateVector& svold) 
		{
			StateVector retval(svold); 
			double amp = 0.2;
			double diff = rng.real(-amp, amp);
			retval[0] += diff;
			return retval;
		};
};



// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType = double>
class SSH
{
	public:
		double m, k, g, dtau, betaQM;
		constexpr static int SymD = 1;
		const std::string name;
		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = SSH_Initializer<StateVector, RNG>;

		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta = 1;
		static constexpr uint Ngamma = 1;
		
		SSH(double m, double k, double g, double bQ, int Ltime) : m(m), k(k), g(g), dtau(bQ/double(Ltime)), betaQM(bQ), name("SSH")
		{
            std::cout<<betaQM<<" "<<Ltime<<std::endl;
			interactions[0] = new SSH_interaction<StateVector>(m, dtau); 
			onsite[0] = new SSH_onsite<StateVector>(m, k, dtau);
			multisite[0] = new SSH_multisite<StateVector*,StateVector>(g, betaQM, dtau);
		}
		
		// instantiate interaction terms (requires pointers)
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, double>* onsite[Nbeta];
		SSH_multisite<StateVector*,  StateVector>* multisite[Ngamma];
	
		// instantiate and choose observables
		SSHMag       obs_m;
		SSHMagSq       obs_msq;
		SSHTwoPointCorrSpace obs_corr;
		auto getobs()	{return std::make_tuple(obs_m, obs_msq, obs_corr);}


		// initialize state space
		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for (int i=0; i<grid.size(); i++)
			{
				if (rng.real() > 0.5) statespace[i][0] = 1;
				else statespace[i][0] = -1;
			}
		}


		// using the Wolff cluster algorithm requires to implement
		// the functions 'wolff_coupling' and 'wolff_flip'

		template <class A = bool>
		inline double wolff_coupling(StateVector& sv1, StateVector& sv2, const A a=0) const 
		{
			if (sv1[0] == sv2[0]) return 0.0;
			else return -1.0;
		}

		template <class A = bool>
		inline void wolff_flip(StateVector& sv, const A a=0) const 
		{
			sv[0] *= -1;
		}


};


#endif
