#ifndef SSH_H
#define SSH_H
#include <array>
#include <tuple>
#include <string>
#include <complex>
#include <functional>
#include <math.h> // for debug, provides "isnan"
#include "../hamparts.h"
#include "../metropolis.h"


// ----------------------------------- OBSERVABLES --------------------------------

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




// only 1D so far!!!
class SSHTwoPointCorrSpace
{
	public:
	std::string name, desc;
	template <class StateSpace, class Grid>
	std::array<double, 64> measure(const StateSpace& statespace, const Grid& grid)
//	std::vector<double> measure(const StateSpace& statespace, const Grid& grid) // not yet supported
	{

		const int Ls = grid.len;
		const int Lt = grid.lentime;

		const double norml = 1.0 / Ls / Ls;

//		std::vector<double> retval;
//		for (int i=0; i<10; i++) retval.push_back(0);

		std::array<double,64> retval;
		for (int i=0; i<64; i++) retval[i] = 0;

		for (int i=0; i<Lt; i++)
		{
			for (int j=0; j<Ls; j++)
			{
				const int idx1 = i*Ls+j;

				for (int k=0; k<Ls; k++)
				{
//					if (k==j) continue;

					const int idx2 = i*Ls+k;

					const double x1 = statespace[idx1][0];
					const double x2 = statespace[idx2][0];

					//// actual coordinates not necessary for 1D
					// const auto c1 = grid.getcrds(idx1)[0];
					// const auto c2 = grid.getcrds(idx2)[0];
					// auto dist = c1-c2;
					// if (fabs(dist) > 0.5) dist = 1.0 - fabs(dist);

					int diff = std::abs(j-k);
					if (diff>Ls/2) diff -= Ls/2;

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



// --------------------------------------------------------------------------------
// ============================ Multi-Site Class ==================================
// --------------------------------------------------------------------------------

template <class StateSpace, class StateVector>
class SSH_multisite
{
	public:
		double k, beta, dtau, g;
		int L;
		int ntau;
		double *const gdat;
        std::complex<double>* dexpk;
        double *const ftexp;
		SSH_multisite(double g, double b, double d, int myL) : k(-0.5*g*g), beta(b), dtau(d), L(myL), ntau(std::round(beta/dtau)), g(g),
#ifndef SSH_2D
		gdat(new double[ntau*L]),
		ftexp(new double[L*L])
#else
        gdat(new double[ntau*L*L]),
		ftexp(new double[L*L*L*L])
#endif
		{
			// ntau should be nothing else than Ltime
            dexpk = new std::complex<double>[L];
            for(int d = 0; d < L; ++d)
            {
                dexpk[d] = std::exp(std::complex<double>(0.0, -2.0*M_PI/L*d ));
                for(int k = 0; k < L; ++k)
                    ftexp[d*L + k] = std::real(std::exp(std::complex<double>(0.0, -2.0*M_PI*k/L*d )));
            }
#ifndef SSH_2D
			for(int j = 0; j < L; ++j)
			{
				double k = (j*2)*M_PI/double(L);
				// 1D disperson relation
				const double mu = 0;
				double eps = -2.0*std::cos(k)-mu;
				for(int dt = 0; dt < ntau; ++dt)
				{
					gdat[dt * L + j] = 0.5*std::exp((-beta/2 + dt*dtau)*eps)/std::cosh(beta*eps/2);
				}
			}
#else
			for(int jx = 0; jx < L; ++jx)
			{
				double kx = (jx*2)*M_PI/double(L);
                for(int jy = 0; jy < L; ++jy)
                {
                    double ky = (jy*2)*M_PI/double(L);
                    // 2D disperson relation
                    const double mu = 0;
                    double eps = -2.0*std::cos(kx) - std::cos(ky) - mu;
                    for(int dt = 0; dt < ntau; ++dt)
                    {
                        gdat[dt * L*L + L*jx + jy] = 0.5*std::exp((-beta/2 + dt*dtau)*eps)/std::cosh(beta*eps/2);
                    }
                }
			}

#endif
		}


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
				auto b2 = suscept(grid, rsite, nbrs[i]);
				auto b3 = s[i];

				auto a1 = dot(b1, mult(b2,b3));
				auto a2 = dot(b3, mult(b2,b1));

				retval = retval + a1 + a2;
			}

			retval += dot(svnew, mult(suscept(grid, rsite, rsite), svnew));
			retval -= dot(svold, mult(suscept(grid, rsite, rsite), svold));

			return retval;
		}

		~SSH_multisite() {delete [] gdat; delete [] dexpk; delete [] ftexp; }

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
	
		
		#ifdef SSH_2D

		// Green's function for 2+1 dimensions
		// takes two coordinate vectors (x,y,t)
		template <typename VertexType>
		double green(const VertexType& c1, const VertexType& c2)
		{
			double retval = 0;
			std::complex<double> jj(0,1);

			// space
			auto distx = std::lrint(c1[0]-c2[0]); // account for p.b.c not relevant? 
			auto disty = std::lrint(c1[1]-c2[1]);
			if (distx < 0) distx = L + distx;
			if (disty < 0) disty = L + disty;

			// time
			const double t1 = c1[2];
			const double t2 = c2[2];
			int t1i = floor(t1);
			int t2i = floor(t2);
			int dti = t1i - t2i;
			double signum = 1;
			if (dti < 0) 
			{
				signum = -1; // anti-symmetric behaviour in delta tau
				dti = ntau+dti;
			}

			// compute Greens function in Fourier space
			std::complex<double> dexpkx = dexpk[distx];
			std::complex<double> dexpky = dexpk[disty];
			std::complex<double> expk = 1.0;

			for (int jx = 0; jx < L; ++jx)
			{
				for (int jy = 0; jy < L; ++jy)
				{
					// dispersion relation
					// do the summation
					retval += expk.real() * gdat[dti*L*L + jx*L + jy];
					// increment Fourier transform
					expk *= dexpky; 
				}
				expk *= dexpkx; // increment Fourier transform
			}

			const double norml = 1.0 / pow(2*L,2); // the number of sites per time slice
			return norml*signum*retval;
		}

		#else


		// Green's function for 1+1 dimensions
		// takes two coordinate vectors (x,t)
		template <typename VertexType>
[[gnu::hot, gnu::optimize("fast-math") ]]
		double green(const VertexType& c1, const VertexType& c2)
		{
			// space
			auto dist = std::lrint(c1[0]-c2[0]);
			if (dist < 0) dist = L + dist;
	
			// time
			const double t1 = c1[1];
			const double t2 = c2[1];
			int t1i = floor(t1);
			int t2i = floor(t2);
			int dti = t1i - t2i;
			double signum = 1;
			if (dti < 0) 
			{
				signum = -1; // anti-symmetric behaviour in delta tau
				dti = ntau+dti;
			}
	
			// account for periodic boundaries of the lattice
			// not needed, says Florian
			// if (dti  > 0.5*ntau) dti = ntau - dti;
			// if (dist > 0.5*L)   dist = L - dist;
	
	
			// compute Greens function in Fourier space
                        const double *const __restrict__ dat = ftexp + dist*L;
                        const double *const __restrict__ gdatloc =  gdat + dti*L;
                        double retval = 0;
			for (int j = 0; j < L; ++j)
			{
//                retval = std::fma(dat[j], gdatloc[j], retval);//currently fastest and simplest on my CPU...
 				retval +=  dat[j] * gdat[j];
			}

			const double norml = 1.0/L;
			return norml*signum*retval;
		}
		#endif




	// performs a Wick decomposition
	template <typename VertexType>
	std::complex<double> wick(VertexType v0,VertexType v1, VertexType v2, VertexType v3)
	{
		std::complex<double> retval = 0;
		if (v0==v1 && v1==v2 && v2==v3 && v3==v0)
		{
			retval = green(v0,v0);
		}
		else
		{
			retval = green(v0,v1)*green(v2,v3) - green(v0,v3)*green(v2,v1);
		}
		return retval;
	}



	// the actual susceptibility
	// takes two indices, representing two bonds in the system
	template <class Lattice>
	double suscept(Lattice& grid, int idx1, int idx2)
	{
		std::complex<double> retval = 0;

		auto w1 = grid.getcrds(idx1); // (xstart, xend, time) in 1D / (xstart, xend, ystart, yend, time) in 2D
		auto w2 = grid.getcrds(idx2);

		#ifdef SSH_2D
			const auto t1 = w1[4];
			const auto t2 = w2[4];

//			const auto c1 = wick({w1[0],w1[1],t1}, {w1[2],w1[3],t1}, {w2[0],w2[1],t2}, {w2[2],w2[3],t2});
//			const auto c2 = wick({w1[0],w1[1],t1}, {w1[2],w1[3],t1}, {w2[2],w2[3],t2}, {w2[0],w2[1],t2});
//			const auto c3 = wick({w1[2],w1[3],t1}, {w1[0],w1[1],t1}, {w2[0],w2[1],t2}, {w2[2],w2[3],t2});
//			const auto c4 = wick({w1[2],w1[3],t1}, {w1[0],w1[1],t1}, {w2[2],w2[3],t2}, {w2[0],w2[1],t2});

			const auto c1 = wick<decltype(w1)>({w1[0],w1[2],t1}, {w1[1],w1[3],t1}, {w2[0],w2[2],t2}, {w2[1],w2[3],t2});
			const auto c2 = wick<decltype(w1)>({w1[0],w1[2],t1}, {w1[1],w1[3],t1}, {w2[1],w2[3],t2}, {w2[0],w2[2],t2});
			const auto c3 = wick<decltype(w1)>({w1[1],w1[3],t1}, {w1[0],w1[2],t1}, {w2[0],w2[2],t2}, {w2[1],w2[3],t2});
			const auto c4 = wick<decltype(w1)>({w1[1],w1[3],t1}, {w1[0],w1[2],t1}, {w2[1],w2[3],t2}, {w2[0],w2[2],t2});

//			const std::complex<double> K1 = green({w1[0],w1[1],w1[4]},{w1[2],w1[3],w1[4]}) + green({w1[2],w1[3],w1[4]},{w1[0],w1[1],w1[4]});
//			const std::complex<double> K2 = green({w2[0],w2[1],w2[4]},{w2[2],w2[3],w2[4]}) + green({w2[2],w2[3],w2[4]},{w2[0],w2[1],w2[4]});

			const std::complex<double> K1 = green<decltype(w1)>({w1[0],w1[2],t1},{w1[1],w1[3],t1}) + green<decltype(w1)>({w1[1],w1[3],t1},{w1[0],w1[2],t1});
			const std::complex<double> K2 = green<decltype(w1)>({w2[0],w2[2],t2},{w2[1],w2[3],t2}) + green<decltype(w1)>({w2[1],w2[3],t2},{w2[0],w2[2],t2});
		#else
			const auto t1 = w1[2];
			const auto t2 = w2[2];

			const auto c1 = wick<decltype(w1)>({w1[0],t1}, {w1[1],t1}, {w2[0],t2}, {w2[1],t2});
			const auto c2 = wick<decltype(w1)>({w1[0],t1}, {w1[1],t1}, {w2[1],t2}, {w2[0],t2});
			const auto c3 = wick<decltype(w1)>({w1[1],t1}, {w1[0],t1}, {w2[0],t2}, {w2[1],t2});
			const auto c4 = wick<decltype(w1)>({w1[1],t1}, {w1[0],t1}, {w2[1],t2}, {w2[0],t2});

			const std::complex<double> K1 = green<decltype(w1)>({w1[0],t1},{w1[1],t1}) + green<decltype(w1)>({w1[1],t1},{w1[0],t1});
			const std::complex<double> K2 = green<decltype(w1)>({w2[0],t2},{w2[1],t2}) + green<decltype(w1)>({w2[1],t2},{w2[0],t2});
		#endif

		return std::real(c1+c2+c3+c4-K1*K2);
	}
};





// ------------------------------ Initializer ------------------------------


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
			double amp = 0.7;
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
		double m, k, g, dtau, betaQM, L;
		constexpr static int SymD = 1;
		const std::string name;
		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = SSH_Initializer<StateVector, RNG>;

		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta = 1;
		static constexpr uint Ngamma = 1;
		
		SSH(double m, double k, double g, double bQ, int Ltime, int L) : m(m), k(k), g(g), dtau(bQ/double(Ltime)), L(L), betaQM(bQ), name("SSH")
		{
			interactions[0] = new SSH_interaction<StateVector>(m, dtau); 
			onsite[0] = new SSH_onsite<StateVector>(m, k, dtau);
			multisite[0] = new SSH_multisite<StateVector*,StateVector>(g, betaQM, dtau, L);
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
