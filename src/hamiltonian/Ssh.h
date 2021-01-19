#ifndef SSH_H
#define SSH_H
#include <array>
#include <tuple>
#include <string>
#include <complex>
#include <functional>
#include "../hamparts.h"
#include "../metropolis.h"

namespace ssh
{


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



// =========================== Multi-Site Class =================================

template <class StateSpace, class StateVector>
class SSH_multisite
{
	public:
		double k, beta, dtau, g;
		int L;
		double* gdat;
		int ntau;
		SSH_multisite(double g, double b, double d, int myL) : k(-0.5*g*g), beta(b), dtau(d), L(myL), g(g)
		{
			ntau = std::round(beta/dtau); // ntau should be nothing else than Ltime
			gdat = new double[ntau*L];

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



			/*
			auto c1 = std::vector<double>{0,0};
			auto c2 = std::vector<double>{11,0};
			auto c3 = std::vector<double>{7,0};
			auto c4 = std::vector<double>{42,0};
			std::vector<std::vector<double>> arg1;
			arg1.push_back(c1);
			arg1.push_back(c2);
			arg1.push_back(c3);
			arg1.push_back(c4);

			std::vector<int> arg2 = {1,1,1,1};

			cout << wick(arg1, arg2) << endl;
			*/

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

		~SSH_multisite() {delete [] gdat;}

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
	
		
		double green(std::vector<double> c1, std::vector<double> c2)
		{
			std::complex<double> retval = 0;
			std::complex<double> jj(0,1);
		
			const double x1 = c1[0];
			const double y1 = c1[1];
			const double x2 = c2[1];
			const double y2 = c2[1];

			double distx = x1-x2;
			double disty = y1-y2;

			if (distx < 0) distx = L + distx;
			if (disty < 0) disty = L + disty;

			const double t1 = c1[4];
			const double t2 = c2[4];
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
			std::complex<double> dexpkx = std::exp(-2*M_PI/L*jj*distx);
			std::complex<double> dexpky = std::exp(-2*M_PI/L*jj*disty);
			std::complex<double> expkx = 1.0;
			std::complex<double> expky = 1.0;


			for (int jx = 0; jx < L; ++jx)
			{
				for (int jy = 0; jy < L; ++jy)
				{
					double disp = - 2*cos(2*M_PI*jx/L) - 2*cos(2*M_PI*jy/L);
					retval += expkx * expky * exp(dti*dtau*disp) * fermi(beta*disp);
					expky *= dexpky;
				}

				expkx *= dexpkx; // increment Fourier transform
			}

			auto retv = signum*retval.real()/double(pow(2*L,2));

			return retv;



		}


		// Green's function for a regular hypercubic lattice in 1D (i.e. a chain)
		/*
		double green(std::vector<double> c1, std::vector<double> c2, int sign1, int sign2)
		{
			std::complex<double> retval = 0;
			std::complex<double> jj(0,1);
		
			double step = 2*M_PI/double(L);
		
			// spatial coordinates
			const double r1 = c1[0] + sign1*0.5; // account for p.b.c not relevant? (we only use relative distances...)
			const double r2 = c2[0] + sign2*0.5;
			double dist = r1-r2;
			if (dist < 0) dist = L + dist;
	
			// temporal coordinates
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
//			if (dti  > 0.5*ntau) dti = ntau - dti;
//			if (dist > 0.5*L)   dist = L - dist;
	
	
			// compute Greens function in Fourier space
			std::complex<double> dexpk = std::exp(-2*M_PI/L*jj*dist);
			std::complex<double> expk = 1.0;
			for (int j = 0; j < L; ++j)
			{
				retval += expk.real()*gdat[dti*L + j];
				expk   *= dexpk; //loop-carried dependency. breaks vectorization.
			}
			auto retv = signum*retval.real()/double(L);

			return retv;
		}
		*/


/*

	 j1/i2      b2     j2/i3      b3      j3/i4      b4
  ------x-----------------x------------------x-----------------x---------- // original model: bonds

  ---------------o-----------------o------------------o----------------    // our model: bonds represented as sites


   e.g.  r(b2(j2)) = r(b3(i3))

*/

	// performs a Wick decomposition
	std::complex<double> wick(std::vector<double> v0, std::vector<double> v1, std::vector<double> v2, std::vector<double> v3)
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


//	std::complex<double> wick(std::vector<std::vector<double>> s, std::vector<int> e)
//	{
//		std::complex<double> retval = 0;
//
//		// the square of the number operator of fermions is the number operator
//		if (s[0]==s[1] && s[1]==s[2] && s[2]==s[3] && e[0]==e[1] && e[1]==e[2] && e[2]==e[3])
//		{
//			retval = green(s[0], s[1], e[0], e[1]);
//		}
//
//		// Wick decomposition
//		else 
//		{
//			retval = green(s[0],s[1],e[0],e[1])*green(s[2],s[3],e[2],e[3]) 
//				  - green(s[0],s[3],e[0],e[3])*green(s[2],s[1],e[2],e[1]);
//		}
//
//		return retval;
//	}



	template <class Lattice>
	double suscept(Lattice& grid, int idx1, int idx2)
	{
		std::complex<double> retval = 0;

		auto w1 = grid.getcrds(idx1);
		auto w2 = grid.getcrds(idx2);
		const int L = grid.len;

		const auto c1 = wick({w1[0],w1[1],w1[4]}, {w1[2],w1[3],w1[4]}, {w2[0],w2[1],w2[4]}, {w2[2],w2[3],w2[4]});
		const auto c2 = wick({w1[0],w1[1],w1[4]}, {w1[2],w1[3],w1[4]}, {w2[2],w2[3],w2[4]}, {w2[0],w2[1],w2[4]});
		const auto c3 = wick({w1[2],w1[3],w1[4]}, {w1[0],w1[1],w1[4]}, {w2[0],w2[1],w2[4]}, {w2[2],w2[3],w2[4]});
		const auto c4 = wick({w1[2],w1[3],w1[4]}, {w1[0],w1[1],w1[4]}, {w2[2],w2[3],w2[4]}, {w2[0],w2[1],w2[4]});

		retval = c1+c2+c3+c4;

		const std::complex<double> K1 = green({w1[0],w1[1],w1[4]},{w1[2],w1[3],w1[4]}) + green({w1[2],w1[3],w1[4]},{w1[0],w1[1],w1[4]});
		const std::complex<double> K2 = green({w2[0],w2[1],w2[4]},{w2[2],w2[3],w2[4]}) + green({w2[2],w2[3],w2[4]},{w2[0],w2[1],w2[4]});

		retval -= K1*K2;

		return std::real(retval);

	}


//	template <class Lattice>
//	double suscept(Lattice& grid, int idx1, int idx2)
//	{
//
//		auto w1 = grid.getcrds(idx1);
//		auto w2 = grid.getcrds(idx2);
//		const int L = grid.len;
//	
//		// signs
//		const int i = -1;
//		const int j = +1;
//
//		std::complex<double> retval = 0;
//
//		/*
//
//		< K(b1,t1) K(b2,t2) > = 	  
//							  < c†(i1) c(j1) c†(i2) c(j2) > 
//							+ < c†(i1) c(j1) c†(j2) c(i2) > 
//							+ < c†(j1) c(i1) c†(i2) c(j2) > 
//							+ < c†(j1) c(i1) c†(j2) c(i2) > // temporal indices supressed; b1 = (i1,j1), b2 = (i2,j2)
//						  =
//						  	  <i1 j1> <i2 j2> - <i1 i2> <i2 i1>
//						  	+ <i1 j1> <j2 i2> - <i1 i2> <j2 j1>
//						  	+ <j1 i1> <i2 j2> - <j1 j2> <i2 i1>
//						  	+ <j1 i1> <j2 i2> - <j1 j2> <j2 j1> // operators supressed
//		*/
//
//		std::vector<std::vector<double>> sites = {w1,w1,w2,w2};
//
//		const auto c1 = wick(sites, {i,j,i,j});
//		const auto c2 = wick(sites, {i,j,j,i});
//		const auto c3 = wick(sites, {j,i,i,j});
//		const auto c4 = wick(sites, {j,i,j,i});
//
//		retval = c1+c2+c3+c4;
//
//		/*       < K(b1,t1) > < K(b2,t2 >    	*/
//		const std::complex<double> K1 = green(w1,w1,i,j)+green(w1,w1,j,i);
//		const std::complex<double> K2 = green(w2,w2,i,j)+green(w2,w2,j,i);
//		retval -= K1*K2;
//
//		return std::real(retval);
//	}

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
