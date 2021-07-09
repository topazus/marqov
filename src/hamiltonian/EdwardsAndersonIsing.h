/* MARQOV - A modern framework for classical spin models on general topologies
 * Copyright (C) 2020-2021, The MARQOV Project
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef EDWARDSANDERSONISING_H
#define EDWARDSANDERSONISING_H
#include <array>
#include <tuple>
#include <string>
#include <complex>
#include <functional>
#include "util/hamparts.h"
#include "../libmarqov/metropolis.h"


// ------------------------------ OBSERVABLES ---------------------------

// Important Notes:

// Spin glass observables typically require local thermal averages
// which means that first, spin are _individually_ averaged over the
// thermal history of the (properly thermalized) system;
// then the average over the system is performed;
// and then a configuration/replica average

// Definitions:
// Katzgraber et. al, PRB 73, 224432 (2006)
// or Parisi, PRL 50.24 (1983): 1946.

// since keeping track of the thermal history of every site 
// or even every pair of sites (in the susceptibilities) is 
// expensive, we implement the absorvables as "running" averages

// note that for the susceptibility is does not match the
// strict definition!

// IMPROVE ME!




/** Edwards-Anderson order parameter */
class EdwardsAndersonOrderParameter
{
	public:
		int counter = 0;
		std::string name;
		std::vector<int> local_sum;

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int size = grid.size();

			// prepare large enough vector where the running thermal average
			// of every spin can be stored
			if (local_sum.size() == 0) 
			{
				local_sum.resize(size);
				for (int i=0; i<size; i++) { local_sum[i] = 0; }
			}

			double retval = 0;

			// compute and store running thermal average of every site
			// and compute current total average 
			counter++;
			for (int i=0; i<size; i++)
			{
				local_sum[i] += statespace[i][0];
				retval += pow(local_sum[i],2);
			}

			// working with running averages requires some normalization
			return retval / double(size) / double(counter) / double(counter);
		}

		EdwardsAndersonOrderParameter() : name("qEA") {}
};



// another attempt at the spin glass susceptibility
// now we actually track everything
// however the phase factors are already included in the system average
// is that correct?
class RealSusceptibility
{
	public:
		double kx;
		std::string name;
		int counter = 0;
		std::vector<std::vector<std::vector<int>>> sisj;
		std::vector<std::vector<std::complex<double>>> phases;

		template <class StateSpace, class Grid>
		double measure(const StateSpace& s, const Grid& grid)
		{
			const int size = grid.size();
			std::complex<double> jj(0,1); 

			// when this function is executed for the first time
			// we do some preparations
			if (this->counter == 0)
			{	
				sisj.resize(size);
				phases.resize(size);
				for (int i=0; i<size; i++)
				{
					sisj[i].resize(size);
					phases[i].resize(size);
				}
			

				// compute and store phase factors
				const int dir = 0; // we consider only the first spatial component
				for (int i=0; i<size; i++) 
				{
					for (int j=0; j<size; j++) 
					{
						const std::vector<double> xi = {grid.crds(i)[dir]};
						const std::vector<double> xj = {grid.crds(j)[dir]};
						auto diff = xi[0] - xj[0];

						if (fabs(diff)>0.5) diff = 1.0 - fabs(diff); // account for PBC
						phases[i][j] = std::exp(kx*diff*jj);
					}
				}
			}


			// store current values
			for (int i=0; i<size; i++) for (int j=0; j<size; j++) sisj[i][j].push_back(s[i][0]*s[j][0]);

			// compute actually observable only every n-th time the function is executed
			std::complex<double> system_average = 0;

			if ((counter-1)%10 == 0)
			{
				// loop over pairs of sites
				for (int i=0; i<size; i++) 
				{
					for (int j=0; j<size; j++) 
					{
						// compute thermal average for current pair s_i*s_j
						double current_site_average = 0;
						const int nmeas = sisj[i][j].size();
						for (int t=0; t<nmeas; t++) current_site_average += sisj[i][j][t];
						current_site_average /= double(nmeas);
						
						system_average += pow(current_site_average,2) * phases[i][j];
					}
				}
			}

			counter++;
			return std::abs(system_average) / double(size);
		}

		RealSusceptibility(double kx, std::string name) : kx(kx), name(name) {}

};
		
//* Spin glass susceptibility */
// compare general remarks above!
class Susceptibility
{
	public:
		int counter = 0;
		double kx;
		std::string name;
		std::vector<int> sum_i;
		std::vector<std::vector<int>> sum_ij;

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int size = grid.size();
			std::complex<double> jj(0,1); 
			const double norml = 1. /  double(size) / double(size) / double(counter) / double(counter);

			// prepare large enough vector where the running thermal average
			// of every pair of spins can be stored
			if (sum_ij.size() == 0) 
			{
				sum_ij.resize(size);
				for (int i=0; i<size; i++)
				{
					sum_ij[i].resize(size);
					for (int j=0; j<size; j++)
					{
						sum_ij[i][j] = 0;
					}
				}
			}

			if (sum_i.size() == 0) 
			{
				sum_i.resize(size);
				for (int i=0; i<size; i++) sum_i[i] = 0;
			}

			std::complex<double> retval = 0;
			counter++;

			for (int i=0; i<size; i++)
			{
				sum_i[i] += statespace[i][0];

				for (int j=0; j<size; j++)
				{
					sum_ij[i][j] += statespace[i][0]*statespace[j][0];
				}
			}


			// perform system average and add proper phase factors
			for (int i=0; i<size; i++)
			{
				for (int j=0; j<size; j++)
				{
					const int dir = 0; // we consider only the first spatial component

					const std::vector<double> xi = {grid.crds(i)[dir]};
					const std::vector<double> xj = {grid.crds(j)[dir]};
					auto diff = xi[0] - xj[0];

					if (fabs(diff)>0.5) diff = 1.0 - fabs(diff); // account for PBC
					std::complex<double> phase = std::exp(kx*diff*jj);
					retval += pow(sum_ij[i][j],2) * phase;
				}
			}

			return norml * std::abs(retval);

			// open questions (TODO):
			// - should the distance vector account for PBC? -> most likely yes
			// - correct normalization -> almost ;)
			// - order of averages correct? -> think so
			// - what to return? absolute value, real part, ...?
		}

		Susceptibility(double kx, std::string name) : kx(kx), name(name) {}
};






class LinkOverlap /// not working so far!!!! TODO
{
	public:
		int counter = 0;
		std::string name;
		std::vector<std::vector<int>> sum_ij;

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int size = grid.size();

			if (sum_ij.size() == 0) 
			{
				sum_ij.resize(size);
				for (int i=0; i<size; i++)
				{
					auto bnds = grid.bnds(0,i);
					auto nbnds = bnds.size();
					sum_ij[i].resize(nbnds);

					for (decltype(nbnds) j=0; j<nbnds; j++) sum_ij[i][j] = 0;
				}
			}
			int nbondstot = 0;

			for (int i=0; i<size; i++)
			{
				auto bnds = grid.bnds(0,i);
				auto nbnds = bnds.size();
				nbondstot += nbnds;
				for (decltype(nbnds) j=0; j<nbnds; j++) sum_ij[i][j] += statespace[i][0]*statespace[j][0];
			}

			counter++;
			const double norml = 1. / double(counter) / double(nbondstot);


			double retval = 0;

			for (int i=0; i<size; i++)
			{
				auto bnds = grid.bnds(0,i);
				auto nbnds = bnds.size();
				for (decltype(nbnds) j=0; j<nbnds; j++) retval += pow(sum_ij[i][j],2);
			}
			return norml * retval;
		}

		LinkOverlap() : name("ql") {}
};




class InternalEnergy /// not working so far!!!! TODO
{
	public:
		int counter = 0;
		std::string name;
		std::vector<std::vector<int>> sum_ij;

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int size = grid.size();

			if (sum_ij.size() == 0) 
			{
				sum_ij.resize(size);
				for (int i=0; i<size; i++)
				{
					auto bnds = grid.bnds(0,i);
					auto nbnds = bnds.size();
					sum_ij[i].resize(nbnds);

					for (decltype(nbnds) j=0; j<nbnds; j++) sum_ij[i][j] = 0;
				}
			}

			for (int i=0; i<size; i++)
			{
				auto bnds = grid.bnds(0,i);
				auto nbnds = bnds.size();
				for (decltype(nbnds) j=0; j<nbnds; j++) sum_ij[i][j] += statespace[i][0]*statespace[j][0];
			}

			const double norml = 1. / double(counter) / double(size);

			counter++;

			double retval = 0;

			for (int i=0; i<size; i++)
			{
				auto bnds = grid.bnds(0,i);
				auto nbnds = bnds.size();
				for (decltype(nbnds) j=0; j<nbnds; j++) retval += sum_ij[i][j] * bnds[j];
			}
			return - norml * retval;
		}

		InternalEnergy() : name("U") {}
};
		



// ------------------------------ INITIALIZER ---------------------------

template <class StateVector, class RNG>
class EdwardsAndersonIsing_Initializer
{
	public:
		EdwardsAndersonIsing_Initializer()   {}
		EdwardsAndersonIsing_Initializer(RNG&) {}

		// specifies how a random new state vector is generated
		// in this case a simple spin flip
		StateVector newsv(const StateVector& svold) 
		{
			StateVector retval(svold); 
			retval[0] = -retval[0];
			return retval;
		};
};


// ------------------------------ HAMILTONIAN ---------------------------

template <class StateVector>
class EdwardsAndersonIsing_interaction
{
public:
    const double& J;
	EdwardsAndersonIsing_interaction(const double& J) : J(J) {}
	StateVector get (const StateVector& phi) {return phi;};
};



/** Hamiltonian for a Edwards-Anderson spin glass */
template <typename SpinType = int>
class EdwardsAndersonIsing
{
	public:

		//  ----  Parameters  ----

		double J;
		static constexpr int SymD = 1;
		const std::string name;


		//  ---- Definitions  -----

		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = EdwardsAndersonIsing_Initializer<StateVector, RNG>;


		//  ----  Hamiltonian terms  ----

		std::array<EdwardsAndersonIsing_interaction<StateVector>*, 1> interactions 
			= {new EdwardsAndersonIsing_interaction<StateVector>(J)};

		EdwardsAndersonIsing(double J) : J(J), 
										 name("EdwardsAndersonIsing"), 
										 obs_chi(0, "chi"), 
										 obs_chiKmin(2.0*M_PI, "chiKmin") {}

		~EdwardsAndersonIsing() {delete interactions[0];}


		//  ----  Observables ----
		
		EdwardsAndersonOrderParameter	obs_qEA;
		RealSusceptibility					obs_chi;
		RealSusceptibility					obs_chiKmin;
		InternalEnergy					obs_U;
        decltype(std::make_tuple(obs_qEA, obs_chi, obs_chiKmin)) observables 
			= {std::make_tuple(obs_qEA, obs_chi, obs_chiKmin)};


		//  ----  Initializer  ----

		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for (decltype(grid.size()) i = 0; i < grid.size(); i++)
			{
				if (rng.real() > 0.5) statespace[i][0] = 1;
				else statespace[i][0] = -1;
			}
		}
};


#endif
