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
#include "../hamparts.h"
#include "../metropolis.h"


// ------------------------------ OBSERVABLES ---------------------------

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

			if (local_sum.size() == 0) 
			{
				local_sum.resize(size);
				for (int i=0; i<size; i++) { local_sum[i] = 0; }
			}

			double retval = 0;

			counter++;
			for (int i=0; i<size; i++)
			{
				local_sum[i] += statespace[i][0];
				retval += pow(local_sum[i],2);
			}

			return retval / double(size) / double(counter) / double(counter);
		}

		EdwardsAndersonOrderParameter() : name("qEA") {}
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
		



//* Spin glass susceptibility */
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

		std::array<EdwardsAndersonIsing_interaction<StateVector>*, 1> interactions = {new EdwardsAndersonIsing_interaction<StateVector>(J)};
        std::array<OnSite<StateVector, int>*, 0> onsite;
        std::array<FlexTerm<StateVector*,  StateVector>*, 0> multisite;

		EdwardsAndersonIsing(double J) : J(J), 
										 name("EdwardsAndersonIsing"), 
										 obs_chi(0, "chi"), 
										 obs_chiKmin(2.0*M_PI, "chiKmin") {}

		~EdwardsAndersonIsing() {delete interactions[0];}


		//  ----  Observables ----
		
		EdwardsAndersonOrderParameter	obs_qEA;
		Susceptibility					obs_chi;
		Susceptibility					obs_chiKmin;
		InternalEnergy					obs_U;
        decltype(std::make_tuple(obs_qEA, obs_chi, obs_chiKmin, obs_U)) observables = {std::make_tuple(obs_qEA, obs_chi, obs_chiKmin, obs_U)};


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





// ------------------------------ SPECIALIZATIONS ---------------------------

#endif
