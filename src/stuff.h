/* This file is part of MARQOV:
 * A modern framework for classical spin models on general topologies
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

#include <vector>
		// ----------------- consistency check ------------------

		// perform consistency check according to 
		// Hasenbusch, J. Phys. A: Math. Gen. 34 8221 (2001)
		
		void perform_consistency_check(std::vector<int>& checkidxs)
		{
			std::vector<std::vector<double>> subcheck;

			for (int k=0; k<checkidxs.size(); k++)
			{
				const int checkidx = checkidxs[k];
				std::vector<double> subsubcheck;

				const auto checksite = statespace[checkidx];
				const auto nbrs = this->grid.getnbrs(0, checkidx);

				for (int i = 0; i < nbrs.size(); ++i)
				{
					const auto currentnbr = statespace[nbrs[i]];
					subsubcheck.push_back(dot(checksite,currentnbr));
				}

				const double selfdot = dot(checksite,checksite);

				subsubcheck.push_back(selfdot);
				subsubcheck.push_back((selfdot-1)*selfdot);

				subcheck.push_back(subsubcheck);
			}

			check.push_back(subcheck);

			// monitor memory consumption
			const int nsites = checkidxs.size();
			const int nmeasure = check.size(); 
			const int ncol = check[0][0].size(); 

			constexpr static double check_GB_limit = 4.0;

			if (nsites*nmeasure*ncol > check_GB_limit*1024*1024*1024/8) 
			{
				throw std::overflow_error(
				"\n Running out of memory during consistency check! \n Decrease either lattice size or number of measurements!");
			}
		}


		// evaluate check and display results
		void finalize_consistency_check()
		{
			const int SymD = std::tuple_size<StateVector>::value;
			const int ncol = 8;
			const int nmeasure = check.size();
			const int nsites = check[0].size();

			std::vector<double> sum(ncol,0);

			// compute averages in each column
			for (int k=0; k<nmeasure; k++)
			{
				for (int i=0; i<nsites; i++)
				{
					for (int j=0; j<ncol; j++)
					{
						sum[j] += check[k][i][j];
					}
				}
			}
			
			for (int j=0; j<ncol; j++) 
			{
				sum[j] = sum[j] / double(nmeasure) / double(nsites);
				std::cout << sum[j] << " ";
			}
			std::cout << std::endl;

			// summation formula
			double retval = 0.5*ham.beta*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]) - sum[6] - 2*ham.lambda*sum[7] + 0.5*SymD;
			std::cout << retval << "\n\n";

		}
					

		// specificy site indices which will enter the check
		// default: all (good statistic, but requires a somewhat large amount of memory)
		void prepare_consistency_check(std::vector<int>& checkidxs)
		{
			for (int i=0; i<this->grid.size(); i++)
			{
				checkidxs.push_back(i);
			}
		}


		
		// ----------------- consistency check end ------------------




		double elementaryMCstep();
	    
        void gameloop()
		{
			prepare_consistency_check(checkidxs);

			double avgclustersize = 0;
			for (int k=0; k < this->mcfg.gli; k++)
			{

				if (this->mcfg.id == 0) std::cout << "." << std::flush;
				for (int i=0; i < this->mcfg.gameloopsteps/10; ++i)
				{
					avgclustersize += elementaryMCstep();
					auto obs = ham.getobs();
					marqov_measure(obs, statespace, this->grid);
					perform_consistency_check(checkidxs);
				}
			}

			if (this->mcfg.id == 0) std::cout << "|\n" << avgclustersize/this->mcfg.gameloopsteps << std::endl;
			finalize_consistency_check();
		}
