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

#ifndef FILTERS_H
#define FILTERS_H
#include <string>
#include <tuple>
	// -------------------- filters --------------------

	// filter to determine output file path and name
	// the filter _must_ set the outname

	auto defaultfilter = [](auto p)
	{
		auto& mp = std::get<1>(p);		// Monte Carlo params
		auto& hp = std::get<2>(p);	// Hamiltonian params
		
		std::string str_repid = std::to_string(mp.repid);
		std::string str_beta  = "beta"+std::to_string(std::get<0>(hp));
		std::string str_extf  = "extf"+std::to_string(std::get<3>(hp));
        
		mp.outname = str_beta+"_"+str_extf+"_"+str_repid;
		return p;
	};

// 	auto defaultfilter_triple = [](auto p)
// 	{
//        	auto& lp = p.first;
//        	auto& mp = p.second;
//        	auto& hp = p.third;
// 
//  		auto str_repid = std::to_string(mp.repid);
// 		auto str_beta  = "beta"+std::to_string(std::get<0>(hp));
// 		auto str_L     = std::to_string(std::get<0>(lp));
// 
// 		mp.outname = str_beta+"_"+str_repid;
// 
// 		return p;
// 	};

	auto xxzfilter = [](auto p)
	{	
		auto& mp = std::get<1>(p);		// Monte Carlo params
		auto& hp = std::get<2>(p);	// Hamiltonian params
	
		std::string str_repid = std::to_string(mp.repid);
		std::string str_beta  = "beta"+std::to_string(std::get<0>(hp));
		std::string str_extf  = "extf"+std::to_string(std::get<1>(hp));
		mp.outname = str_beta+"_"+str_extf+"_"+str_repid;

		return p;
	};


	auto sshfilter = [](auto p)
	{
		auto& mp = p.first;		// Monte Carlo params
		auto& hp = p.second;	// Hamiltonian params
		
		std::string str_repid = std::to_string(mp.repid);
		std::string str_k     = "k"+std::to_string(std::get<2>(hp));
		std::string str_dtau  = "dtau"+std::to_string(std::get<3>(hp));
		mp.outname = mp.outname+"_"+str_k+"_"+str_dtau+"_"+str_repid;
		return p;
	};

#endif
