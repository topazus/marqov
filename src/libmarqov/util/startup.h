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

#ifndef STARTUP_H
#define STARTUP_H
#include <vector>
#include <tuple>
//#include "../core.h"
//#include "registry.h"

void create_config_select(std::string HamiltonianName)
{
	const auto filename = std::string{"./config/select.ini"};
	std::ofstream select(filename);
	select << "[General]" << '\n' << "Hamiltonian = " << HamiltonianName << '\n' << "[END]" << std::endl;
	select.close();

}

void print_welcome_message()
{
	std::cout << std::endl << "WELCOME TO MARQOV!" << std::endl << std::endl;
}

template <class RegistryType>
void check_registry_availability(RegistryType& reg, std::string HamiltonianName, std::string dir = "./config")
{
	try
	{
		reg.init("./config", "ini");
	}

	catch(Registry_Exception& re)
	{
		std::cout << "[MARQOV::main] No configuration directory found! ";
		std::cout << "Assuming you're starting this MARQOV binary for the first time!" << std::endl;
		std::cout << "[MARQOV::main] To get you going we will generate and populate ";
		std::cout << "a configuration directory locally under ./config" << std::endl;

		makeDir("./config"); 
		const auto filename = std::string{"./config/select.ini"};

		if(!fileexists(filename))
		{
			create_config_select(HamiltonianName);
			reg.init("./config", "ini");
		}
		else
		{
			std::cout << "[MARQOV::main] "<< filename <<" already exists, but is not usable.";
			std::cout << "I would overwrite its content, hence I'm terminating now"<<std::endl;
			throw;
		}
	}
}


template <class RegistryType>
void check_registry_file_exists(RegistryType& reg, std::string name, std::string dir = "./config")
{

	std::string fn{name + ".ini"};
	try
	{
		auto dim  = reg.Get<int>(fn, name, "dim" );
	}
	catch(Registry_cfgfile_not_found_Exception& e)
	{
		std::cout<<"[MARQOV::main] Unable to find Heisenberg config! Generating new one in ./config/"<<fn<<std::endl;
		std::ofstream cfg("./config/" + fn);
		cfg << "[Heisenberg]\n" << "L = 12\n" << "rep = 4\n" << "dim = 3\n";
		cfg << "beta = 0.67,0.68\n" << "J = -1.0\n\n";
		cfg << "[MC]\n" << "nmetro = " << 2 << "\nnclusteramp = "<< 0.5 << "\nnclusterexp = " << 1; 
		cfg << "\nwarmupsteps = " << 500 << "\nmeasuresteps = "<< 5000 << "\n";
		cfg << "[IO]\n" << "outdir = ./out\n" << "[END]" << std::endl;
		cfg.close();

		reg.init("./config"); 
	}
}





#endif
