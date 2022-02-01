/* This file is part of MARQOV:
 * A modern framework for classical spin models on general topologies
 * Copyright (C) 2021-2022, The MARQOV Project
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
#include <string>
#include <fstream>
#include <iostream>
#include "registry.h"
#include "query_threads.h"
#include "systemtools.h"

/** Creates the file select.ini
*
* @param HamiltonianName content of the file
* @param dir directory where the file is created
*/
void create_config_select(const std::string& HamiltonianName, const std::string& dir)
{
	const auto filename = std::string{dir+"/select.ini"};
	std::ofstream select(filename);
	select << "[General]" << '\n' << "Hamiltonian = " << HamiltonianName << '\n' << "[END]" << std::endl;
	select.close();

}


/** Print the default MARQOV startup message
*/
void print_startup_message()
{
	std::cout<<"MARQOV Copyright (C) 2020-2022, The MARQOV Project contributors"<<std::endl;
	std::cout<<"This program comes with ABSOLUTELY NO WARRANTY."<<std::endl;
	std::cout<<"This is free software, and you are welcome to redistribute it under certain conditions."<<std::endl;
}


/** Print a friendly welcome message
*/
void print_welcome_message()
{
	std::cout << std::endl << "WELCOME TO MARQOV!" << std::endl << std::endl;
}


/** Check whether the folder "config" exists and whether there is a file "select.ini"
*
* @param reg a reference to the registry
* @param name the name of the Hamiltonian
* @param dir the desired config directory
*/
void check_registry_availability(RegistryDB& reg, const std::string& name, const std::string dir = "./config")
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

		makeDir(dir); 
		const auto filename = std::string{dir+"/select.ini"};

		if(!fileexists(filename))
		{
			create_config_select(name,dir);
			reg.init(dir, "ini");
		}
		else
		{
			std::cout << "[MARQOV::main] "<< filename <<" already exists, but is not usable.";
			std::cout << "I would overwrite its content, hence I'm terminating now"<<std::endl;
			throw;
		}
	}
}


/** Check whether config file exists for specific Hamiltonian.
*
* This function checks whether the registry has configuration information 
* for the specified Hamiltonian. If there is none, and a suitable "rule" 
* is available, it will be created.
*
* @param reg a reference to the registry
* @param name the name of the Hamiltonian
* @param dir the desired config directory
*/
void check_registry_file_exists(RegistryDB& reg, const std::string& name, const std::string dir = "./config")
{

	std::string fn{name + ".ini"};
	try
	{
		auto dim  = reg.template Get<int>(fn, name, "dim" );
	}
	catch(Registry_cfgfile_not_found_Exception& e)
	{

		if (name == "Heisenberg")
		{
			std::cout<<"[MARQOV::main] Unable to find Heisenberg config! Generating new one in ./config/"<<fn<<std::endl;
			std::ofstream cfg(dir + fn);
			cfg << "[Heisenberg]\n" << "L = 12\n" << "rep = 4\n" << "dim = 3\n";
			cfg << "beta = 0.67,0.68\n" << "J = -1.0\n\n";
			cfg << "[MC]\n" << "nmetro = " << 2 << "\nnclusteramp = "<< 0.5 << "\nnclusterexp = " << 1; 
			cfg << "\nwarmupsteps = " << 500 << "\nmeasuresteps = "<< 5000 << "\n";
			cfg << "[IO]\n" << "outdir = ./out\n" << "[END]" << std::endl;
			cfg.close();
		}
		else if (name == "Ising")
		{
			std::cout<<"[MARQOV::main] Unable to find Ising config! Generating new one in ./config/" <<fn<< std::endl;
			std::ofstream cfg(dir + fn);
			cfg << "[Ising]\n" << "L = 32\n" << "rep = 4\n" << "dim = 2\n";
			cfg << "beta = 0.67,0.68\n" << "J = -1.0\n\n";
			cfg << "[MC]\n" << "nmetro = " << 2 << "\nnclusteramp = "<< 0.5 << "\nnclusterexp = " << 1; 
			cfg << "\nwarmupsteps = " << 500 << "\nmeasuresteps = "<< 5000 << "\n";
			cfg << "[IO]\n" << "outdir = ./out\n" << "[END]" << std::endl;
			cfg.close();
		}
		else
		{
			std::cout << "[MARQOV::main] Unable to find suitable config file!" << std::endl;
			std::cout << "[MARQOV::main] No rule available to auto-create one!" << std::endl;
			throw;
		}

		reg.init("./config"); 
	}
}


/** Find out if a string starts with something.
 * @param longword we search in this string
 * @param shortword we look for this
 * @return truen if longword strarts with shortword, else false
 */
bool startswith(const std::string& longword, const std::string& shortword) noexcept
{
    return longword.find(shortword) == 0;
}


/** Check the vailidity of the replica configuration.
 * throws if invalid.
 * 
 * @param nr number of replicas
 * @param nL amount of lattice simulations
 */
void checkreplicaconfig(int nr, int nL)
{
    if ((nr != nL) && (nr != 1)) throw std::invalid_argument("[MARQOV] Invalid replica configuration!");
}


/** Removes previous simulations.
 * 
 * @param outbasedir The folder that we remove entirely
 */
void tidyupoldsims(const std::string& outbasedir)
{
	// delete previous output // fixme: don't do that by default!
#ifdef MPIMARQOV
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if(myrank == 0)
    {
#endif
        std::cout<<"[MARQOV::main] Erasing previous data!!!!"<<std::endl;
        std::string command = "rm -r " + outbasedir;
        system(command.c_str());
        makeDir(outbasedir);
#ifdef MPIMARQOV
    }
#endif
}


/** Print some nice information
 * 
 * @param registry Where to get the data from.
 * @param ham the hamiltonian
 */
void printInfoandcheckreplicaconfig(RegistryDB& registry, const std::string& ham)
{
    const auto dim 	      = registry.Get<int>(ham+".ini", ham, "dim" );
	const auto nreplicas  = registry.Get<std::vector<int>>(ham+".ini", ham, "rep" );
	const auto nreplicas_str = registry.Get<std::string>(ham+".ini", ham, "rep" );
	const auto nL  	      = registry.Get<std::vector<int>>(ham+".ini", ham, "L" );
	const auto nLs 	      = registry.Get<std::string>(ham+".ini", ham, "L" );

	std::cout << std::endl;
	std::cout << "Hamiltonian: \t" << ham << std::endl;
	std::cout << "Dimension: \t" << dim << std::endl;
	std::cout << "Lattice sizes:\t" << nLs << std::endl;
	std::cout << "Replicas:\t" << nreplicas_str << std::endl;
    
    checkreplicaconfig(nreplicas.size(), nL.size());
}


/** Utility function to write the common part of the config file
 * @param os the stream associated with the new config file.
 * @param nmetro 
 * @param nclusteramp
 * @param nclusterexp
 * @param warmupsteps
 * @param measuresteps
 */
void createcfgfooter(std::ostream& os, int nmetro, double nclusteramp, int nclusterexp, int warmupsteps, int measuresteps)
{
    os << "[MC]\n"<<"nmetro = " << nmetro;
	os << "\nnclusteramp = " << nclusteramp << "\nnclusterexp = " << nclusterexp;
	os << "\nwarmupsteps = " << warmupsteps << "\nmeasuresteps = " << measuresteps << "\n";
    os << "[IO]\n"<<"outdir = ./out\n" << "[END]" << std::endl;
}


#endif
