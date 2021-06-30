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

#include <array>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <tuple>
#include <iomanip>

using std::cout;
using std::endl;
using std::flush;
using std::ofstream;

// MARQOV
#include "libmarqov/libmarqov.h"

// Lattices
#include "lattice/graph_from_csv.h"

// Hamiltonians
#include "hamiltonian/MassiveScalarField.h"
#include "hamiltonian/MassiveScalarField2.h"

using namespace MARQOV;


int main()
{

    std::cout<<"MARQOV Copyright (C) 2020-2021, The MARQOV Project contributors"<<std::endl;
    std::cout<<"This program comes with ABSOLUTELY NO WARRANTY."<<std::endl;
    std::cout<<"This is free software, and you are welcome to redistribute it under certain conditions."<<std::endl;
    

    
	// remove old output and prepare new one
	std::string outbasedir = "../out";
//	std::string command;
//	command = "rm -r " + outbasedir;
//	system(command.c_str());
//	makeDir(outbasedir);

	

	// Parameters
	std::vector<double> beta = {1.0};
	std::vector<double> K = {0.50315223}; // geometric factor
	std::vector<double> sqrtg = {1.0471975511}; // geometric factor
	std::vector<double> mass = { -10, -1, -0.1, 0, 0.1, 1, 10};	// mass squared

	auto hp = cart_prod(beta, K, sqrtg, mass);
	
	const int nthreads = 2;	
	const int nreplicas = 1;

	
	// Typedefs
	typedef MassiveScalarField Hamiltonian;
	typedef GraphFromCSV Lattice;

	// Lattice
	int nlayers = 8; 
	std::string latfile = "/home/schrauth/hyperbolic-7-3-"+std::to_string(nlayers)+".csv";
    GraphFromCSV lat(latfile);



	// Scheduler
    typedef typename std::tuple<Lattice&, MARQOV::Config, std::tuple<double,double,double,double> > ParameterType;
	typedef typename GetSchedulerType<Hamiltonian, Lattice, ParameterType>::MarqovScheduler SchedulerType;
 	SchedulerType sched(1,nthreads);

	// Output
	std::string outpath = outbasedir+"/test-"+std::to_string(nlayers)+"/";
	makeDir(outpath);
	
	// Monte Carlo
	MARQOV::Config mp(outpath);
	mp.setnmetro(10);
	mp.setncluster(0);
	mp.setwarmupsteps(0);
	mp.setgameloopsteps(200);

	// lattice parameters
	// form parameter triple and replicate
	auto params  = finalize_parameter(lat, mp, hp);
    auto rparams = replicator(params, nreplicas);

	// schedule simulations
 	for (auto p: rparams) sched.createSimfromParameter(p, hyperfilter);
 	
	// run!
	sched.start();

}
