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

#include "timetracker.h"
#include "helpers.h"
#include "vectorhelpers.h"
#include "registry.h"
#include "systemtools.h"
#include "replicate.h"
#include "svmath.h"
#include "filters.h"
#include "embedder.h"
#include "marqovscheduler.h"
#include "util.h"

// Geometry
#include "geometry/regular_lattice.h"
#include "geometry/ssh_lattice.h"
#include "geometry/grid.h"
#include "geometry/neighbourclass.h"
#include "geometry/io.h"


// Hamiltonians
#include "hamiltonian/Phi4.h"

using namespace MARQOV;


//#define STATIC_BOUNDARY


int main()
{
    
    
    std::cout<<"MARQOV Copyright (C) 2020-2021, The MARQOV Project contributors"<<std::endl;
    std::cout<<"This program comes with ABSOLUTELY NO WARRANTY."<<std::endl;
    std::cout<<"This is free software, and you are welcome to redistribute it under certain conditions."<<std::endl;

	// read config files
	RegistryDB registry("../src/config", "ini");

	// remove old output and prepare new one
	std::string outbasedir = registry.Get<std::string>("mc.ini", "IO", "outdir" );
	std::string logbasedir = registry.Get<std::string>("mc.ini", "IO", "logdir" );

	//FIXME: NEVER DELETE USER DATA
	std::string command;
	command = "rm -r " + outbasedir;
	system(command.c_str());
	command = "rm -r " + logbasedir;
	system(command.c_str());

	makeDir(outbasedir);
	makeDir(logbasedir);


	std::string ham = "Phi4";


	auto beta   = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
	auto lambda = registry.Get<std::vector<double> >("mc.ini", ham, "lambda");
	auto mass   = registry.Get<std::vector<double> >("mc.ini", ham, "mass");
	auto hp = cart_prod(beta, lambda, mass);
	
	// Parameters
	const auto name = registry.Get<std::string>("mc.ini", "General", "Hamiltonian" );
	auto nreplicas  = registry.Get<std::vector<int>>("mc.ini", name, "rep" );
	const auto nL   = registry.Get<std::vector<int>>("mc.ini", name, "L" );
	const auto dim  = registry.Get<int>("mc.ini", name, "dim" );
	
	
	// Typedefs
	typedef Phi4<double,double> Hamiltonian;
	typedef Hyperbolic Lattice;

    Hyperbolic lat(42,42);
    typedef typename std::tuple<Lattice&, MARQOV::Config, std::tuple<double,double,double> > ParameterType;
	typedef typename GetSchedulerType<Hamiltonian, Lattice, ParameterType>::MarqovScheduler SchedulerType;

 	SchedulerType sched(1,1);

	// prepare output
	std::string outpath = outbasedir+"/hyperbolic/";
	makeDir(outpath);
	
	// Monte Carlo parameters
	MARQOV::Config mp(outpath);
	mp.setnmetro(1);
	mp.setncluster(0);
	mp.setwarmupsteps(0);
	mp.setgameloopsteps(1000);

	// lattice parameters

	// form parameter triple and replicate
	auto params  = finalize_parameter(lat, mp, hp);
    auto rparams = replicator(params, 1);

	// schedule simulations
 	for (auto p: rparams) sched.createSimfromParameter(p, defaultfilter);
 	
	sched.start(); // run!
}
