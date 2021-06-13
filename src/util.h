#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <vector>

#include "registry.h"
#include "lattice/regular_hypercubic.h"


/** Helper to execute a series of simulations on regular hypercubic lattices
*
* @tparam Hamiltonian the type of the Hamiltonian
* @tparam Params the type of the parameter space
* @tparam Callable the type of the filter
* @param reg the registry
* @param outbasedir	the base directory of the simulation output
* @param hp the Hamiltonian parameters
* @param filter the filter
*/
template <class Hamiltonian, class Params, class Callable>
void RegularLatticeLoop(RegistryDB& reg, const std::string outbasedir, const std::vector<Params>& hp, Callable filter)
{
	
	// Typedefs
	typedef typename std::tuple<RegularHypercubic&, MARQOV::Config, Params> ParameterType;
	typedef typename MARQOV::GetSchedulerType<Hamiltonian, RegularHypercubic, ParameterType>::MarqovScheduler SchedulerType;


	// Parameters
	const auto name = reg.Get<std::string>("mc.ini", "General", "Hamiltonian" );
	auto nreplicas  = reg.Get<std::vector<int>>("mc.ini", name, "rep" );
	const auto nL   = reg.Get<std::vector<int>>("mc.ini", name, "L" );
	const auto dim  = reg.Get<int>("mc.ini", name, "dim" );


	// Number of threads
	int nthreads = 0;
	try 
	{
		nthreads = reg.template Get<int>("mc.ini", "General", "threads_per_node" );
	}
	catch (const Registry_Key_not_found_Exception&) 
	{
		std::cout<<"threads_per_node not set -> automatic"<<std::endl;
	}

	// Prepare Geometry
	if (nreplicas.size() == 1) { for (decltype(nL.size()) i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }
	std::vector<RegularHypercubic> latts;
	for (std::size_t j=0; j<nL.size(); j++) latts.emplace_back(nL[j], dim);

	// Init Scheduler
	SchedulerType sched(1, nthreads);

	for (std::size_t j=0; j<nL.size(); j++)
	{
		// prepare
		int L = nL[j];
		cout << endl << "L = " << L << endl << endl;

		std::string outpath = outbasedir+"/"+std::to_string(L)+"/";

		MARQOV::Config mp(outpath);
		mp.setnmetro(5);
		mp.setncluster(int(L/2));
		mp.setwarmupsteps(500);
		mp.setgameloopsteps(3000);

		makeDir(mp.outpath);
		
		auto params = finalize_parameter(latts[j], mp, hp);
 		auto rparams = replicator(params, nreplicas[j]);
 
		for(auto p : rparams) sched.createSimfromParameter(p, filter);
	}
	sched.start();
}

#endif
