#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <vector>

#include "query_threads.h"
#include "registry.h"
#include "../../lattice/regular_hypercubic.h"

template <class Hamiltonian, class Params, class Callable, int dim>
void FixedRegularLatticeLoop(const std::vector<Params>& hp, Callable filter, const std::vector<int>& nL, const std::vector<int>& nreplicas, const MARQOV::Config& mcdefault, double nclusteramp, int nclusterexp, int nthreads)
{
    // Typedefs
    using Lattice = FixedDimRegularHypercubic<dim>;
	typedef typename std::tuple<Lattice&, MARQOV::Config, Params> ParameterType;
	typedef typename MARQOV::GetSchedulerType<Hamiltonian, Lattice, ParameterType>::MarqovScheduler SchedulerType;
	std::vector<Lattice> latts;
	for (std::size_t j=0; j<nL.size(); j++) latts.emplace_back(nL[j]);

	// Init Scheduler
	SchedulerType sched(1, nthreads);

	for (std::size_t j=0; j<nL.size(); j++)
	{
		// prepare
		int L = nL[j];
		std::cout << std::endl << "L = " << L << std::endl << std::endl;

		MARQOV::Config mp(mcdefault);
        mp.outpath += "/"+std::to_string(L)+"/";
		mp.setncluster(int(nclusteramp*pow(L,nclusterexp)));

		makeDir(mp.outpath);
		
		auto params = finalize_parameter(latts[j], mp, hp);
 		auto rparams = replicator(params, nreplicas[j]);
 
		for(auto p : rparams) sched.createSimfromParameter(p, filter);
	}
	sched.start();
}

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
void RegularLatticeLoop(RegistryDB& reg, std::string configfile, std::string name, const std::string outbasedir, const std::vector<Params>& hp, Callable filter)
{
	// Parameters
	auto nreplicas  = reg.Get<std::vector<int>>(configfile, name, "rep" );
	const auto nL   = reg.Get<std::vector<int>>(configfile, name, "L" );
	const auto dim  = reg.Get<int>(configfile, name, "dim" );

	const auto nclusteramp   = reg.Get<double>(configfile, "MC", "nclusteramp");
	const auto nclusterexp   = reg.Get<int>(configfile, "MC", "nclusterexp");
	const auto nmetro        = reg.Get<int>(configfile, "MC", "nmetro");
	const auto warmupsteps   = reg.Get<int>(configfile, "MC", "warmupsteps");
	const auto measuresteps  = reg.Get<int>(configfile, "MC", "measuresteps");

	const auto nthreads      = number_of_threads_per_node(reg, configfile);

	int loglevel = DEBUG;
	try
	{
	    loglevel  = reg.Get<int>(configfile, "MC", "loglevel");
	}
	catch(Registry_Exception& re)
	{
	;//log level is debug if not explicitly set.
	}
	MARQOV::Config mcdefault(outbasedir);
		mcdefault.setnmetro(nmetro);
		mcdefault.setwarmupsteps(warmupsteps);
		mcdefault.setgameloopsteps(measuresteps);
		mcdefault.setloglevel(loglevel);



	// Prepare Geometry
	if (nreplicas.size() == 1) { for (decltype(nL.size()) i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }

	switch(dim)
    {
        case 1:
            FixedRegularLatticeLoop<Hamiltonian, Params, Callable, 1>(hp, filter, nL, nreplicas, mcdefault, nclusteramp, nclusterexp, nthreads);
            break;
        case 2:
            FixedRegularLatticeLoop<Hamiltonian, Params, Callable, 2>(hp, filter, nL, nreplicas, mcdefault, nclusteramp, nclusterexp, nthreads);
            break;
        case 3:
            FixedRegularLatticeLoop<Hamiltonian, Params, Callable, 3>(hp, filter, nL, nreplicas, mcdefault, nclusteramp, nclusterexp, nthreads);
            break;
        case 4:
            FixedRegularLatticeLoop<Hamiltonian, Params, Callable, 4>(hp, filter, nL, nreplicas, mcdefault, nclusteramp, nclusterexp, nthreads);
            break;
        default:
        {
            // Typedefs
            typedef typename std::tuple<RegularHypercubic&, MARQOV::Config, Params> ParameterType;
            typedef typename MARQOV::GetSchedulerType<Hamiltonian, RegularHypercubic, ParameterType>::MarqovScheduler SchedulerType;
            std::vector<RegularHypercubic> latts;
            for (std::size_t j=0; j<nL.size(); j++) latts.emplace_back(nL[j], dim);

            // Init Scheduler
            SchedulerType sched(1, nthreads);

            for (std::size_t j=0; j<nL.size(); j++)
            {
                // prepare
                int L = nL[j];
                std::cout << std::endl << "L = " << L << std::endl << std::endl;

                MARQOV::Config mp(mcdefault);
                mp.outpath += "/"+std::to_string(L)+"/";
                mp.setncluster(int(nclusteramp*pow(L,nclusterexp)));

                makeDir(mp.outpath);

                auto params = finalize_parameter(latts[j], mp, hp);
                auto rparams = replicator(params, nreplicas[j]);
 
                for(auto p : rparams) sched.createSimfromParameter(p, filter);
            }
            sched.start();
            break;
        }
    }
}

#endif
