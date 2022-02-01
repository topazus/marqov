#ifndef MARQOVSCHEDULER_H
#define MARQOVSCHEDULER_H
/* MIT License
 * 
 * Copyright (c) 2020 - 2022 Florian Goth
 * fgoth@physik.uni-wuerzburg.de
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <vector>
#include <string>
#include <mutex>
#include <algorithm>
#include <tuple>
#include <chrono>

#include "marqovqueue.h"
#include "core.h"
#ifdef MPIMARQOV
#include <type_traits>
#include <mpi.h>
#endif
#include "util/filters.h"

namespace MARQOV
{
    template <typename FPType>
    constexpr auto calcMHratio(FPType en)
    {
        FPType mhratio = 0;
        //prevent under/overflow
        if (en > std::log(std::numeric_limits<FPType>::min()))
        {
            if (en > std::log(std::numeric_limits<FPType>::max()))
                mhratio = std::numeric_limits<FPType>::max();
            else
                mhratio = std::exp(en);
        }
        return mhratio;
    }
    /** The Marqov internal scheduler.
     * 
     * It encapsulates the creation of simulations, the parallel tempering
     * and the distribution across nodes/cores.
     * @tparam Sim a fully specified Marqov type
     */
    template <class Sim>
    class CXX11Scheduler
    {
    public:
        /** Create a full simulation from a parameter.
         * 
         * This gives us the parameters of a simulation and we are responsible for setting everything up.
         * It has a template parameter, but of course all used parameters have to resolve to the same underlying MarqovType.
         * @param p The full set of parameters that are relevant for your problem.
         * @param filter A filter that can be applied before the actual creation of MARQOV.
         */
        template <typename ParamType, typename Callable = decltype(defaultfilter)>
        void createSimfromParameter(ParamType& p, Callable filter = defaultfilter)
        {
            mlogstate.reset();
            auto t = filter(p);//FIXME: I think the filter may not modify the type of parameters anymore.
            bool needswarmup = !Sim::dumppresent(std::get<1>(t));
            int idx = gamekernels.size();
            std::get<1>(t).id = idx+1000;

            auto loadkernel = [&, t]()
            {
                    return makeCore<typename Sim::Lattice, typename Sim::HamiltonianType, typename Sim::RNGT>(t, mutexes.hdf);
            };

            std::function<void(Simstate, int)> gamekernel = [&, t](Simstate mywork, int npt)
            {
                mlogstate.reset();
                 std::cout << "("<<mywork.id<<")>>" << std::endl;
                {
//                     auto sim = makeCore<typename Sim::Lattice, typename Sim::HamiltonianType>(t, mutexes.hdf);
                    auto sim = kernelloaders[mywork.id]();
                    // We loop until the next PT step
                    for(; mywork.npt < npt; ++mywork.npt)
                    {
                        //std::cout<<"Gamelooping on item "<<mywork.id<<" "<<mywork.npt<<std::endl;
                        sim.gameloop();
                    }
                    mlogstate.reset();
                    MLOGRELEASEVERBOSE<<"Final action of id "<<mywork.id<<" : "<<sim.calcAction(sim.statespace)<<std::endl;
                }
                if (mywork.npt < maxpt) // determine whether this itm needs more work
                {
                    //                 std::cout<<"putting item again into workloop"<<std::endl;
                    workqueue.push_back(mywork);
                }
                else
                {
                    //                 std::cout<<"no more work required on "<<mywork.id<<std::endl;
                    masterwork.notify_all();//trigger those waiting for signals from the taskqueue. since we don't push_back anything they would not be notified.
                }
//                 std::cout<<"finished gamekernel closing file"<<std::endl;
            };

            gamekernelmutex.lock();
            gamekernels.push_back(gamekernel);
            kernelloaders.push_back(loadkernel);
            gamekernelmutex.unlock();
            if(needswarmup)
            {
                std::function<void()> warmupkernel = [&, t, idx]
                {
//                     std::cout<<"Beginning warmup of "<<idx<<std::endl;
                 		std::cout << "("<<idx<<")w" << std::flush;
                    {
                        auto sim = makeCore<typename Sim::Lattice, typename Sim::HamiltonianType>(t, mutexes.hdf);
                        sim.init();
                        sim.wrmploop();
                    }
                    //enqueue the next full work item into the workqueue immediately
                    workqueue.push_back(Simstate(idx));
//                     std::cout<<"finished warmup closing file"<<std::endl;
                };
                taskqueue.enqueue(warmupkernel);
            }
            else
            {
                MLOGRELEASEVERBOSE<<"[MARQOV::CXX11Scheduler] Previous step found! Restarting!"<<std::endl;
                workqueue.push_back(Simstate(idx));
            }
        }
        /** This registers an already allocated simulation with us.
         * 
         * DEPRECATED: Untested code path as of now.
         * 
         * @param sim A reference to the sim that already exists.
         * @param warmup a boolean to select whether we require warmup.
         */
        void enqueuesim(Sim& sim, bool warmup = true)
        {
            std::function<void(Simstate, int)> gamekernel = [&](Simstate mywork, int npt)
            {
                // We loop until the next PT step
                for(; mywork.npt < npt; ++mywork.npt)
                {
                    //    std::cout<<"Gamelooping on item "<<mywork.id<<" "<<mywork.npt<<std::endl;
                    sim.gameloop();
                }
                if (mywork.npt < maxpt) // determine whether this itm needs more work
                {
                    //                 std::cout<<"putting item again into workloop"<<std::endl;
                    workqueue.push_back(mywork);
                }
                else
                {
                    //                 std::cout<<"no more work required on "<<mywork.id<<std::endl;
                    masterwork.notify_all();//trigger those waiting for signals from the taskqueue. since we don't push_back anything they would not be notified.
                }
            };
            
            int idx = gamekernels.size();
            gamekernelmutex.lock();
            gamekernels.push_back(gamekernel);
            gamekernelmutex.unlock();

            if (warmup)
            {
                std::function<void()> warmupkernel = [&, idx]
                {
                    sim.init();
                    sim.wrmploop();
                    //enqueue the next full work item into the workqueue immediately
                    workqueue.push_back(Simstate(idx));
                };
                taskqueue.enqueue(warmupkernel);
            }
            else
            {
                workqueue.push_back(Simstate(idx));
            }
        }

        //FIXME...
        static bool exchangepossible(Sim& sima, Sim& simb)
        {
           return true;
        }

        template <class F>
        void createPTplan(F f)
        {
            auto len = gamekernels.size();
            auto adapter = [&f, &len](){auto t = f(); return std::make_pair(t.first%len, t.second%len);};
            std::generate_n(std::back_inserter(ptplan), maxpt, adapter);
            
            for (int i = 0; i < maxpt; ++i)
                std::cout<<ptplan[i].first<<" "<<ptplan[i].second<<std::endl;
        }
        
        /** Start the simulations! GoGoGo...!
         */
        void start()
        {
//             class Seq
//             {
//                 int i = 0;
//                 auto operator()(){return std::make_pair(i,i++);}
//             } seq;
            auto ran = [&] {
                return std::make_pair(rng.integer(), rng.integer());
            };
            mlogstate.reset();
            createPTplan(ran);

            std::cout<<"Starting up master"<<std::endl;
            Simstate itm;
            
            while(!masterstop)
            {
                //              std::cout<<"Master waiting for work"<<std::endl;
                bool busy = false;
                //FIXME: The following wait_for construct hides a bug that occurs if the last notify in the gameloop triggers the master, but the associated task is still running.
                masterwork.wait_for(std::chrono::seconds(10), [&]{
                    busy = workqueue.pop_front(itm);
                    if (!busy)
                        masterstop = nowork();
                    return busy || masterstop;
                });
                if(busy) //there really is sth. to do
                {
//                    std::cout<<"dealing with work"<<std::endl;
                    // check if this sim is selected for PT in this time step. This should usually be the case since we do as many steps as necessary.
                    // in the first time step no pt happens. Hence npt = 1 should disable PT.
                    if((itm.npt > 0) && item_needs_pt(itm.npt, itm.id))
                    {
                        ptstep(itm);
                    }
                    else
                    {//usually triggered at the beginning
                        movesimtotaskqueue(itm);
                    }
                }
                else
                {
                    //test whether there is work lying around somewhere
                    masterstop = nowork();
                }
            }
            //          std::cout<<"Master stopped"<<std::endl;
            //};
            //      taskqueue.enqueue(master);
        }
        void waitforall() {}
        /** Construct Scheduler.
         * 
         * @param maxptsteps How many parallel tempering steps do we do. Defaults to just a single PTstep and hence disables it.
         * @param nthreads how many threads should be used. If not specified defaults to what is reported by the OS.
         * @param id An integer id. This is used by MPI to pass down the rank.
         */
        CXX11Scheduler(int maxptsteps = 1, uint nthreads = 0, int mid = 0) : myid(mid), mlogstate(DEBUG, "marqovmaster_rank"+std::to_string(myid)+".mlog"), maxpt(maxptsteps), masterstop(false), masterwork{},
        workqueue(masterwork),
        taskqueue(((nthreads == 0)?std::thread::hardware_concurrency():nthreads)), rng(time(0) + std::random_device{}())
        {}
        /** Tidy up scheduler.
         * 
         * This frees all resources and waits until all threads have finished.
         */
        ~CXX11Scheduler() {
            if (!nowork() && !masterstop && (taskqueue.tasks_enqueued() > 0) )
            {
                masterwork.wait([&]{
                    return workqueue.is_empty() && (taskqueue.tasks_enqueued() == 0) && (taskqueue.tasks_assigned() == 0);
                });
            }
            mlogstate.reset();
            masterstop = true;
            MLOGVERBOSE<<"PT acceptance: "<<acceptedmoves/static_cast<double>(maxpt)<<std::endl;
        }

//        Scheduler() = delete; //FIXME: If this would be present, the call to Scheduler() would be ambiguous
        CXX11Scheduler(const CXX11Scheduler&) = delete;
        /** Move Constructor
         * 
         * The other object over whose resources we take ownership.
         * Mutexes are a bit odd here. We don't reuse the other mutexes but use and create our own.
         * 
         * @param rhs the other object.
         */
        
        CXX11Scheduler(CXX11Scheduler&& rhs) : myid(rhs.myid), mutexes{}, mlogstate(rhs.mlogstate), maxpt(rhs.maxpt), ptqueue(std::move(rhs.ptqueue)), ptplan{std::move(rhs.ptplan)}, masterstop(rhs.masterstop),
        masterwork{}, workqueue(masterwork), taskqueue{std::move(rhs.taskqueue)}, gamekernels{}
        {
            std::swap(gamekernels, rhs.gamekernels);
            if (rhs.taskqueue.tasks_enqueued() > 0 || rhs.taskqueue.tasks_assigned() > 0)
            {
                mlogstate.reset();
                MLOGRELEASE<<"CXX11Scheduler: invalid assignment\n";
                throw std::runtime_error("[MARQOV::CXX11Scheduler] invalid assignment");
            }
        }
        CXX11Scheduler& operator=(const CXX11Scheduler&) = delete;
        CXX11Scheduler& operator=(CXX11Scheduler&& ) = delete;
    private:
        /**
         * Simstate helper class
         * This class encapsulates the parallel tempering state of a single sim.
         */
        struct Simstate
        {
           Simstate() = default;
           Simstate(int i) : id(i), npt(0) {}
           Simstate(int i, int np) : id(i), npt(np) {}
           int id = -1;///< my id
           int statespacesize = 0;
           int npt = -100;///< which will be my next parallel tempering step.
        };

        /** Find the parallel tempering exchange partner of the given id.
         * 
         * @param id find the next partner that this id has.
         */
        auto findpartner(uint id) const
        {
            return std::find_if(ptqueue.cbegin(), ptqueue.cend(), [&id](const Simstate& itm){return itm.id == static_cast<int>(id);});
        }
        
        /** Check whether a work item needs parallel tempering in this time step.
         * 
         * @param pttime the parallel tempering time step that the item is at.
         * @param id id for this work item
         * @return true if the plan requires a pt exchange with someone at this point in time, else false.
         */
        bool item_needs_pt(int pttime, int id) const
        {
            return (ptplan[pttime].first == id || ptplan[pttime].second == id);
        }
        
        /** Test whether there is work available.
         * 
         * @return true if no task is working and no work is to be executed by a task and no sim is to moved to the taskqueue.
         */
        bool nowork() {return workqueue.is_empty() && taskqueue.tasks_assigned() == 0 && taskqueue.tasks_enqueued() == 0;}
        
        /** Do a parallel tempering step. 
         * 
         * This function is called when the current simulation is up for a parallel tempering (PT) step.
         * If its partner is already waiting we do the parallel tempering, if not we get moved into 
         * a queue and wait for a partner.
         * @param itm The Sim which is chosen for PT
         */
        void ptstep(Simstate itm) {
                mlogstate.reset();
                MLOGRELEASEVERBOSE<<"Parallel Tempering!"<<std::endl;
                MLOGRELEASEVERBOSE<<"itm.id "<<itm.id<<" itm.npt "<<itm.npt<<std::endl;
                MLOGRELEASEVERBOSE<<"Expected pairing for this time step: "<<ptplan[itm.npt].first<<" "<<ptplan[itm.npt].second<<std::endl;
                MLOGRELEASEVERBOSE<<"ptstep begin"<<std::endl;
            if (ptplan[itm.npt].first == ptplan[itm.npt].second)
            {//this can happen and is not prevented
                MLOGRELEASEVERBOSE<<"[MARQOV::Scheduler] self swap in PT..."<<std::endl;
                movesimtotaskqueue(itm);
                return;
            }
            int partner = ptplan[itm.npt].first;
            if (partner == itm.id) partner = ptplan[itm.npt].second;//it must be the other. no exchanges with myself
            auto partnerinfo = findpartner(partner);

            if ((partnerinfo != ptqueue.cend()) && (itm.npt == partnerinfo->npt) )
            {// partner is at the same stage, hence we can PT exchange

                MLOGRELEASEVERBOSE<<"Partner "<<partner<<" found in queue"<<std::endl;
                {
                    auto sima = kernelloaders[itm.id]();
                    auto simb = kernelloaders[partnerinfo->id]();
                    double mhratio = calcprob(sima, simb);
                    if(rng.real() < std::min(1.0, mhratio)){
                        exchange_statespace(sima, simb);
                        ++acceptedmoves;
                    }
                }
                ptqueue.erase(partnerinfo);
                //put both sims back into the taskqueue for more processing until their next PT step
                movesimtotaskqueue(itm);
                movesimtotaskqueue(Simstate(partner, itm.npt));
            }
            else
            {//we have to wait for the PT partner
                MLOGRELEASEVERBOSE<<"Partner "<<partner<<" not in queue"<<std::endl;
                ptqueue.push_back(itm);
            }
        }
        /** Take simulation and move it to the workqueue.
         * 
         * This determines how many steps have to be done until the next PTstep
         * and moves the simulation into the taskqueue where the gameloop is executed.
         * @param itm The simulation that gets further worked on.
         */
        void movesimtotaskqueue(Simstate itm)
        {
            int newnpt = findnextnpt(itm.id, itm.npt);
            mlogstate.reset();
            MLOGRELEASEVERBOSE<<"Putting a new item with id "<<itm.id<<" with npt = "<<itm.npt <<" until npt = "<< newnpt<<" into the taskqueue"<< std::endl;
            taskqueue.enqueue(
                [&, itm, newnpt]{gamekernels[itm.id](itm, newnpt);} //Get the required kernel from the array of gamekernels and execute it.
            );
        }
        /** Determine the next PT step.
         *
         * @param idx simulation id to check
         * @param curnpt current PT time
         * @return the next PT step where this simulation is selected for PT.
         */
        uint findnextnpt(int idx, uint curnpt) const
        {
            uint retval = curnpt+1;
            while ((retval < static_cast<uint>(maxpt)) && (ptplan[retval].first != idx) && (ptplan[retval].second != idx))
            {
                ++retval;
            }
            return retval;
        }
        template <class T>
        double calcprob (T& sima, T& simb) const
        {
//            auto sima = kernelloaders[ida]();
//            auto simb = kernelloaders[idb]();
            double actionaa = sima.calcAction(sima.statespace);
            double actionbb = simb.calcAction(simb.statespace);
            double actionab = sima.calcAction(simb.statespace);
            double actionba = simb.calcAction(sima.statespace);
            std::cout<<"Action of id 1"<<" : "<<actionaa<<std::endl;
            std::cout<<"Action of id 2"<<" : "<<actionbb<<std::endl;
            
            std::cout<<"Action of id 1"<<" with other statespace"<<" : "<<actionab<<std::endl;
            std::cout<<"Action of id 2"<<" with other statespace"<<" : "<<actionba<<std::endl;
            //calculate actual metropolis transition ratio from the propability densities
//             double mhratio = std::exp(-actionab)*std::exp(-actionba)/std::exp(-actionaa)/std::exp(-actionbb);
            double endiff = actionaa - actionab + actionbb - actionba;
            double mhratio = calcMHratio(endiff);
//             double mhratio = std::exp(actionaa - actionab + actionbb - actionba);
            std::cout<<"M-H Ratio: "<<mhratio<<std::endl;
            return mhratio;
        }
        template <class T>
        void exchange_statespace(T& sima, T& simb) 
        {
    	    swap(sima.statespace, simb.statespace);
        }
        /**
         * This class collects mutexes that synchronize I/O.
         */
        struct GlobalMutexes
        {
            std::mutex hdf;///< Lock for the HDF5 I/O since the library for C++ is not thread-safe.
            std::mutex io;///< Lock for the rest?
        } mutexes;
        int myid{0};///< An idea to distinguish giles from multiple schedulers.
        MLogState mlogstate;
        int acceptedmoves = 0;
        int maxpt; ///< how many pt steps do we do
        std::vector<Simstate> ptqueue; ///< here we collect who is waiting for its PT partner
        std::vector<std::pair<int, int> > ptplan; ///< An array of who exchanges with whom in each step
        bool masterstop; ///< A global flag to denote that the master has decided to stop.
        ThreadPool::Semaphore masterwork; ///< The semaphore that triggers the master process
        ThreadPool::ThreadSafeQueue<Simstate> workqueue; ///< This is the queue where threads put their finished work and the master does PT.
        std::mutex simvectormutex; ///< A mutex to protect accesses to the simvector which could be invalidated by the use of push_back
        std::mutex gamekernelmutex; ///< A mutex to protect accesses to the gamekernels which could be invalidated by the use of push_back
        ThreadPool::Queue taskqueue; ///< This is the queue where threads pull their work from.
        std::vector<std::function<void(Simstate, int)> > gamekernels; ///< prefabricated workitems that get executed to move a simulation forward.
        std::vector<std::function<Sim(void)> > kernelloaders;
        RNGCache<std::mt19937_64> rng{static_cast<std::mt19937_64::result_type>(0)};
    };

#ifndef MPIMARQOV
    template <class MarqovType>
    using Scheduler = CXX11Scheduler<MarqovType>;
#else
    template <class Sim>
    class MPIScheduler
    {
    public:
        /** This gives us the parameters of a simulation and we are responsible for setting everything up.
         * It has a template parameter, but of course all used parameters have to resolve to the same underlying MarqovType.
         * FIXME: It is expected that all MPI ranks execute the same code until here!!! That makes it easier to have valid data on every node...
         * @param p The full set of parameters that are relevant for your Problem
         * @param filter A filter that can be applied before the actual creation of MARQOV
         */
        template <typename ParamType, typename Callable = decltype(defaultfilter)>
        void createSimfromParameter(ParamType& p, Callable filter = defaultfilter)
        {
            if (myrank == rrctr)
                myScheduler.createSimfromParameter(p, filter);
            if (myrank == MASTER)
            {
                rrctr = (rrctr + 1) % nr_nodes;
            }
            MPI_Bcast(&rrctr, 1, MPI_INT, MASTER, marqov_COMM);
        }
        void enqueuesim(Sim& sim)
        {
            if (myrank == rrctr)
                myScheduler.enqueuesim(sim);
            if (myrank == MASTER)
            {
                rrctr = (rrctr + 1) % nr_nodes;
            }
            MPI_Bcast(&rrctr, 1, MPI_INT, MASTER, marqov_COMM);
        }
        /** Start execution on all nodes.
         */
        void start()
        {
            myScheduler.start();
        }
        void waitforall() {}
        /** Construct MPI Scheduler.
         * 
         * We expect MPI to be initialized beforehand. The user should feel that he is writing MPI code.
         * @param maxptsteps How many parallel tempering steps do we do
         * @param nthreads how many threads should be used. If not specified defaults to what is reported by the OS.
         */
        MPIScheduler(int maxptsteps = 1, uint nthreads = 0) : rrctr(0), maxpt(maxptsteps), marqov_COMM(std::move(createMPICommunicator())), nr_nodes(getavailableMPIhosts()), myrank(getmyMPIrank()), myScheduler(maxptsteps, nthreads, myrank)
        {}
         /** Move Copy Constructor
         * 
         * The other object over whose resources we take ownership.
         * 
         * @param rhs the other object
         */
        MPIScheduler(MPIScheduler&& rhs) : rrctr(rhs.rrctr), marqov_COMM(rhs.marqov_COMM), myScheduler(std::move(rhs.myScheduler)), myrank(rhs.myrank), nr_nodes(rhs.nr_nodes), maxpt(rhs.maxpt) {}
        MPIScheduler(const MPIScheduler&) = delete;
        ~MPIScheduler() = default;
        MPIScheduler& operator=(const MPIScheduler&) = delete;
        MPIScheduler& operator=(MPIScheduler&& ) = delete;
    private:
        static auto createMPICommunicator()
        {
            int mpi_inited;
            MPI_Initialized(&mpi_inited);
            if (!mpi_inited)
                throw("MPI not initialized!");
            MPI_Comm tmp_COMM;///< our own MPI communicator.
            MPI_Comm_dup(MPI_COMM_WORLD, &tmp_COMM);
            return tmp_COMM;
        }
        int getavailableMPIhosts(){
            int t;
            MPI_Comm_size(marqov_COMM, &t);
            return t;
        }
        int getmyMPIrank(){
            int t;
            MPI_Comm_rank(marqov_COMM, &t);
            return t;
        }
        int rrctr;
        int maxpt; ///< how many pt steps do we do
        MPI_Comm marqov_COMM;///< our own MPI communicator.
        static constexpr int MASTER = 0;
        int nr_nodes; ///< how many nodes are we actually executed on.
        int myrank; ///< my MPI rank
        CXX11Scheduler<Sim> myScheduler;//MPI starts the program parallely on every node(if properly executed), hence every node needs his CXX11scheduler.
    };

    template <class MarqovType>
    using Scheduler = MPIScheduler<MarqovType>;
#endif
    /** A helper class to figure out the type of the scheduler.
     * 
     * @tparam Hamiltonian the type of the Hamiltonian.
     * @tparam Lattice The type of the lattice.
     * @tparam Parameters The type of the parameters. We instantiate Hamiltonian
     *                    and probably lattice in Marqov, hence the params.
     * @tparam RNGType The type of the Random Number Generator. Must conform to C++11 Interface.
     *                 If you intend to use it for a non-STL RNG, have a look at @see rngcache.h
     *                 to have a proper name dumped for it in the HDF5 files.
     */
    template <class Hamiltonian, class Lattice, class Parameters, class RNGType = std::mt19937_64>
    struct GetSchedulerType
    {
        typedef std::mutex& mtxref;
        typedef decltype(makeCore<Lattice, Hamiltonian, RNGType>(std::declval<Parameters>(), std::declval<mtxref>())) MarqovType;
        typedef Scheduler<MarqovType> MarqovScheduler; ///< Holds the type of a scheduler for these simulations.
    };

    /** A helper function to figure out the type of the scheduler and actually create it.
     *
     * @tparam Hamiltonian The Hamiltonian to use for the simulations.
     * @tparam Lattice The lattice to use for the simulations.
     * @tparam RNGType The Random Number Generator to use.
     * @tparam Parameters the tuple that makes up the arguments for a simulation.
     * @tparam Ts a parameter pack of the remaining arguements to the scheduler.
     *
     * @param args A dummy parameter to get rid of template magic on the user end.
     * @param ts The remaining arguments for the scheduler.
     *
     */
    template <class Hamiltonian, class Lattice, class RNGType = std::mt19937_64, class Parameters, typename ...Ts>
    auto makeScheduler(const Parameters& args, Ts&&... ts)
    {
        return typename GetSchedulerType<Hamiltonian, Lattice, Parameters>::MarqovScheduler(std::forward<Ts>(ts)...);
    }
};
#endif
