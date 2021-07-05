#ifndef MARQOVSCHEDULER_H
#define MARQOVSCHEDULER_H
/* MIT License
 * 
 * Copyright (c) 2020-2021 Florian Goth
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
#include "util/filters.h"

namespace MARQOV
{
    /** The Marqov internal scheduler.
     * 
     * It encapsulates the creation of simulations, the parallel tempering
     * and the distribution across nodes/cores.
     * @tparam Sim a fully specified Marqov type
     */
    template <class Sim>
    class Scheduler
    {
    private:
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
            auto t = filter(p);//FIXME: I think the filter may not modify the type of parameters anymore.
            bool needswarmup = !Sim::dumppresent(std::get<1>(t));
            int idx = gamekernels.size();
            
            auto loadkernel = [&, t]()
            {
                    return makeCore<typename Sim::Lattice, typename Sim::HamiltonianType>(t, mutexes.hdf);
            };
            
            std::function<void(Simstate, int)> gamekernel = [&, t](Simstate mywork, int npt)
            {
//                 std::cout<<"Beginning gamekernel"<<std::endl;
                {
//                     auto sim = makeCore<typename Sim::Lattice, typename Sim::HamiltonianType>(t, mutexes.hdf);
                    auto sim = kernelloaders[mywork.id]();
                    // We loop until the next PT step
                    for(; mywork.npt < npt; ++mywork.npt)
                    {
                        //    std::cout<<"Gamelooping on item "<<mywork.id<<" "<<mywork.npt<<std::endl;
                        sim.gameloop();
                    }
                    std::cout<<"Final action of id "<<mywork.id<<" : "<<sim.calcAction(sim.statespace)<<std::endl;
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
                std::cout<<"[MARQOV::Scheduler] Previous step found! Restarting!"<<std::endl;
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
        /** Start the simulations! GoGoGo...!
         */
        void start()
        {
            //create dummy data for the ptplan
            for (int i = 0; i < maxpt; ++i)
                ptplan.emplace_back(i, i+1 );
            
            std::cout<<"Starting up master"<<std::endl;
            Simstate itm;
            
            for (int i = 0; i < maxpt; ++i)
                std::cout<<ptplan[i].first<<" "<<ptplan[i].second<<std::endl;
            
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
                    std::cout<<"dealing with work"<<std::endl;
                    if(ptplan[itm.npt].first == itm.id || ptplan[itm.npt].second == itm.id) // check if this sim is selected for PT in this time step. This should usually be the case since we do as many steps as necessary
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
        /** Construct Scheduler
         * 
         * @param maxptsteps How many parallel tempering steps do we do.
         * @param nthreads how many threads should be used. If not specified defaults to what is reported by the OS.
         */
        Scheduler(int maxptsteps, uint nthreads = 0) : maxpt(maxptsteps), masterstop(false), masterwork{},
        workqueue(masterwork),
        taskqueue(((nthreads == 0)?std::thread::hardware_concurrency():nthreads))
        {}
        /** Tidy up scheduler.
         * 
         * This frees all resources and waits until all threads have finished.
         */
        ~Scheduler() {
            if (!nowork() && !masterstop && (taskqueue.tasks_enqueued() > 0) )
            {
                masterwork.wait([&]{
                    return workqueue.is_empty() && (taskqueue.tasks_enqueued() == 0) && (taskqueue.tasks_assigned() == 0);
                });
            }
            masterstop = true;
        }
    private:
        /**
         * Simstate helper class
         * This class encapsulates the parallel tempering state of a single sim.
         */
        struct Simstate
        {
           Simstate() : id(-1), npt(-100) {}
           Simstate(int i) : id(i), npt(0) {}
           Simstate(int i, int np) : id(i), npt(np) {}
           int id;///< my id
           int npt;///< which will be my next parallel tempering step.
        };
        
        /**
         * This class collects mutexes that synchronize I/O.
         */
        struct GlobalMutexes
        {
            std::mutex hdf;///< Lock for the HDF5 I/O since the library for C++ is not thread-safe.
            std::mutex io;///< Lock for the rest?
        } mutexes;
        /** Find the parallel tempering exchange partner of the given id.
         * 
         * @param id find the next partner that this id has.
         */
        auto findpartner(uint id)
        {
            return std::find_if(ptqueue.cbegin(), ptqueue.cend(), [&id](const Simstate& itm){return itm.id == static_cast<int>(id);});
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
                    std::cout<<"Parallel Tempering!"<<std::endl;
                    std::cout<<"itm.id "<<itm.id<<" itm.npt "<<itm.npt<<std::endl;
                    std::cout<<"Expected pairing for this time step: "<<ptplan[itm.npt].first<<" "<<ptplan[itm.npt].second<<std::endl;
            
            int partner = ptplan[itm.npt].first;
            if (partner == itm.id) partner = ptplan[itm.npt].second;//it must be the other. no exchanges with myself
            auto partnerinfo = findpartner(partner);
            
            if ((partnerinfo != ptqueue.cend()) && (itm.npt == partnerinfo->npt) )
            {// partner is at the same stage, hence we can PT exchange
                            std::cout<<"Partner "<<partner<<" found in queue"<<std::endl;
                ptqueue.erase(partnerinfo);
                double mhratio = calcprob(itm.id, partnerinfo->id);
                exchange();
                //put both sims back into the taskqueue for more processing until their next PT step
                movesimtotaskqueue(itm);
                movesimtotaskqueue(Simstate(partner, itm.npt));
            }
            else
            {//we have to wait for the PT partner
                std::cout<<"Partner "<<partner<<" not in queue"<<std::endl;
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
            std::cout<<"Putting a new item with id "<<itm.id<<" with npt = "<<itm.npt <<" until npt = "<< newnpt<<" into the taskqueue"<< std::endl;
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
        uint findnextnpt(int idx, uint curnpt)
        {
            uint retval = curnpt+1;
            while ((retval < static_cast<uint>(maxpt)) && (ptplan[retval].first != idx) && (ptplan[retval].second != idx))
            {
                ++retval;
            }
            return retval;
        }
        
        int maxpt; ///< how many pt steps do we do
        std::vector<Simstate> ptqueue; ///< here we collect who is waiting for its PT partner
        std::vector<std::pair<int, int> > ptplan; ///< An array of who exchanges with whom in each step
        bool masterstop; ///< A global flag to denote that the master has decided to stop.
        ThreadPool::Semaphore masterwork; ///< The semaphore that triggers the master process
        ThreadPool::ThreadSafeQueue<Simstate> workqueue; ///< This is the queue where threads put their finished work and the master does PT.
        std::mutex simvectormutex; ///< A mutex to protect accesses to the simvector which could be invalidated by the use of push_back
        std::mutex gamekernelmutex; ///< A mutex to protect accesses to the gamekernels which could be invalidated by the use of push_back
        std::vector<Sim*> simvector; ///< An array for the full state of the simulations.
        ThreadPool::Queue taskqueue; ///< This is the queue where threads pull their work from.
        std::vector<std::function<void(Simstate, int)> > gamekernels; ///< prefabricated workitems that get executed to move a simulation forward.
        std::vector<std::function<Sim(void)> > kernelloaders;
        double calcprob(int ida, int idb)
        {
            auto sima = kernelloaders[ida]();
            auto simb = kernelloaders[idb]();
            double enaa = sima.calcAction(sima.statespace);
            double enbb = simb.calcAction(simb.statespace);
            double enab = sima.calcAction(simb.statespace);
            double enba = simb.calcAction(sima.statespace);
            std::cout<<"Action of id "<<ida<<" : "<<enaa<<std::endl;
            std::cout<<"Action of id "<<idb<<" : "<<enbb<<std::endl;
            
            std::cout<<"Action of id "<<ida<<" with other statespace"<<" : "<<enab<<std::endl;
            std::cout<<"Action of id "<<idb<<" with other statespace"<<" : "<<enba<<std::endl;
            //calculate actual metropolis transition ratio from the propability densities
            double mhratio = std::exp(-enab)*std::exp(-enba)/std::exp(-enaa)/std::exp(-enbb);
            std::cout<<"M-H Ratio: "<<mhratio<<std::endl;
            return mhratio;
        }
        //FIXME fill this function for proper PT
        void exchange() {}
    };

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
};
#endif
