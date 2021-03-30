#ifndef MARQOVSCHEDULER_H
#define MARQOVSCHEDULER_H
/*
 * MIT License
 * 
 * Copyright (c) 2020 Florian Goth
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

namespace MARQOV
{
    
    template <class Cont, class Tuple1, class Tuple2, std::size_t... I>
    constexpr auto emplace_from_tuple_impl(Cont&& cont, Tuple1&& t1, MARQOV::Config&& mc, std::mutex& mtx, Tuple2&& t2, std::index_sequence<I...> )
    {
        return cont.emplace_back(
            std::forward<Tuple1>(t1), 
                                 std::forward<MARQOV::Config>(mc), mtx,
                                 std::get<I>(std::forward<Tuple2>(t2))...);
    }
    
    //c++17 make_from_tuple from cppreference adapted for emplace
    template <class Cont, class Tuple1, class Tuple2>
    constexpr auto emplace_from_tuple(Cont&& cont, Tuple1&& t1, MARQOV::Config&& mc, std::mutex& mtx, Tuple2&& t2 )
    {
        return emplace_from_tuple_impl(
            cont, 
            std::forward<Tuple1>(t1), 
                                       std::forward<MARQOV::Config>(mc), mtx,
                                       std::forward<Tuple2>(t2),
                                       std::make_index_sequence<std::tuple_size<std::remove_reference_t<Tuple2>>::value>{});
    }
    
    template <class T, class Tuple1, class Tuple2, std::size_t... I>
    constexpr T* ptr_from_tuple_impl(Tuple1&& t1, MARQOV::Config&& mc, std::mutex& mtx, Tuple2&& t2, std::index_sequence<I...> )
    {
        return new T(std::forward<Tuple1>(t1),
                     std::forward<MARQOV::Config>(mc), 
                     mtx,
                     std::get<I>(std::forward<Tuple2>(t2))...
        );
    }
    
    template <class T, class Tuple1, class Tuple2>
    constexpr T* ptr_from_tuple(Tuple1&& t1, MARQOV::Config&& mc, std::mutex& mtx, Tuple2&& t2)
    {
        return ptr_from_tuple_impl<T>(std::forward<Tuple1>(t1),
                                      std::forward<MARQOV::Config>(mc), 
                                      mtx,
                                      std::forward<Tuple2>(t2),
                                      std::make_index_sequence<std::tuple_size<std::remove_reference_t<Tuple2>>::value>{});
    }
    
    template<class ... Ts> struct sims_helper {};
    
    template <class H,  class L, class HArgstuple, size_t... S>
    struct sims_helper<H, L, HArgstuple, std::index_sequence<S...> >
    {
        typedef decltype(MARQOV::makeCore<H>(std::declval<L>(),
                                               std::declval<MARQOV::Config>(), std::declval<std::mutex>(),
                                               std::declval<typename std::tuple_element<S, HArgstuple>::type>()...
        )) MarqovType;
    };
    
    template <class ... Ts>
    struct sims_helper2 {};
    
    template <class Hamiltonian, class Lattice, class LArgs, class HArgs>
    struct sims_helper2<Hamiltonian, Lattice, Triple<LArgs, MARQOV::Config, HArgs> >
    {
        typedef decltype(MARQOV::makeCore<Hamiltonian, Lattice>(std::declval<MARQOV::Config>(), std::declval<std::mutex>(),
                                                                  std::declval<std::pair<LArgs, HArgs>& >()
        )) MarqovType;
        template <typename T>
        static void emplacer(std::vector<MarqovType>& sims, T&  t, std::mutex& mtx)
        {
            emplace_from_tuple(sims, t.first, std::forward<MARQOV::Config>(t.second), mtx, t.third);
        }
        
        template <typename T>
        static MarqovType* creator(std::mutex& mtx, T&  t)
        {
            return ptr_from_tuple<MarqovType>(t.first, std::forward<MARQOV::Config>(t.second), std::forward<decltype(mtx)>(mtx), t.third);
        }
    };
    
    template <class Hamiltonian, class Lattice, class HArgs>
    struct sims_helper2<Hamiltonian, Lattice, std::pair<MARQOV::Config, HArgs> >
    {
        static constexpr std::size_t tsize = std::tuple_size<typename std::remove_reference<HArgs>::type>::value;
        typedef std::make_index_sequence<tsize> HArgSequence;
        typedef typename sims_helper<Hamiltonian, Lattice, HArgs, HArgSequence>::MarqovType MarqovType;
        
        template <typename T>
        static void emplacer(std::vector<MarqovType>& sims, T& t, std::mutex& mtx)
        {
            emplace_from_tuple(sims, 
                               std::forward<decltype(std::get<0>(t))>(std::get<0>(t)), 
                               std::forward<MARQOV::Config>(std::get<1>(t)), mtx, std::get<2>(t));
        }
        
        template <typename T>
        static MarqovType* creator(std::mutex& mtx, T& t)
        {
            return ptr_from_tuple<MarqovType>(std::forward<decltype(std::get<0>(t))>(std::get<0>(t)), 
                                              std::forward<MARQOV::Config>(std::get<1>(t)), std::forward<decltype(mtx)>(mtx), std::get<2>(t));
        }
    };
    
    template <class Sim>
    class Scheduler
    {
    private:
    public:
        /** This gives us the parameters of a simulation and we are responsible for setting everything up.
         * It has a template parameter, but of course all used parameters have to resolve to the same underlying MarqovType.
         * @param p The full set of parameters that are relevant for your Problem
         * @param filter A filter that can be applied before the actual creation of MARQOV
         */
        template <typename ParamType, typename Callable>
        void createSimfromParameter(ParamType& p, Callable filter)
        {
            auto t = filter(p);
            auto simptr = sims_helper2<typename Sim::HamiltonianType, typename Sim::Lattice, ParamType>::template creator(mutexes.hdf, t);
            oursims.push_back(simptr);
            this->enqueuesim(*simptr);
        }
        std::vector<Sim*> oursims; ///< Collects the sims that we have created and for which we feel repsonsible.
        /** This registers an already allocated simulation with us.
         */
        void enqueuesim(Sim& sim)
        {
            int idx = simvector.size();
            simvectormutex.lock();
            simvector.push_back(&sim);//NOTE: We only take care about sims that go through addsims.
            simvectormutex.unlock();
            taskqueue.enqueue([&, idx]{
                //work here
                simvectormutex.lock();
                auto mysim = simvector[idx];
                simvectormutex.unlock();
                mysim->init();
                mysim->wrmploop();
                //enqueue the next full work item into the workqueue immediately
                workqueue.push_back(Simstate(idx));
            });//Put some warmup into the taskqueue
        }
        void start()
        {
            //create dummy data for the ptplan
            for (int i = 0; i < maxpt; ++i)
                ptplan.emplace_back(-1, -1 );
            
            //         std::cout<<"Starting up master"<<std::endl;
            Simstate itm;
            
            while(!masterstop)
            {
                //              std::cout<<"Master waiting for work"<<std::endl;
                bool busy = false;
                //The following wait_for construct hides a bug that occurs if the last notify in the gameloop triggers the master, but the associated task is still running.
                masterwork.wait_for(std::chrono::seconds(10), [&]{
                    busy = workqueue.pop_front(itm);
                    if (!busy)
                        masterstop = nowork();
                    return busy || masterstop;
                });
                if(busy) //there really is sth. to do
                {
                    //                  std::cout<<"dealing with work"<<std::endl;
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
        Scheduler(int maxptsteps) : maxpt(maxptsteps), masterstop(false), masterwork{},
        workqueue(masterwork),
        taskqueue(std::thread::hardware_concurrency())
        {}
        ~Scheduler() {
            if (!nowork() && !masterstop && (taskqueue.tasks_enqueued() > 0) )
            {
                masterwork.wait([&]{
                    return workqueue.is_empty() && (taskqueue.tasks_enqueued() == 0) && (taskqueue.tasks_assigned() == 0);
                });
            }
            masterstop = true;
            for (auto sim : oursims)
                delete sim;
            //         std::cout<<"Deleting Scheduler"<<std::endl;
        }
    private:
        struct Simstate
        {
            Simstate() : id(-1), npt(-100) {}
            Simstate(int i) : id(i), npt(0) {}
            Simstate(int i, int np) : id(i), npt(np) {}
            int id;
            int npt;
        };
        struct GlobalMutexes
        {
            std::mutex hdf;//lock for the HDF5 I/O since the library for C++ is not thread-safe
            std::mutex io;// lock for the rest?
        } mutexes;
        auto findpartner(uint id)
        {
            return std::find_if(ptqueue.cbegin(), ptqueue.cend(), [&id](const Simstate& itm){return itm.id == id;});
        }
        
        /** Test whether there is work available.
         * @return true if no task is working and no work is to be executed by a task and no sim is to moved to the taskqueue.
         */
        bool nowork() {return workqueue.is_empty() && taskqueue.tasks_assigned() == 0 && taskqueue.tasks_enqueued() == 0;}
        
        /** This function is called when the current simulation is up for a parallel tempering (PT) step.
         * If its partner is already waiting we do the parallel tempering, if not we got moved into 
         * a queue and wait for a partner
         * @param itm The Sim which is chosen for PT
         */
        void ptstep(Simstate itm) {
            //         std::cout<<"Parallel Tempering!"<<std::endl;
            //         std::cout<<"itm.id "<<itm.id<<" itm.npt "<<itm.npt<<std::endl;
            //         std::cout<<"Expected pairing for this time step: "<<ptplan[itm.npt].first<<" "<<ptplan[itm.npt].second<<std::endl;
            
            int partner = ptplan[itm.npt].first;
            if (partner == itm.id) partner = ptplan[itm.npt].second;//it must be the other. no exchanges with myself
            auto partnerinfo = findpartner(partner);
            
            if (partnerinfo != ptqueue.cend())
            {// partner is at the same stage, hence we can PT exchange
                //             std::cout<<"Partner found in queue"<<std::endl;
                ptqueue.erase(partnerinfo);
                calcprob();
                exchange();
                //put both sims back into the taskqueue for more processing until their next PT step
                movesimtotaskqueue(itm);
                movesimtotaskqueue(Simstate(partner, itm.npt));
            }
            else
            {//we have to wait for the PT partner
                //             std::cout<<"Partner not in queue"<<std::endl;
                ptqueue.push_back(itm);
            }
        }
        /** This determines how many steps have to be done until the next PTstep
         * and moves the simulation into the taskqueue where the gameloop is executed.
         * @param itm The simulation that gets further worked on.
         */
        void movesimtotaskqueue(Simstate itm)
        {
            auto gameloop = [&](Simstate mywork, int npt)//This defines the actual workitem that a task executes
            {
                // We loop until the next PT step
                for(; mywork.npt < npt; ++mywork.npt)
                {
                    //                 std::cout<<"Gamelooping on item "<<mywork.id<<" "<<mywork.npt<<std::endl;
                    simvector[mywork.id]->gameloop();
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
            int newnpt = findnextnpt(itm.id, itm.npt);
            //         std::cout<<"Putting a new item "<<itm.id<<" with "<<itm.npt <<" until npt = "<< newnpt<<" into the taskqueue"<< std::endl;
            taskqueue.enqueue(
                [itm, newnpt, gameloop]{gameloop(itm, newnpt);}
            );
        }
        /** Determine the next PT step
         * @param idx simulation id to check
         * @param curnpt current PT time
         * @return the next PT step where this simulation is selected for PT.
         */
        uint findnextnpt(int idx, uint curnpt)
        {
            uint retval = curnpt+1;
            while ((retval < maxpt) && (ptplan[retval].first != idx) && (ptplan[retval].second != idx))
            {
                ++retval;
            }
            return retval;
        }
        
        int maxpt; ///< how many pt steps do we do
        std::vector<Simstate> ptqueue; ///< here we collect who is waiting for its PT partner
        std::vector<std::pair<int, int> > ptplan;///< who exchanges with whom in each step
        bool masterstop;
        ThreadPool::Semaphore masterwork; ///< The semaphore that triggers the master process
        ThreadPool::ThreadSafeQueue<Simstate> workqueue; ///< this is the queue where threads put their finished work and the master does PT
        std::mutex simvectormutex; ///< A mutex to protect accesses to the simvector which could be invalidated by the use of push_back
        std::vector<Sim*> simvector; ///< An array for the full state of the simulations
        ThreadPool::Queue taskqueue; ///< this is the queue where threads pull their work from
        
        //FIXME fill those functions for proper PT
        void calcprob() {}
        void exchange() {}
    };

    //Helpers to determine if the interactions are container-like
    template <class Cont, class = void, class = void>
    struct Is_Container
    {
        static constexpr bool value = false;
    };
    
    template <class Cont>
    struct Is_Container<Cont,
    MARQOV::type_sink_t<decltype(std::declval<Cont>().size())>,
    MARQOV::type_sink_t<decltype(std::declval<Cont>().operator[](std::declval<std::size_t>()))>
    >
    {
        static constexpr bool value = true;
    };
    /** A helper class to figure out the type of the scheduler
     */
    template <class Hamiltonian, class Lattice, class Parameters>
    struct GetSchedulerType
    {
//        decltype(Hamiltonian::interactions) a = "ert";
        static_assert(Is_Container<decltype(std::declval<Hamiltonian>().interactions)>::value, "[MARQOV::Scheduler] COMPILATION FAILED: interactions are not a container.");
        typedef typename sims_helper2<Hamiltonian, Lattice, Parameters >::MarqovType MarqovType;
        typedef Scheduler<MarqovType> MarqovScheduler;
    };
    
};
#endif
