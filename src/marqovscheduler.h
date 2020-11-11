#ifndef MARQOVSCHEDULER_H
#define MARQOVSCHEDULER_H
/*
MIT License

Copyright (c) 2020 Florian Goth
fgoth@physik.uni-wuerzburg.de

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <vector>
#include <string>
#include <mutex>
#include <algorithm>
#include <chrono>

#include "marqovqueue.h"

template <class Sim>
class Scheduler
{
private:
public:
    void enqueuesim(Sim& sim) ///< the entry point for the user. Currently it's undecided whether the sim is instantiated by the user or by the scheduler
    {
        int idx = simvector.size();
        simvectormutex.lock();
        simvector.push_back(&sim);//FIXME: Currently I don't know how a sim terminates...
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
    uint findnextnpt(int idx, uint curnpt)
    {
        uint retval = curnpt+1;
        while ((ptplan[retval].first != idx) && (ptplan[retval].second != idx) && (retval < maxpt))
        {
            ++retval;
        }
        return retval;
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
//             std::cout<<"Master waiting for work"<<std::endl;
            bool busy = false;
            //The following wait_for construct hides a bug that occurs if the last notify in the gameloop triggers the master, but the associated task is still running.
            masterwork.wait_for(std::chrono::seconds(10), [&]{
                busy = workqueue.pop_front(itm);
                masterstop = nowork();
                return busy || masterstop;
            });
            if(busy) //there really is sth. to do
            {
//                 std::cout<<"dealing with work"<<std::endl;
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
//         std::cout<<"Master stopped"<<std::endl;
        //};
        //      taskqueue.enqueue(master);
    }
    bool nowork() {return workqueue.is_empty() && taskqueue.tasks_assigned() == 0 && taskqueue.tasks_enqueued() == 0;}
    void waitforall() {}
    Scheduler(int maxptsteps) : maxpt(maxptsteps), masterstop(false), masterwork{},
    workqueue(masterwork),
    taskqueue(std::thread::hardware_concurrency())
    {}
    ~Scheduler() {
//         std::cout<<"Entering dtor of sched"<<std::endl;
        if (!nowork() && !masterstop)
        {
            masterwork.wait([&]{
                return workqueue.is_empty() && (taskqueue.tasks_enqueued() == 0) && (taskqueue.tasks_assigned() == 0);
            });
        }
            masterstop = true;
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
    auto findpartner(uint id)
    {
        return std::find_if(ptqueue.cbegin(), ptqueue.cend(), [&id](const Simstate& itm){return itm.id == id;});
    }
    
    int maxpt; ///< how many pt steps do we do
    std::vector<Simstate> ptqueue; ///< here we collect who is waiting for its PT partner
    std::vector<std::pair<int, int> > ptplan;///< who exchanges with whom in each step
    bool masterstop;
    Semaphore masterwork; ///< The semaphore that triggers the master process
    ThreadSafeQueue<Simstate> workqueue; ///< this is the queue where threads put their finished work and the master does PT
    std::mutex simvectormutex; ///< A mutex to protect accesses to the simvector which could be invalidated by the use of push_back
    std::vector<Sim*> simvector; ///< An array for the full state of the simulations
    MARQOVQueue taskqueue; ///< this is the queue where threads pull their work from
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
    //FIXME fill those functions for proper PT
    void calcprob() {}
    void exchange() {}
};

#endif
