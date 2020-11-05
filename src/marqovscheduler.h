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
                std::cout<<"Warmuplooping on item "<<idx<<" "<<simvector.at(idx)<<" "<<simvector.size()<<std::endl;
                //work here
                simvectormutex.lock();
                auto mysim = simvector[idx];
                simvectormutex.unlock();
                mysim->init();
                mysim->wrmploop();
                //enqueue the next full work item into the workqueue immediately
                workqueue.push_back(Simstate(idx));
                std::cout<<"Finished Warmuplooping"<<std::endl;
    });//Put some warmup into the taskqueue
 }
 void start()
 {
    //create dummy data for the ptplan
     for (int i = 0; i < 2*maxpt; i += 2)
         ptplan.emplace_back(i%maxpt, (i+1)%maxpt);
     auto master = [&] /*the master thread is a lambda function since by that it captures the variables of the Scheduler*/
     {
         std::cout<<"Starting up master"<<std::endl;
            auto gameloop = [&](Simstate mywork)
            {
                std::cout<<"Gamelooping on item"<<std::endl;
                // work
                simvector[mywork.id]->gameloop();
                mywork.npt = mywork.npt + 1;
                if (mywork.npt < maxpt) // determine whether this itm needs more work
                {
                    std::cout<<"putting item again into workloop"<<std::endl;
                    workqueue.push_back(mywork);
                }
                else
                {
                    std::cout<<"no more work required"<<std::endl;
                    masterwork.notify_all();//trigger those waiting for signals from the taskqueue. since we don't push_back anything they would not be notified.
                }
            };
            Simstate itm;
            std::cout<<"Master 1 "<<masterstop<<std::endl;
            while(!masterstop)
            {
                std::cout<<"Master waiting for work"<<std::endl;
                bool busy = false;
                masterwork.wait([&]{
                 busy = workqueue.pop_front(itm);
                 return busy || masterstop;
                });
                if(busy) //there really is sth. to do
                {
                    std::cout<<"dealing with work"<<std::endl;
//                     if(ptplan[itm.npt].first == itm.id || ptplan[itm.npt].second == itm.id) // check if this sim is selected for PT
//                     {
//                         std::cout<<"Parallel Tempering!"<<std::endl;
//                         std::cout<<"itm.id "<<itm.id<<" itm.npt "<<itm.npt<<std::endl;
//                         std::cout<<ptplan[itm.npt].first<<" "<<ptplan[itm.npt].second<<std::endl;
//                         ptstep();
//                     }
//                     else
//                     {
                        std::cout<<"Putting a new item "<<itm.id<<" with"<<itm.npt <<" into the taskqueue"<< std::endl;
                        taskqueue.enqueue(
                            [itm, gameloop]{gameloop(itm);}
                        );
   //                 }
                }
            }
            std::cout<<"Master stopped"<<std::endl;
    };
     taskqueue.enqueue(master);
 }
 void waitforall() {}
 Scheduler(int maxptsteps) : maxpt(maxptsteps), masterstop(false), masterwork{},
  workqueue(masterwork),
 taskqueue(std::thread::hardware_concurrency() + 1)/*a space for the master thread*/
 {}
 ~Scheduler() {
     std::cout<<"Entering dtor of sched"<<std::endl;
     masterwork.wait([&]{
         std::cout<<"Dtor check: "<<taskqueue.tasks_enqueued()<<" "<<taskqueue.tasks_assigned()<<" "<<(workqueue.is_empty()?42:-42)<<std::endl;
        return workqueue.is_empty() && (taskqueue.tasks_enqueued() == 0) && (taskqueue.tasks_assigned() < 2); //The master thread is also enqueued in the taskqueue, so we need to account for that
     });
     masterstop = true;
     std::cout<<"Deleting Scheduler"<<std::endl;}
private:
    struct Simstate
    {
        Simstate() : id(-1), npt(-100) {}
        Simstate(int i) : id(i), npt(0) {}
        int id;
        int npt;
    };

int maxpt; ///< how many pt steps do we do
std::vector<Simstate> ptqueue; ///< here we collect who is waiting for its PT partner
std::vector<std::pair<int, int> > ptplan;///< who exchanges with whom in each step
bool masterstop;
Semaphore masterwork; ///< the semaphore that triggers the master process
ThreadSafeQueue<Simstate> workqueue; ///< this is the queue where threads put their finished work and the master does PT
std::mutex simvectormutex; ///< A mutex t protect accesses to the simvectorwhich could be invalidated by the use of push_back
std::vector<Sim*> simvector; ///< An array for the full state of the simulations
MARQOVQueue taskqueue; ///< this is the queue where threads pull their work from
void ptstep() {}
void calcprob();
void exchange();
};

#endif
