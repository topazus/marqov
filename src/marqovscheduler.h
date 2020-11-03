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

#include "marqovqueue.h"

template <class Sim>
class Scheduler
{
private:
public:
 void enqueuesim(Sim&& sim) ///< the entry point for the user. Currently it's undecided whether the sim is instantiated by the user or by the scheduler
 {
     auto warmuploop = [&](int id)
            {
                std::cout<<"Warmuplooping on item"<<std::endl;
                //work here
                simvector[id].wrmploop();
                //enqueue the next full work item into the workqueue immediately
                workqueue.push_back(Simstate(id));
            };
     int idx = simvector.size();
     simvector.push_back(sim);//FIXME: Currently I don't know how a sim terminates...
     taskqueue.enqueue([idx, warmuploop]{warmuploop(idx);});//Put some warmup into the taskqueue
 }
 void start()
 {
     auto master = [&] /*the master thread is a lambda function since by that it captures the variables of the Scheduler*/
     {
            auto gameloop = [&](Simstate mywork)
            {
                std::cout<<"Gamelooping on item"<<std::endl;
                // work
                simvector[mywork.id].gameloop();
                mywork.npt = mywork.npt + 1;
                if (mywork.npt < maxpt) // determine whether this itm needs more work
                    workqueue.push_back(mywork);
            };
            Simstate itm;
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
                    if(ptplan[itm.npt].first == itm.id || ptplan[itm.npt].second == itm.id) // check if this sim is selected for PT
                        ptstep();
                    else
                    {
                        std::cout<<"Putting a new item "<<itm.id<<" with"<<itm.npt <<" into the taskqueue"<<std::endl;
                        taskqueue.enqueue([&]{gameloop(itm);});
                    }
                }
            }
    };
     taskqueue.enqueue(master);
 }
 void waitforall() {}
 Scheduler(int maxptsteps) : taskqueue(std::thread::hardware_concurrency() + 1),/*a space for the master thread*/
 masterwork{},
 workqueue(masterwork),
 maxpt(maxptsteps),
 masterstop(false)
 {
     //create dummy data for the ptplan
     for (int i = 0; i < 2*maxpt; i += 2)
         ptplan.emplace_back(i%maxpt, (i+1)%maxpt);
}
 ~Scheduler() {masterstop = true;}
private:
    struct Simstate
    {
        Simstate() : id(-1), npt(-100) {}
        Simstate(int i) : id(i), npt(0) {}
        int id;
        int npt;
    };
MARQOVQueue taskqueue; ///< this is the queue where threads pull their work from
Semaphore masterwork; ///< the semaphore that triggers the master process
ThreadSafeQueue<Simstate> workqueue; ///< this is the queue where threads put their finished work and the master does PT
std::vector<Sim> simvector; ///< An array for the full state of the simulations
std::vector<Simstate> ptqueue; ///< here we collect who is waiting for its PT partner
std::vector<std::pair<int, int> > ptplan;///< who exchanges with whom in each step
int maxpt; ///< how many pt steps do we do
bool masterstop;
void ptstep() {}
void calcprob();
void exchange();
};

#endif
