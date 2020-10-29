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
    class Master
    {
        bool selectedforpt(Sim sim)
        {
            return ptplan[sim.npt].first == sim.id || ptplan[sim.npt].second == sim.id
        }
        void operator()()
        {
            Simstate itm;
            while(true)
            {
                masterwork.wait([&]{
                 workqueue.pop_fron(itm);
                });
                if(selectedforpt(itm))
                {
                    // do pt
                }
                else
                {
                    taskqueue.enqueue(GameLoop());
                }
            }
        }
    };
public:
 void enqueuesim(Sim&& sim) ///< the entry point for the user. Currently it's undecided whether the sim is instantiated by the user or by the scheduler
 {
     
 }
 void start()
 {
     taskqueue.enqueue(Master());
 }
 void waitforall();
 Scheduler(int maxptsteps) : taskqueue(std::thread::hardware_concurrency() + 1),/*a space for the master thread*/
 masterwork{},
 workqueue(masterwork),
 maxpt(maxptsteps)
 {
     //create dummy data for the ptplan
     for (int i = 0; i < 2*maxpt; i += 2)
         ptplan.emplace_back(i%maxpt, (i+1)%maxpt);
}
 ~Scheduler() {waitforall();}
private:
    struct Simstate
    {
        int id;
        int npt;
    };
    
    class GameLoop
    {
        void operator()(){}
    };
    class WarmUpLoop
    {
        void operator()(){}
    };
    
    
MARQOVQueue taskqueue; ///< this is the queue where threads pull their work from
Semaphore masterwork; ///< the semaphore that triggers the master process
ThreadSafeQueue<Simstate> workqueue; ///< this is the queue where threads put their finished work and the master does PT
std::vector<Sim> simvector; ///< An array for the full state of the simulations
std::vector<int> ptqueue; ///< here we collect who is waiting for its PT partner
std::vector<std::pair<int, int> > ptplan;///< who exchanges with whom in each step
int maxpt; ///< how many pt steps do we do
void ptstep();
void calcprob();
void exchange();
};

#endif
