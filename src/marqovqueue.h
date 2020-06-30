#ifndef MARQOVQUEUE_H
#define MARQOVQUEUE_H
#include <memory>
#include <functional>
#include <future>
#include <utility>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include "marqovqueue_utils.h"

class MARQOVQueue
{
public:
    using Task = std::function<void()>;
    /** Construct the work queue
     * @param initc number of threads. defaults to the number of hardware threads
     */
    MARQOVQueue(const uint initc = std::thread::hardware_concurrency()) : workers(initc), semaphores{}, queue(semaphores.work)
    {
        resize(initc);
    }
    ~MARQOVQueue()
    {
        flags.stop.store(true);
        sync();
        while (workers.count); // spin until all workers have exited
    }
    /** Resize the pool of threads
     * @param cnt how many threads there should be 
     */
    void resize(const uint cnt)
    {
        if(flags.stop) return;
        workers.target_count.store(cnt);
        flags.prune.store((workers.count > workers.target_count));
        while(workers.count < workers.target_count)
        {
            std::thread::id id(add_worker());
            while(workers.busy[id]);
        }
    }
    /** The entry point function where we register our tasks
     * We do not support arguments
     * @param f The task
     * @return A future to obtain a possible result
     */
    template <typename F>
    auto enqueue(F&& f)
    {
        using Ret = typename std::result_of<F&(void)>::type;
        //wrap f into a packaged_task structure and afterwards pack it into a shared_pointer
        auto task(std::make_shared<std::packaged_task<void()>>(std::forward<F>(f)));
        std::future<Ret> future(task->get_future()); //connect a future to the task so that we may obtain a result
        if (flags.stop) return future;
        queue.emplace([=]{(*task)();});//put the work into the queue
        ++stats.enqueued;
        return future;
    }
private:
	using auint = std::atomic<uint>;
    using toggle = std::atomic<bool>;
    ThreadSafeQueue<Task> queue;//the actual task queue
    struct Flags
    {
        toggle stop;
        toggle prune;
        Flags() : stop(false), prune(false) {}
    } flags;
    
    struct Semaphores//The signals
    {
        Semaphore work;// This signals that the threads my begin processing
        Semaphore sync;// This signals that the threads should finish
    } semaphores;
    
    struct Stats// A helper structure for various statistics
    {
        auint enqueued;
        Stats() : enqueued(0) {}
    } stats;
    
    struct Workers// A helper structure to bundle the bookkeeping of workers
    {
        std::mutex busy_mtx;
        std::unordered_map<std::thread::id, toggle> busy; //This provides a map between the thread and whether it is busy.
        auint count; //current count
        auint target_count; //intended size of pool
        Workers(const uint tc) : count(0), target_count(tc) {}
    } workers;
    /** trigger synchronization
     */
    void sync()
    {
        semaphores.work.notify_all();
        semaphores.sync.wait([&]
        {
            return !stats.enqueued;
        });
    }
    /** This function adds a worker
     * @return the id of the worker
     */
    std::thread::id add_worker()
    {
        std::promise<std::thread::id> ready_promise;
        std::future<std::thread::id> id(ready_promise.get_future());

        //let's create the actual thread and immediately detach it from the main thread
        std::thread(
            [&, rp=std::move(ready_promise)]() mutable
            {
                Task task;//This will be the place holder of the actual work item
                //let's do the set up of the book keeping structures
                uint count(++workers.count);
//                 std::cout<<"\tWorker "<< count<< " in thread "<< std::this_thread::get_id()<< " ready"<<std::endl;
                {
                    std::lock_guard<std::mutex> lk(workers.busy_mtx);
                    workers.busy.emplace(std::this_thread::get_id(), true);
                }
                auto& busy(workers.busy[std::this_thread::get_id()]);
                rp.set_value(std::this_thread::get_id());
                busy.store(false);
                while(true)// We will basically loop forever in this loop
                {
                    //Wait until we have something to do
                    semaphores.work.wait([&]
                    {
                        busy.store(queue.pop_front(task));
//                        std::cout<<"Thread "<<std::this_thread::get_id()<<" waiting for work."<<std::endl;
                        return (busy || flags.stop || flags.prune);
                    });
                    // We got something to do
                    if (busy)
                    {
                        --stats.enqueued;
                        task();// execute the work item
                        busy.store(false);
                        
                        if(!stats.enqueued)
                        {
                            semaphores.sync.notify_all();
                        }
                    }
                    else if (flags.prune || flags.stop){break;}
                }
                //We were told to stop -> tidy up
                {
                    std::lock_guard<std::mutex> lk(workers.busy_mtx);
                    workers.busy.erase(std::this_thread::get_id());
                }
                --workers.count;
                flags.prune.store(workers.count > workers.target_count);
//                 std::cout<<"\tWorker "<<count<<" in thread "<< std::this_thread::get_id()<<" exiting..."<<std::endl;
            }).detach();
        return id.get();
    }
};
#endif
