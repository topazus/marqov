#ifndef MARQOVQUEUE_H
#define MARQOVQUEUE_H
#include <memory>
#include <functional>
#include <future>
#include <utility>
#include <atomic>
#include <mutex>
#include <unordered_map>

#include "marqovqueue_helpers.h"

namespace ThreadPool
{
    /** Task Queue.
     * 
     * A queue of the tasks that are about to be done and then get distributed
     * to workers.
     */
    class Queue
    {
    public:
        using Task = std::function<void()>;///< The type of the functions that can be enqueued in our task pool. We do not support function arguments.
        /** Construct the work queue.
         * 
         * @param initc number of threads. Defaults to the number of hardware threads.
         */
        Queue(const uint initc = std::thread::hardware_concurrency()) : semaphores{}, queue(semaphores.work), workers(initc)
        {
            resize(initc);
        }
        /** Move constructor of queue.
         *
         */
        Queue(Queue&& rhs) : semaphores{}, queue(semaphores.work), workers{std::move(rhs.workers)} {}
        ~Queue()
        {
            flags.stop.store(true);
            sync();
            while (workers.count); // spin until all workers have exited
        }
        /** Resize the pool of threads.
         * 
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
        /** The entry point function where we register our tasks.
         * 
         * We do not support arguments.
         * @tparam F A placeholder for the function.
         * 
         * @param f The task.
         * @return A future to obtain a possible result.
         */
        template <typename F>
        auto enqueue(F&& f)
        {
            using Ret = typename std::result_of<F&(void)>::type;
            // wrap f into a packaged_task structure and afterwards pack it into a shared_pointer
            auto task(std::make_shared<std::packaged_task<void()>>(std::forward<F>(f)));
            std::future<Ret> future(task->get_future()); //connect a future to the task so that we may obtain a result
            if (flags.stop) return future;
            queue.emplace([=]{(*task)();});//put the work into the queue
            ++stats.enqueued;
            return future;
        }
        /** Obtain the number of currently enqueued tasks.
         * 
         * @return The number of currently enqueued tasks.
         */
        uint tasks_enqueued() const
        {
            return stats.enqueued;
        }
        /** Obtain the number of currently assigned tasks.
         * 
         * @return The number tasks currently being processed.
         */
        uint tasks_assigned() const
        {
            return stats.assigned;
        }
    private:
        using auint = std::atomic<uint>; ///< The type of a counter with atomic access.
        using Toggle = std::atomic<bool>; ///< A short cut for using atomic bools.
        /** Semaphores that trigger behaviour of threads.
         */
        struct Semaphores
        {
            Semaphore work;///< This signals that the threads may begin processing.
            Semaphore sync;///< This signals that the threads should finish.
        } semaphores;///< The two semaphores that control the behaviour of our threads.
        ThreadSafeQueue<Task> queue;///< The actual task queue.
        /** A helper structure to encapsulate some useful flags.
         */
        struct Flags
        {
            Toggle stop;///< A flag for denoting that the threads should stop.
            Toggle prune;///< A flag for denoting that the threads should terminate.
            Flags() noexcept : stop(false), prune(false) {}
            Flags(const Flags& rhs) noexcept : stop(rhs.stop.load()), prune(rhs.prune.load()) {}
        } flags;///< flags
        
        /** A helper structure for some threadpool statistics.
         */
        struct Stats
        {
            auint enqueued; ///< count how many threads are in the queue.
            auint assigned; ///< count how many threads are worked on.
            /** Construct our statistics with initial zero values.
             */
            Stats() noexcept : enqueued(0), assigned(0) {}
        } stats;///< Actual instance for our statistics.

        /** A helper structure to bundle the bookkeeping of workers.
         */
        struct Workers
        {
            std::mutex busy_mtx;///< This mutex protects access to the below map.
            std::unordered_map<std::thread::id, Toggle> busy; ///< This provides a map between the thread and whether it is busy.
            auint count; ///< current count of workers
            auint target_count; ///< intended size of worker threads in the pool.
            /** Create the worker bookkeeping.
             * 
             * @param tc How many workers do we intend to have.
             */
            Workers(const uint tc) noexcept : busy_mtx{}, busy{}, count(0), target_count(tc) {}
            /** Move Constructor
             *
             * Takes ownership of other threads.
             */
            Workers(Workers&& rhs) : busy_mtx{}, busy{}, count(rhs.count.load()), target_count(rhs.target_count.load())
            {
        	std::lock_guard<std::mutex> lk(rhs.busy_mtx);
        	std::swap(busy, rhs.busy);
        	rhs.count.store(0);
            }
        } workers; ///< An instance of the worker bookkeeping structure.

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
            
            // let's create the actual thread and immediately detach it from the main thread
            std::thread(
                [&, rp=std::move(ready_promise)]() mutable
                {
                    Task task;//This will be the place holder of the actual work item
                    // let's do the set up of the book keeping structures
                    ++workers.count;
                    //                std::cout<<"\tWorker "<< workers.count<< " in thread "<< std::this_thread::get_id()<< " ready"<<std::endl;
                    {
                        std::lock_guard<std::mutex> lk(workers.busy_mtx);
                        workers.busy.emplace(std::this_thread::get_id(), true);
                    }
                    auto& busy(workers.busy[std::this_thread::get_id()]);
                    rp.set_value(std::this_thread::get_id());
                    busy.store(false);
                    while(true)// We will basically loop forever in this loop
                    {
                        // Wait until we have something to do
                        semaphores.work.wait([&]
                        {
                            busy.store(queue.pop_front(task));
                            //                        std::cout<<"Thread "<<std::this_thread::get_id()<<" waiting for work."<<std::endl;
                            return (busy || flags.stop || flags.prune);
                        });
                        // We got something to do
                        if (busy)
                        {
                            ++stats.assigned;
                            --stats.enqueued;
                            
                            task();// execute the work item
                            busy.store(false);
                            
                            --stats.assigned;
                            
                            if(!stats.enqueued)
                            {
                                semaphores.sync.notify_all();
                            }
                        }
                        else if (flags.prune || flags.stop){break;}
                    }
                    // We were told to stop -> tidy up
                    {
                        std::lock_guard<std::mutex> lk(workers.busy_mtx);
                        workers.busy.erase(std::this_thread::get_id());
                    }
                    --workers.count;
                    flags.prune.store(workers.count > workers.target_count);
                    //                std::cout<<"\tWorker "<<count<<" in thread "<< std::this_thread::get_id()<<" exiting..."<<std::endl;
                }).detach();
                return id.get();
        }
    };
};
#endif
