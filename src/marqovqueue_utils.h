#ifndef MARQOVQUEUE_UTILS_H
#define MARQOVQUEUE_UTILS_H
#include <atomic>
#include <deque>
#include <mutex>
#include <condition_variable>
#include <utility>
//heavily inspired by https://gitlab.com/cantordust/threadpool

namespace ThreadPool
{
    /** Spinning lock from cppreference.com
     * Future directions for optimizing the Spin-Lock: https://rigtorp.se/spinlock/
     * C++17 has a static function where we can statically check whether atomic<bool> is also lock-free
     */
    class SpinLock
    {
    public:
        SpinLock(std::atomic_flag& f) noexcept : locked(f)
        {
            lock();
        }
        
        ~SpinLock() noexcept
        {
            unlock();
        }
        void lock() noexcept {
            while (locked.test_and_set(std::memory_order_acquire)) { ; }
        }
        void unlock() noexcept {
            locked.clear(std::memory_order_release);
        }
    private:
        std::atomic_flag& locked;
    };
    
    /** A Semaphore: This is a C++11 condition variable where all accesses are protected with mutexes.
     */
    struct Semaphore
    {
        std::mutex mtx;
        std::condition_variable cv;
        
        void wait()
        {
            std::unique_lock<std::mutex> lk(mtx);
            cv.wait(lk);
        }
        
        void wait(std::function<bool()> _check)
        {
            std::unique_lock<std::mutex> lk(mtx);
            cv.wait(lk, _check);
        }
        
        template <typename T>
        void wait_for(T dt, std::function<bool()> _check)
        {
            std::unique_lock<std::mutex> lk(mtx);
            cv.wait_for(lk, dt, _check);
        }    
        
        void notify_one()
        {
            cv.notify_one();
        }
        
        void notify_all()
        {
            cv.notify_all();
        }
    };
    
    /** A class that wraps all access to a queue
     * with locks so that it is thread safe.
     * @tparam T the type of the content of the queue.
     */
    template <class T>
    class ThreadSafeQueue
    {
    public:
    private:
        std::atomic_flag lockflag = ATOMIC_FLAG_INIT;
        Semaphore& s; // the semaphore that we use for notifications
        std::deque<T> queue;
    public:
        /** Construct the queue.
         * @param _s We take the semaphore as argument that we will use to notify waiting threads about changes to the queue.
         */
        ThreadSafeQueue(Semaphore& _s) : s(_s) {}
        /** Add something to the end of the queue
         */
        void push_back(const T& t)
        {
            SpinLock lk(lockflag);
            queue.push_back(t);
            s.notify_one();
        }
        /** Get something from the front of the queue
         */
        bool pop_front(T& t)
        {
            SpinLock lk(lockflag);
            if(queue.empty()) return false;
            t = std::move(queue.front());
            queue.pop_front();
            return true;
        }
        /** Checks whether the queue is empty.
         * @return true if the queue is empty, else false
         */
        bool is_empty()
        {
            SpinLock lk(lockflag);
            return queue.empty();
        }
        /** create something at the end of the queue
         */
        template <typename ... Args>
        void emplace(Args&& ... args)
        {
            SpinLock lk(lockflag);
            queue.emplace_back(std::forward<Args>(args)...);
            s.notify_one();
        }
        /** Clear the queue
         */
        void clear()
        {
            SpinLock lk(lockflag);
            queue.clear();
        }
    };
};
#endif
