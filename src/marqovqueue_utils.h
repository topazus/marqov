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
    /** Spinning lock from cppreference.com .
     * 
     * @note Future directions for optimizing the Spin-Lock: https://rigtorp.se/spinlock/
     * C++17 has a static function where we can statically check whether atomic<bool> is also lock-free
     */
    class SpinLock
    {
    public:
        /** Set up the spinning lock.
         * 
         * @param f The flag on which we are spinning.
         */
        SpinLock(std::atomic_flag& f) noexcept : locked(f)
        {
            lock();
        }
        
        /** Unlock and destroy the lock.
         */
        ~SpinLock() noexcept
        {
            unlock();
        }
        /** Spin until WE can lock the flag.
         */
        void lock() noexcept {
            while (locked.test_and_set(std::memory_order_acquire)) { ; }
        }
        /** Unlock the spinning lock.
         */
        void unlock() noexcept {
            locked.clear(std::memory_order_release);
        }
    private:
        std::atomic_flag& locked;///< the flag that we use for locking. 
    };
    
    /** A Semaphore: This is a C++11 condition variable where all accesses are protected with mutexes.
     */
    struct Semaphore
    {
        std::mutex mtx; ///< The mutex that we internally use for setting up a lock.
        std::condition_variable cv; ///< A condition variable for notifxing waiting threads.
        /** Blocks the current thread until the condition variable is woken up.
         */
        void wait()
        {
            std::unique_lock<std::mutex> lk(mtx);
            cv.wait(lk);
        }
        
        /** Wait until the condition variable is woken up and check the predicate.
         * 
         * @param _check A predicate that is to be checked.
         */
        void wait(std::function<bool()> _check)
        {
            std::unique_lock<std::mutex> lk(mtx);
            cv.wait(lk, _check);
        }
        
        /** Wait until the condition variable is woken up, check the predicate, or just wait for a certain amount of time.
         * 
         * @tparam T The type of the time difference.
         * 
         * @param dt the length of time that we wait.
         * @param _check A predicate that is to be checked.
         */
        template <typename T>
        void wait_for(T dt, std::function<bool()> _check)
        {
            std::unique_lock<std::mutex> lk(mtx);
            cv.wait_for(lk, dt, _check);
        }    
        
        /** Notify one thread that is waiting for the condition variable.
         */
        void notify_one()
        {
            cv.notify_one();
        }
        
        /** Notify all threads that are waiting for the condition variable.
         */
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
        std::atomic_flag lockflag = ATOMIC_FLAG_INIT;///< This flag is used for protecting access to the queue.
        Semaphore& s; ///< the semaphore that we use for notifications
        std::deque<T> queue; ///< This internal queue serves as the storage for the items in the queue. Access is protected by a lock flag.
    public:
        /** Construct the queue.
         * 
         * @param _s We take the semaphore as argument that we will use to notify waiting threads about changes to the queue.
         */
        ThreadSafeQueue(Semaphore& _s) : s(_s) {}
        /** Add something to the end of the queue.
         * 
         * @param t the new item that we put into the queue.
         */
        void push_back(const T& t)
        {
            SpinLock lk(lockflag);
            queue.push_back(t);
            s.notify_one();
        }
        /** Get something from the front of the queue.
         * 
         * @param t the location into which we store the received content.
         * @return false if the queue is empty, else true.
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
         * 
         * @return true if the queue is empty, else false.
         */
        bool is_empty()
        {
            SpinLock lk(lockflag);
            return queue.empty();
        }
        /** Create something at the end of the queue.
         * 
         * @tparam Args the parameter pack of arguments
         * 
         * @param args the arguments of the thing we want to construct.
         */
        template <typename ... Args>
        void emplace(Args&& ... args)
        {
            SpinLock lk(lockflag);
            queue.emplace_back(std::forward<Args>(args)...);
            s.notify_one();
        }
        /** Clear the queue.
         */
        void clear()
        {
            SpinLock lk(lockflag);
            queue.clear();
        }
    };
};
#endif
