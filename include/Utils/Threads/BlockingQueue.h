#pragma once

#include <Utils/Patterns/MoveOnlyBase.h>

#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>

namespace util::mthrd {

    template <typename T>
    class BlockingQueue : public pattern::MoveOnlyBase {
    public:
        using OnOverflow = std::function<void(size_t)>;
        using OnBecomeNormal = std::function<void(size_t)>;

        BlockingQueue(size_t max_size, OnOverflow on_overflow = OnOverflow(), OnBecomeNormal on_become_normal = OnBecomeNormal());

        template <typename Duration>
        bool try_pop(T& item, const Duration duration);
        bool try_pop(T& item);
        void pop(T& item);

        bool try_push(T item);
        void push(T item);        

        void clear();

        size_t size() const;
        void set_max_size(size_t size);

    private:
        bool pop_impl(T& item, std::function<bool(std::unique_lock<std::mutex>&)> check_wait) {
            if (std::unique_lock<std::mutex> lock(m_mutex); check_wait(lock)) {
                item = std::move(m_queue.front());
                m_queue.pop();
                const auto new_size = m_queue.size();
                const auto max_size = m_max_size;
                lock.unlock();

                if (new_size == max_size - 1 && m_on_become_normal) {
                    m_on_become_normal(new_size);
                }
                return true;
            }
            return false;
        }

    private:
        std::queue<T> m_queue;
        mutable std::mutex m_mutex;
        std::condition_variable_any m_cv;

        size_t m_max_size;
        const OnOverflow m_on_overflow;
        const OnBecomeNormal m_on_become_normal;
    };

} // namespace util::mthrd

#include "BlockingQueue.tpp"
