#pragma once

template<class T>
fem::BlockingQueue<T>::BlockingQueue(size_t max_size, OnOverflow on_overflow, OnBecomeNormal on_become_normal)
    : m_max_size(max_size)
    , m_on_overflow(std::move(on_overflow))
    , m_on_become_normal(std::move(on_become_normal))
{
}

template<class T>
void fem::BlockingQueue<T>::pop(T& item) {
    pop_impl(item, [this](std::unique_lock<std::mutex>& lock) {
        m_cv.wait(lock, [this]() { return !m_queue.empty(); });
        return true;
        }
    );
}

template<class T>
bool fem::BlockingQueue<T>::try_pop(T& item) {
    return pop_impl(item, [this](auto& /*lock*/) { return !m_queue.empty(); });
}

template<class T>
template <typename Duration>
bool fem::BlockingQueue<T>::try_pop(T& item, const Duration duration) {
    return pop_impl(item, [this, duration](std::unique_lock<std::mutex>& lock) {
        return !m_cv.wait_for(lock, duration, [this]() { return !m_queue.empty(); });
        }
    );
}

template<class T>
void fem::BlockingQueue<T>::push(T item) {
    std::unique_lock mlock(m_mutex);
    const auto old_size = m_queue.size();
    const auto max_size = m_max_size;

    m_queue.push(std::move(item));
    m_cv.notify_one();
    mlock.unlock();

    if (old_size == max_size && m_on_overflow) {
        m_on_overflow(old_size + 1);
    }
}

template<class T>
bool fem::BlockingQueue<T>::try_push(T item) {
    std::scoped_lock mlock(m_mutex);
    if (m_queue.size() < m_max_size) {
        m_queue.push(std::move(item));
        m_cv.notify_one();
        return true;
    }
    return false;
}

template<class T>
size_t fem::BlockingQueue<T>::size() const {
    std::scoped_lock lock(m_mutex);
    return m_queue.size();
}

template<class T>
void fem::BlockingQueue<T>::clear() {
    std::queue<T> q;
    std::scoped_lock lock(m_mutex);
    m_queue.swap(q);
}

template<class T>
void fem::BlockingQueue<T>::set_max_size(size_t size) {
    std::scoped_lock lock(m_mutex);
    m_max_size = size;
}
