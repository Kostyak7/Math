#include <Threads/Barrier.h>

#include <stdexcept>

mthread::Barrier::Barrier(size_t expected) noexcept
    : m_total(expected)
    , m_arrived(0)
    , m_generation(0)
    , m_broken(false)
{
}

void mthread::Barrier::arrive_and_wait() {
    std::unique_lock<std::mutex> lock(m_mutex);
    const auto gen = m_generation.load();  // Текущий "поколение" 

    if (++m_arrived == m_total) {
        // Все потоки прибыли — пробуждаем все и начинаем новое поколение
        m_generation.fetch_add(1);
        m_arrived.store(0);
        m_cv.notify_all();
        return;
    }
    m_cv.wait(lock, [this, gen] { return gen != m_generation.load(); });
}

void mthread::Barrier::arrive_and_drop() {
    std::unique_lock<std::mutex> lock(m_mutex);
    if (m_broken) {
        throw std::runtime_error("Barrier is broken");
    }

    if (m_total == 0) {
        throw std::runtime_error("No threads expected in barrier");
    }

    --m_total;
    if (++m_arrived == m_total) {
        m_arrived = 0;
        ++m_generation;
        m_cv.notify_all();
    }
}

void mthread::Barrier::arrive() {
    std::unique_lock<std::mutex> lock(m_mutex);
    if (m_broken) {
        throw std::runtime_error("Barrier is broken");
    }

    if (++m_arrived == m_total) {
        m_arrived = 0;
        ++m_generation;
        m_cv.notify_all();
    }
}

void mthread::Barrier::wait() {
    std::unique_lock<std::mutex> lock(m_mutex);
    if (m_broken) {
        throw std::runtime_error("Barrier is broken");
    }

    const size_t gen = m_generation.load();
    if (m_arrived == m_total) {
        return;  // Барьер уже достигнут
    }

    m_cv.wait(lock, [this, gen] {
        return gen != m_generation.load() || m_broken;
        });
    if (m_broken) {
        throw std::runtime_error("Barrier is broken");
    }
}

void mthread::Barrier::break_barrier() {
    std::unique_lock<std::mutex> lock(m_mutex);
    m_broken = true;
    m_cv.notify_all();
}

size_t mthread::Barrier::get_expected() const {
    return m_total;
}

size_t mthread::Barrier::get_arrived() const {
    return m_arrived.load();
}