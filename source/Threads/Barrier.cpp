#include <Utils/Threads/Barrier.h>

#include <stdexcept>

fem::Barrier::Barrier(const size_t expected) noexcept
    : m_total(expected)
    , m_arrived(0)
    , m_generation(0)
{
}

void fem::Barrier::arrive_and_wait() noexcept {
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
