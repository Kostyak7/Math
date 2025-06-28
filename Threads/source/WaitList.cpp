#include <Threads/WaitList.h>

void mthread::WaitList::add(std::future<void>&& future) {
    std::unique_lock lock(m_mutex);
    m_list.emplace_back(std::move(future));
}

void mthread::WaitList::wait() {
    std::unique_lock lock(m_mutex);
    for (const auto& future : m_list) {
        future.wait();
    }
    m_list.clear();
}
