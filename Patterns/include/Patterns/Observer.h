#pragma once

#include <mutex>
#include <vector>

namespace pattern {

    template <typename... Args>
    class Observer {
    public:
        virtual ~Observer() = default;
        virtual void update(Args... args) = 0;
    };

    template <typename... Args>
    class Observable {
    public:
        void add_observer(Observer<Args...>* observer) {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_observers.push_back(observer);
        }

        void remove_observer(Observer<Args...>* observer) {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_observers.erase(
                std::remove(m_observers.begin(), m_observers.end(), observer),
                m_observers.end()
            );
        }

        void notify(Args... args) {
            std::vector<Observer<Args...>*> observers_copy;
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                observers_copy = m_observers;
            }
            for (auto* observer : observers_copy) {
                observer->update(args...);
            }
        }

    private:
        std::vector<Observer<Args...>*> m_observers;
        mutable std::mutex m_mutex;
    };

} // namespace pattern
