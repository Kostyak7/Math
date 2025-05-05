#pragma once

#include <functional>
#include <mutex>
#include <vector>

namespace pattern {

    template <typename EventT>
    class EventDispatcher {
    public:
        using Handler = std::function<void(const EventT&)>;

        void subscribe(Handler handler) {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_handlers.push_back(std::move(handler));
        }

        void dispatch(const EventT& event) {
            std::vector<Handler> handlers_copy;
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                handlers_copy = m_handlers;
            }
            for (auto& handler : handlers_copy) {
                handler(event);
            }
        }

    private:
        std::vector<Handler> m_handlers;
        mutable std::mutex m_mutex;
    };

} // namespace pattern
