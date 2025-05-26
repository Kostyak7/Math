#pragma once

#include <functional>
#include <memory>
#include <chrono>
#include <string>

namespace fem {

    class ICallbackQueue;

    class ITimersFactory {
    protected:
        class ITimerHandlerImpl {
        public:
            virtual void stop_timer_from_handler() const = 0;
            virtual ~ITimerHandlerImpl() = default;
        };

    public:
        class TimerHandler {
        public:
            TimerHandler() = default;
            TimerHandler(std::weak_ptr<const ITimerHandlerImpl> impl);
            void stop_timer();

        private:
            std::weak_ptr<const ITimerHandlerImpl> m_impl;
        };

    public:
        virtual TimerHandler start_timer(std::chrono::milliseconds delay, std::weak_ptr<ICallbackQueue> cbq, std::function<void()> callback) = 0;
        virtual ~ITimersFactory() = default;
    };

    std::shared_ptr<ITimersFactory> make_timers_factory(std::string thread_name);

} // namespace fem
