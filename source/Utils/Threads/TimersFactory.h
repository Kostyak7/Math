#pragma once

#include "ITimersFactory.h"

#include <set>
#include <thread>
#include <condition_variable>

namespace fem {

    class TimersFactory : public ITimersFactory {
    public:
        TimersFactory(std::string thread_name);
        ~TimersFactory() override;

        TimerHandler start_timer(std::chrono::milliseconds delay, std::weak_ptr<ICallbackQueue> cbq, std::function<void()> callback) override;        

    private:
        struct TimeOutedCallback {
            TimeOutedCallback(TimersFactory& papa, std::chrono::steady_clock::time_point exp_time, std::weak_ptr<ICallbackQueue> cbq, std::function<void()> cb, std::multiset<TimeOutedCallback>::iterator it);
            bool operator<(const TimeOutedCallback& other) const;
            bool expired(std::chrono::steady_clock::time_point now) const;
            bool before(std::chrono::steady_clock::time_point when) const;

            class Impl;
            std::shared_ptr<const Impl> impl;
        };

        using TimeOutedCallbacksSet = std::multiset<TimeOutedCallback>;

    private:
        void main_loop();
        void stop_timer_impl(const TimeOutedCallback::Impl* hndl);

    private:
        const std::string m_thread_name;
        TimeOutedCallbacksSet m_timers;
        std::mutex m_mutex;
        std::thread m_main_thread;
        std::condition_variable m_cv;
        std::atomic<bool> m_stopped = false;
    };

} // namespace fem
