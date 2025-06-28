#include <Threads/TimersFactory.h>
#include <Threads/ICallbackQueue.h>

mthread::TimersFactory::TimersFactory(std::string thread_name)
    : m_thread_name(std::move(thread_name))
{
}

mthread::TimersFactory::~TimersFactory() {
    {
        std::lock_guard guard(m_mutex);
        m_stopped = true;
        m_cv.notify_one();
    }

    if (m_main_thread.joinable()) {
        m_main_thread.join();
    }
}

class mthread::TimersFactory::TimeOutedCallback::Impl : public ITimersFactory::ITimerHandlerImpl {
public:
    TimersFactory& owner;
    const std::chrono::steady_clock::time_point expiration_time;
    const std::weak_ptr<ICallbackQueue> callback_queue;
    mutable std::function<void()> callback;
    mutable TimeOutedCallbacksSet::iterator iter;

    Impl(TimersFactory& papa,
        std::chrono::steady_clock::time_point expTime,
        std::weak_ptr<ICallbackQueue> cbq,
        std::function<void()> cb, TimeOutedCallbacksSet::iterator it)
        : owner(papa)
        , expiration_time(expTime)
        , callback_queue(std::move(cbq))
        , callback(std::move(cb))
        , iter(it)
    {
    }

    void stop_timer_from_handler() const override {
        owner.stop_timer_impl(this);
    }
};

mthread::TimersFactory::TimeOutedCallback::TimeOutedCallback(TimersFactory& papa,
    std::chrono::steady_clock::time_point expTime,
    std::weak_ptr<ICallbackQueue> cbq,
    std::function<void()> cb,
    TimeOutedCallbacksSet::iterator it)
    : impl(std::make_shared<Impl>(papa, expTime, cbq, std::move(cb), it))
{
}

bool mthread::TimersFactory::TimeOutedCallback::operator<(const TimeOutedCallback& another) const {
    return impl->expiration_time < another.impl->expiration_time;
}

bool mthread::TimersFactory::TimeOutedCallback::expired(std::chrono::steady_clock::time_point now) const {
    return impl->expiration_time <= now;
}

bool mthread::TimersFactory::TimeOutedCallback::before(std::chrono::steady_clock::time_point when) const {
    return impl->expiration_time < when;
}

mthread::ITimersFactory::TimerHandler mthread::TimersFactory::start_timer(std::chrono::milliseconds delay, std::weak_ptr<ICallbackQueue> cbq, std::function<void()> callback) {
    if (cbq.expired()) {
        throw std::runtime_error("Attempt to start timer on non shared cbq");
    }
    if (m_stopped || !callback) {
        return {};
    }

    if (delay < std::chrono::milliseconds::zero()) {
        throw std::runtime_error("Negative timeout");
    }

    std::scoped_lock<std::mutex> lock(m_mutex);
    if (!m_stopped) {
        if (!m_main_thread.joinable()) {
            m_main_thread = std::thread(&TimersFactory::main_loop, this);
        }
        auto inserted_it = m_timers.emplace(*this, std::chrono::steady_clock::now() + delay, cbq, std::move(callback), m_timers.end());
        inserted_it->impl->iter = inserted_it;
        m_cv.notify_one();
        return TimerHandler(inserted_it->impl);
    }
    return {};
}

void mthread::TimersFactory::stop_timer_impl(const TimeOutedCallback::Impl* hndl) {
    std::scoped_lock lock(m_mutex);
    if (hndl->iter != m_timers.end()) {
        m_timers.erase(hndl->iter);
    }
}


void mthread::TimersFactory::main_loop() {
    std::unique_lock<std::mutex> lock(m_mutex);
    while (!m_stopped) {
        const auto now = std::chrono::steady_clock::now();
        while (!m_timers.empty() && m_timers.begin()->expired(now)) {
            if (m_stopped) {
                return;
            }
            auto front = m_timers.begin();
            auto& expired = front->impl;
            if (auto cbq = expired->callback_queue.lock()) {
                cbq->add_callback(std::move(expired->callback));
            }
            expired->iter = m_timers.end();
            m_timers.erase(front);
        }
        try {
            if (m_timers.empty()) {
                m_cv.wait(lock, [this]() {
                    return !m_timers.empty() || m_stopped;
                    });
            }
            else {
                const auto until = m_timers.begin()->impl->expiration_time;
                m_cv.wait_until(lock, until,
                    [this, until]() { return m_stopped || (!m_timers.empty() && m_timers.begin()->before(until)); });
            }
        }
        catch (const std::exception& ex) {
            ex.what();            
        }
    }
}

std::shared_ptr<mthread::ITimersFactory> mthread::make_timers_factory(std::string thread_name) {
    return std::make_shared<TimersFactory>(std::move(thread_name));
}
