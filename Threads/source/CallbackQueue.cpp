#include <Threads/CallbackQueue.h>

namespace {
    constexpr size_t DEFAULT_CBQ_MAX_SIZE = 1000;
} // namespace 

mthread::CallbackQueue::CallbackQueue(std::string name)
    : CallbackQueue(std::move(name), DEFAULT_CBQ_MAX_SIZE)
{
}

mthread::CallbackQueue::CallbackQueue(std::shared_ptr<ITimersFactory> timers_factory, std::string name)
    : CallbackQueue(std::move(timers_factory), std::move(name), DEFAULT_CBQ_MAX_SIZE)
{
}

mthread::CallbackQueue::CallbackQueue(std::string name, size_t max_size, OnOverflow on_overflow_callback, OnBecomeNormal on_become_normal_callback)
    : CallbackQueue(make_timers_factory("TF" + name), std::move(name), max_size, std::move(on_overflow_callback), std::move(on_become_normal_callback))
{
}
   

mthread::CallbackQueue::CallbackQueue(std::shared_ptr<ITimersFactory> timers_factory, std::string name, size_t max_size, OnOverflow on_overflow_callback, OnBecomeNormal on_become_normal_callback)
    : m_name(std::move(name))
    , m_queue(max_size, std::move(on_overflow_callback), std::move(on_become_normal_callback))
    , m_timers_factory(std::move(timers_factory))
{
    start();
}

mthread::CallbackQueue::~CallbackQueue() {
    if (m_queue_thread.get_id() != std::this_thread::get_id()) {
        throw std::runtime_error{ "Try dstory callback queue inside itself: " + get_name() };
    }
    m_timers_factory.reset();
    stop();
}

void mthread::CallbackQueue::start() {
    m_queue_thread = std::thread([this]() mutable {
        main_loop();
        });
    m_therad_id.store(m_queue_thread.get_id());
}

void mthread::CallbackQueue::stop() {
    m_stopped.store(true);

    if (m_queue_thread.joinable()) {
        m_queue.push(nullptr);
        m_queue_thread.join();
        m_therad_id.store(std::thread::id());
    }
}

void mthread::CallbackQueue::add_callback(std::function<void()> callback) {
    if (m_stopped) {
        return;
    }

    if (!callback) {
        throw std::runtime_error("No callback");
    }

    m_queue.push(std::move(callback));
}

void mthread::CallbackQueue::add_callback(std::function<void()> callback, std::weak_ptr<const void> lifetime) {
    add_callback(make_safe_callback(std::move(callback), std::move(lifetime)));
}

bool mthread::CallbackQueue::try_add_callback(std::function<void()> callback) {
    if (m_stopped) {
        return false;
    }

    if (!callback) {
        throw std::runtime_error("No callback");
    }

    return m_queue.try_push(std::move(callback));
}

bool mthread::CallbackQueue::try_add_callback(std::function<void()> callback, std::weak_ptr<const void> lifetime) {
    return try_add_callback(make_safe_callback(std::move(callback), std::move(lifetime)));
}

mthread::ICallbackQueue::DelayedHandler mthread::CallbackQueue::add_delayed_callback(std::function<void()> callback, std::chrono::milliseconds time_out) {
    if (m_stopped) {
        return {};
    }
    return m_timers_factory->start_timer(time_out, weak_from_this(), std::move(callback));
}

mthread::ICallbackQueue::DelayedHandler mthread::CallbackQueue::add_delayed_callback(std::function<void()> callback, std::chrono::milliseconds time_out, std::weak_ptr<const void> lifetime) {
    return add_delayed_callback(make_safe_callback(std::move(callback), std::move(lifetime)), time_out);
}

void mthread::CallbackQueue::remove_delayed_callback(ICallbackQueue::DelayedHandler& hndl) {
    hndl.stop_timer();
}

void mthread::CallbackQueue::wait(AwatingType awating_type) {
    if (m_queue_thread.get_id() == std::this_thread::get_id()) {
        throw std::logic_error("DEAD LOCK in wait() call: can't awaiting inside CallbackQueue " + get_name());
    }

    while (!m_stopped) {
        std::mutex mutex;
        std::condition_variable cv;
        bool ready = false;
        add_callback([&] {
            std::scoped_lock lock(mutex);
            ready = true;
            cv.notify_all();
            });
        std::unique_lock lock(mutex);
        cv.wait(lock, [&] { return m_stopped || ready; });

        if (awating_type == AwatingType::SCHEDULED_ONLY || m_queue.size() == 0) {
            break;
        }
    };

    if (awating_type == AwatingType::COMPREHENSIVE) {
        wait(AwatingType::SCHEDULED_ONLY);
    }
}

void mthread::CallbackQueue::main_loop() {
    std::function<void()> callback;
    do {
        callback = nullptr;
        m_queue.pop(callback); 
        if (callback) {
            try {
                callback();
            }
            catch (...) {
                // TODO std::current_exception();
            }
        }
    } while (callback);
}

std::thread::id mthread::CallbackQueue::get_thread_id() const noexcept {
    return m_therad_id.load();
}

bool mthread::CallbackQueue::is_working_thread() const noexcept {
    return get_thread_id() == std::this_thread::get_id();
}

std::string mthread::CallbackQueue::get_name() const {
    return m_name;
}

size_t mthread::CallbackQueue::size() const {
    return m_queue.size();
}

void mthread::CallbackQueue::set_max_size(size_t size) {
    m_queue.set_max_size(size);
}

void mthread::CallbackQueue::clear() {
    m_queue.clear();
}

std::shared_ptr<mthread::ITimersFactory> mthread::CallbackQueue::get_timers_factory() {
    return m_timers_factory;
}
