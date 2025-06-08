#include <Utils/Threads/CallbackQueue.h>

namespace {
    constexpr size_t DEFAULT_CBQ_MAX_SIZE = 1000;
} // namespace 

util::mthrd::CallbackQueue::CallbackQueue(std::string name)
    : CallbackQueue(std::move(name), DEFAULT_CBQ_MAX_SIZE)
{
}

util::mthrd::CallbackQueue::CallbackQueue(std::shared_ptr<ITimersFactory> timers_factory, std::string name)
    : CallbackQueue(std::move(timers_factory), std::move(name), DEFAULT_CBQ_MAX_SIZE)
{
}

util::mthrd::CallbackQueue::CallbackQueue(std::string name, size_t max_size, OnOverflow on_overflow_callback, OnBecomeNormal on_become_normal_callback)
    : CallbackQueue(make_timers_factory("TF" + name), std::move(name), max_size, std::move(on_overflow_callback), std::move(on_become_normal_callback))
{
}
   

util::mthrd::CallbackQueue::CallbackQueue(std::shared_ptr<ITimersFactory> timers_factory, std::string name, size_t max_size, OnOverflow on_overflow_callback, OnBecomeNormal on_become_normal_callback)
    : m_name(std::move(name))
    , m_queue(max_size, std::move(on_overflow_callback), std::move(on_become_normal_callback))
    , m_timers_factory(std::move(timers_factory))
{
    start();
}

util::mthrd::CallbackQueue::~CallbackQueue() {
    if (m_queue_thread.get_id() != std::this_thread::get_id()) {
        throw std::runtime_error{ "Try dstory callback queue inside itself: " + get_name() };
    }
    m_timers_factory.reset();
    stop();
}

void util::mthrd::CallbackQueue::start() {
    m_queue_thread = std::thread([this]() mutable {
        main_loop();
        });
    m_therad_id.store(m_queue_thread.get_id());
}

void util::mthrd::CallbackQueue::stop() {
    m_stopped.store(true);

    if (m_queue_thread.joinable()) {
        m_queue.push(nullptr);
        m_queue_thread.join();
        m_therad_id.store(std::thread::id());
    }
}

void util::mthrd::CallbackQueue::add_callback(std::function<void()> callback) {
    if (m_stopped) {
        return;
    }

    if (!callback) {
        throw std::runtime_error("No callback");
    }

    m_queue.push(std::move(callback));
}

void util::mthrd::CallbackQueue::add_callback(std::function<void()> callback, std::weak_ptr<const void> lifetime) {
    add_callback(make_safe_callback(std::move(callback), std::move(lifetime)));
}

bool util::mthrd::CallbackQueue::try_add_callback(std::function<void()> callback) {
    if (m_stopped) {
        return false;
    }

    if (!callback) {
        throw std::runtime_error("No callback");
    }

    return m_queue.try_push(std::move(callback));
}

bool util::mthrd::CallbackQueue::try_add_callback(std::function<void()> callback, std::weak_ptr<const void> lifetime) {
    return try_add_callback(make_safe_callback(std::move(callback), std::move(lifetime)));
}

util::mthrd::ICallbackQueue::DelayedHandler util::mthrd::CallbackQueue::add_delayed_callback(std::function<void()> callback, std::chrono::milliseconds time_out) {
    if (m_stopped) {
        return {};
    }
    return m_timers_factory->start_timer(time_out, weak_from_this(), std::move(callback));
}

util::mthrd::ICallbackQueue::DelayedHandler util::mthrd::CallbackQueue::add_delayed_callback(std::function<void()> callback, std::chrono::milliseconds time_out, std::weak_ptr<const void> lifetime) {
    return add_delayed_callback(make_safe_callback(std::move(callback), std::move(lifetime)), time_out);
}

void util::mthrd::CallbackQueue::remove_delayed_callback(ICallbackQueue::DelayedHandler& hndl) {
    hndl.stop_timer();
}

void util::mthrd::CallbackQueue::wait(AwatingType awating_type) {
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

void util::mthrd::CallbackQueue::main_loop() {
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

std::thread::id util::mthrd::CallbackQueue::get_thread_id() const noexcept {
    return m_therad_id.load();
}

bool util::mthrd::CallbackQueue::is_working_thread() const noexcept {
    return get_thread_id() == std::this_thread::get_id();
}

std::string util::mthrd::CallbackQueue::get_name() const {
    return m_name;
}

size_t util::mthrd::CallbackQueue::size() const {
    return m_queue.size();
}

void util::mthrd::CallbackQueue::set_max_size(size_t size) {
    m_queue.set_max_size(size);
}

void util::mthrd::CallbackQueue::clear() {
    m_queue.clear();
}

std::shared_ptr<util::mthrd::ITimersFactory> util::mthrd::CallbackQueue::get_timers_factory() {
    return m_timers_factory;
}
