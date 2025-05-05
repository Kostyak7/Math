#pragma once

#include "ITimersFactory.h"

#include <future>
#include <functional>
#include <optional>
#include <thread>
#include <type_traits>

namespace fem {

    class ICallbackQueue {
    public:
        using Callback = std::function<void()>;
        using OnOverflow = std::function<void(size_t)>;
        using OnBecomeNormal = std::function<void(size_t)>;
        using DelayedHandler = ITimersFactory::TimerHandler;

        struct Params {
            size_t max_queue_size{ 1000 };
            OnOverflow on_overflow{ nullptr };
            OnBecomeNormal on_become_normal{ nullptr };
        };

    public:
        virtual ~ICallbackQueue() = default;

        enum class Priority : size_t {
            HIGH,
            MEDIUM,
            LOW
        };
        virtual void add_callback(Callback callback, Priority priority = Priority::MEDIUM) = 0;
        virtual void add_callback(Callback callback, Priority priority, std::weak_ptr<const void> lifetime) = 0;

        virtual bool try_add_callback(Callback callback, Priority priority = Priority::MEDIUM) = 0;
        virtual bool try_add_callback(Callback callback, Priority priority, std::weak_ptr<const void> lifetime) = 0;

        virtual DelayedHandler add_delayed_callback(Callback callback, std::chrono::milliseconds time_out, Priority priority = Priority::MEDIUM) = 0;
        virtual DelayedHandler add_delayed_callback(Callback callback, std::chrono::milliseconds time_out, Priority priority, std::weak_ptr<const void> lifetime) = 0;
        virtual void remove_delayed_callback(DelayedHandler& hndl) = 0;

        enum class AwatingType {
            SCHEDULED_ONLY,
            COMPREHENSIVE,
        };
        virtual void wait(AwatingType = AwatingType::COMPREHENSIVE) = 0;

        virtual std::thread::id get_thread_id() const noexcept = 0;
        virtual bool is_working_thread() const noexcept = 0;
        virtual std::string get_name() const = 0;
        virtual size_t size() const = 0;
            
        virtual void set_max_size(size_t size) = 0;
        virtual void clear() = 0;

        virtual std::shared_ptr<ITimersFactory> get_timers_factory() = 0;
    };

    template <typename F>
    auto make_safe_callback(F callback, std::weak_ptr<const void> lifetime);

    template <typename F>
    auto make_safe_callback(F callback, std::weak_ptr<const void> lifetime, std::weak_ptr<ICallbackQueue> callback_queue);

    template <typename F>
    void call_in_worker(F callback, std::weak_ptr<ICallbackQueue> callback_queue, std::optional<std::weak_ptr<const void>> lifetime = std::nullopt);

    enum class ExecutionType {
        SYNCHRONE,
        ASYNCHRONE,
    };

    template <typename BaseClass>
    class SyncCaller {
    public:
        template<typename DerivedClass, typename R, typename... MArgs, typename... Args>
        std::future<R> sync(R(DerivedClass::* method)(MArgs...), Args&&... args);

        template <typename F>
        std::future<std::invoke_result_t<F>> add_sync(F&& callback, std::optional<std::weak_ptr<const void>> lifetime = std::nullopt);
    };

    template <typename BaseClass>
    class AsyncCaller {
    public:
        template<typename DerivedClass, typename R, typename... MArgs, typename... Args>
        std::future<R> async(R(DerivedClass::* method)(MArgs...), Args&&... args);

        template <typename F>
        std::future<std::invoke_result_t<F>> add_async(F&& callback, std::optional<std::weak_ptr<const void>> lifetime = std::nullopt);
    };

    template <typename BaseClass>
    class TypedCaller : public SyncCaller<BaseClass>,
                        public AsyncCaller<BaseClass> {
    public:
        template<typename DerivedClass, typename R, typename... MArgs, typename... Args>
        std::future<R> typed_call(ExecutionType type, R(DerivedClass::* method)(MArgs...), Args&&... args);

        template <typename F>
        std::future<std::invoke_result_t<F>> add_command(ExecutionType type, F&& callback, std::optional<std::weak_ptr<const void>> lifetime = std::nullopt);
    };

} // namespace fem

#include "ICallbackQueue.hpp"
