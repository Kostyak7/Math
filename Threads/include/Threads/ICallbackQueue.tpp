#pragma once

#include <future>
#include <optional>
#include <type_traits>

template<typename F>
auto mthread::make_safe_callback(F callback, std::weak_ptr<const void> lifetime) {
    return [callback = std::move(callback), lifetime = std::move(lifetime)](auto&&... args) -> void {
        if (auto lock = lifetime.lock()) {
            callback(std::forward<decltype(args)>(args)...);
        }
        };
}

template<typename F>
auto mthread::make_safe_callback(F callback, std::weak_ptr<const void> lifetime, std::weak_ptr<ICallbackQueue> callback_queue) {
    return [callback = std::move(callback), lifetime = std::move(lifetime), callback_queue = std::move(callback_queue)](auto&&... args) -> void {
        if (auto cq = callback_queue.lock()) {
            if (cq->is_working_thread()) {
                if (auto lock = lifetime.lock()) {
                    callback(std::forward<decltype(args)>(args)...);
                }
            }
            else {
                cq->add_callback([callback, args...]() mutable { callback(std::forward<decltype(args)>(args)...); }, std::move(lifetime));
            }
        }
        };
}

template <typename F>
void mthread::call_in_worker(F callback, std::weak_ptr<ICallbackQueue> callback_queue, std::optional<std::weak_ptr<const void>> lifetime) {
    auto cb = callback_queue.lock();
    if (!cb) {
        return;
    }

    if (cb->is_working_thread()) {
        callback();
    }
    else {
        if (lifetime.has_value()) {
            cb->add_callback(std::move(callback), std::move(lifetime.value()));
        }
        else {
            cb->add_callback(std::move(callback));
        }
    }
}

template <typename BaseClass>
template<typename DerivedClass, typename R, typename... MArgs, typename... Args>
std::future<R> mthread::SyncCaller<BaseClass>::sync(R(DerivedClass::* method)(MArgs...), Args&&... args) {
    static_assert(std::is_base_of_v<BaseClass, DerivedClass>, "Mismatch deriving");

    return add_sync([self = static_cast<BaseClass*>(this)->shared_from_this(), method, args = std::make_tuple(std::forward<Args>(args)...)]() mutable -> R {
        return std::apply([derived_self = std::static_pointer_cast<DerivedClass>(self), method](auto&&... unpacked_args) -> R {
            (derived_self.get()->*method)(std::forward<decltype(unpacked_args)>(unpacked_args)...);
            }, std::move(args));
        });
}

template <typename BaseClass>
template <typename F>
std::future<std::invoke_result_t<F>> mthread::SyncCaller<BaseClass>::add_sync(F&& callback, std::optional<std::weak_ptr<const void>> lifetime) {
    static_assert(std::is_base_of_v<ICallbackQueue, BaseClass>, "Caller does not belong to ICallbackQueue");

    using result_t = std::invoke_result_t<F>;
    auto promise = std::make_shared<std::promise<result_t>>();
    auto future = promise->get_future();

    auto to_cbq = [callback = std::move(callback), promise]() mutable {
        try {
            if constexpr (std::is_void_v<result_t>) {
                callback();
                promise->set_value();
            }
            else {
                promise->set_value(callback());
            }
        }
        catch (...) {
            promise->set_exception(std::current_exception());
        }
        };

    auto self = static_cast<BaseClass*>(this);
    if (lifetime.has_value()) {
        self->add_callback(std::move(to_cbq), std::move(*lifetime));
    }
    else {
        self->add_callback(std::move(to_cbq));
    }

    return future;
}

template <typename BaseClass>
template<typename DerivedClass, typename R, typename... MArgs, typename... Args>
std::future<R> mthread::AsyncCaller<BaseClass>::async(R(DerivedClass::* method)(MArgs...), Args&&... args) {
    static_assert(std::is_base_of_v<BaseClass, DerivedClass>, "Mismatch deriving");

    return add_async([self = static_cast<BaseClass*>(this)->shared_from_this(), method, args = std::make_tuple(std::forward<Args>(args)...)]() mutable -> R {
        return std::apply([derived_self = std::static_pointer_cast<DerivedClass>(self), method](auto&&... unpacked_args) -> R {
            return (derived_self.get()->*method)(std::forward<decltype(unpacked_args)>(unpacked_args)...);
            }, std::move(args));
        });
}

template <typename BaseClass>
template <typename F>
std::future<std::invoke_result_t<F>> mthread::AsyncCaller<BaseClass>::add_async(F&& callback, std::optional<std::weak_ptr<const void>> lifetime) {
    static_assert(std::is_base_of_v<ICallbackQueue, BaseClass>, "Caller does not belong to ICallbackQueue");

    auto self = static_cast<BaseClass*>(this)->shared_from_this();
    if (self->is_working_thread()) {
        return std::async(std::launch::async, std::move(callback));
    }

    using result_t = std::invoke_result_t<F>;
    auto future_ptr = std::make_shared<std::future<result_t>>();

    call_in_worker([callback = std::move(callback), future_ptr]() mutable {
        *future_ptr = std::async(std::launch::async, std::move(callback));
        }, self, std::move(lifetime));

    return std::async(std::launch::deferred, [future_ptr]() -> result_t {
        return future_ptr->get();
        });
}

template <typename BaseClass>
template<typename DerivedClass, typename R, typename... MArgs, typename... Args>
std::future<R> mthread::TypedCaller<BaseClass>::typed_call(ExecutionType type, R(DerivedClass::* method)(MArgs...), Args&&... args) {
    if (type == ExecutionType::SYNCHRONE) {
        return sync(method, std::forward<Args>(args)...);
    }
    else {
        return async(method, std::forward<Args>(args)...);
    }
}

template <typename BaseClass>
template <typename F>
std::future<std::invoke_result_t<F>> mthread::TypedCaller<BaseClass>::add_command(ExecutionType type, F&& callback, std::optional<std::weak_ptr<const void>> lifetime) {
    if (type == ExecutionType::SYNCHRONE) {
        return add_sync(std::move(callback), std::move(lifetime));
    }
    else {
        return add_async(std::move(callback), std::move(lifetime));
    }
}
