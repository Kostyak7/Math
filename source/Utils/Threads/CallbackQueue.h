#pragma once

#include "BlockingQueue.h"
#include "IStoppableCallbackQueue.h"

namespace fem {

    class CallbackQueue : public IStoppableCallbackQueue,
                          public std::enable_shared_from_this<CallbackQueue> {
    public:
        explicit CallbackQueue(std::string name);
        CallbackQueue(std::shared_ptr<ITimersFactory> timers_factory, std::string name);
        CallbackQueue(std::string name, size_t max_size, OnOverflow on_overflow_callback = nullptr, OnBecomeNormal on_become_normal_callback = nullptr);
        CallbackQueue(std::shared_ptr<ITimersFactory> timers_factory, std::string name, size_t max_size, OnOverflow on_overflow_callback = nullptr, OnBecomeNormal on_become_normal_callback = nullptr);
        ~CallbackQueue() override;

        template <typename T>
        std::shared_ptr<T> shared_from_base() {
            return std::static_pointer_cast<T>(shared_from_this());
        }

        template <typename T>
        std::weak_ptr<T> weak_from_base() {
            return shared_from_base<T>();
        }

        void add_callback(std::function<void()> callback) override;
        void add_callback(std::function<void()> callback, std::weak_ptr<const void> lifetime) override;

        bool try_add_callback(std::function<void()> callback) override;
        bool try_add_callback(std::function<void()> callback, std::weak_ptr<const void> lifetime) override;

        DelayedHandler add_delayed_callback(std::function<void()> callback, std::chrono::milliseconds time_out) override;
        DelayedHandler add_delayed_callback(std::function<void()> callback, std::chrono::milliseconds time_out, std::weak_ptr<const void> lifetime) override;
        void remove_delayed_callback(DelayedHandler& hndl) override;

        void wait(AwatingType = AwatingType::COMPREHENSIVE) override;

        std::thread::id get_thread_id() const noexcept override;
        bool is_working_thread() const noexcept override;
        std::string get_name() const override;
        size_t size() const override;

        void set_max_size(size_t size) override;
        void clear() override;

        void stop() override;

        std::shared_ptr<ITimersFactory> get_timers_factory() override;

    protected:
        void start();

    private:
        void main_loop();

    protected: // Подсматривая в этот флаг, наследник будет срочно завершать все текущие задачи
        std::atomic<bool> m_stopped = false; ///< Флаг взводится, когда нужно как можно скорее "убить" поток.

    private:
        const std::string m_name;
        BlockingQueue<std::function<void()>> m_queue;
        std::shared_ptr<ITimersFactory> m_timers_factory;
        std::thread m_queue_thread;
        std::atomic<std::thread::id> m_therad_id;
    };

} // namespace fem
