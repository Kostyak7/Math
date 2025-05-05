#pragma once

#include "MoveOnlyBase.h"

#include <atomic>
#include <mutex>

namespace pattern {

    template <typename DerivedT>
    class ThreadSafeSingleton : public MoveOnlyBase {
    public:
        static DerivedT& Instance() {
            DerivedT* instance = s_instance.load(std::memory_order_acquire);
            if (!instance) {
                std::lock_guard<std::mutex> lock(s_mutex);
                instance = s_instance.load(std::memory_order_relaxed);
                if (!instance) {
                    instance = new DerivedT;
                    s_instance.store(instance, std::memory_order_release);
                }
            }
            return *instance;
        }

    protected:
        ThreadSafeSingleton() = default;
        virtual ~ThreadSafeSingleton() = default;

    private:
        inline static std::atomic<DerivedT*> s_instance = nullptr;
        inline static std::mutex s_mutex;
    };


} // namespace pattern
