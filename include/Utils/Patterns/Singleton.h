#pragma once

#include "NonCopyable.h"

namespace pattern {

    template <typename DerivedT>
    class Singleton : public NonCopyable {
    public:
        static DerivedT& Instance() {
            thread_local DerivedT instance;
            return instance;
        }

    protected:
        Singleton() = default;
        virtual ~Singleton() = default;
    };

} // namespace pattern
