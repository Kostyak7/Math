#pragma once

#include "MoveOnlyBase.h"

namespace pattern {

    template <typename DerivedT>
    class Singleton : public MoveOnlyBase {
    public:
        static DerivedT& Instance() {
            static DerivedT instance;
            return instance;
        }

    protected:
        Singleton() = default;
        virtual ~Singleton() = default;
    };

} // namespace pattern
