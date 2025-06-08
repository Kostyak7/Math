#pragma once

#include "ICallbackQueue.h"

namespace util::mthrd {
    
    class IStoppableCallbackQueue : public ICallbackQueue {
    public:
        virtual void stop() = 0;
    };

} // namespace util::mthrd
