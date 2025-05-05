#pragma once

#include "ICallbackQueue.h"

namespace fem {
    
    class IStoppableCallbackQueue : public ICallbackQueue {
    public:
        virtual void stop() = 0;
    };

} // namespace fem
