///@file IStoppableCallbackQueue.h
#pragma once

#include "ICallbackQueue.h"

namespace mthread {
    
    class IStoppableCallbackQueue : public ICallbackQueue {
    public:
        virtual void stop() = 0;
    };

} // namespace mthread
