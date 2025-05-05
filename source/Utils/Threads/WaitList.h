#pragma once

#include <Utils/Patterns/MoveOnlyBase.h>

#include <future>
#include <list>

namespace fem {

    class WaitList : public MoveOnlyBase {
    public:
        void add(std::future<void>&& future);
        void wait();

    private:
        std::mutex m_mutex;
        std::list<std::future<void>> m_list;
    };

} // namespace fem
