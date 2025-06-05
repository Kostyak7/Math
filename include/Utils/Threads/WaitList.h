#pragma once

#include <Utils/Patterns/MoveOnlyBase.h>

#include "threads_export.hpp"

#include <future>
#include <list>

namespace util::mthrd {

    class THREADS_EXPORT WaitList : public MoveOnlyBase {
    public:
        void add(std::future<void>&& future);
        void wait();

    private:
        std::mutex m_mutex;
        std::list<std::future<void>> m_list;
    };

} // namespace util::mthrd
