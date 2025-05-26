#pragma once 

#include <Utils/Patterns/MoveOnlyBase.h>
#include <condition_variable>

namespace util::mt {

    class Barrier : public MoveOnlyBase {
    public:
        explicit Barrier(const size_t expected) noexcept;

        void arrive_and_wait() noexcept;

    private:
        const size_t m_total;              ///<  Сколько должно дойти
        std::atomic<size_t> m_arrived;     ///<  Количество потоков, которые достигли барьера
        std::atomic<size_t> m_generation;  ///<  Номер текущего поколения барьера
        std::mutex m_mutex;                
        std::condition_variable m_cv; 
    };

} // namespace util::mt
