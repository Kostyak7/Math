///@file Barrier.h
#pragma once 

#include <Patterns/MoveOnlyBase.h>

#include "threads_export.hpp"

#include <condition_variable>

namespace mthread {

    /**
        @brief Барьер — инструмент синхронизации потоков.
        @details Барьер представляет собой механизм синхронизации, который заставляет
            несколько потоков ожидать друг друга перед продолжением выполнения.
            Все потоки блокируются при вызове `arrive_and_wait()`, пока
            заданное количество потоков не достигнет барьера.
 
        @note Класс является move-only.
        @note Потокобезопасен.
        @note Аналог `std::barrier` из C++20, но реализован на C++17.
    */
    class THREADS_EXPORT Barrier : public pattern::MoveOnlyBase {
    public:
        explicit Barrier(size_t expected) noexcept;

        void arrive_and_wait();             ///< Блокирует текущий поток до тех пор, пока все потоки не достигнут барьера
        void arrive_and_drop();             ///< Уменьшает счетчик ожидаемых потоков (поток покидает барьер без ожидания)
        void arrive();                      ///< Поток приходит, но не блокируется
        void wait();                        ///< Поток ожидает, пока барьер не будет достигнут

        void break_barrier();               ///< Разрушает барьер, пробуждает все ожидающие потоки с ошибкой

        size_t get_expected() const;        ///< Возвращает текущее количество потоков, которые должны достичь барьера
        size_t get_arrived() const;         ///< Возвращает количество потоков, которые уже достигли барьера

    private:
        const size_t m_total;               ///<  Сколько должно дойти
        std::atomic<size_t> m_arrived;      ///<  Количество потоков, которые достигли барьера
        std::atomic<size_t> m_generation;   ///<  Номер текущего поколения барьера
        std::mutex m_mutex;                
        std::condition_variable m_cv; 
        bool m_broken;                      ///< Флаг, указывающий, что барьер сломан
    };

} // namespace mthread
