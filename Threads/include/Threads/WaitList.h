///@file WaitList.h
#pragma once

#include <Patterns/MoveOnlyBase.h>

#include "threads_export.hpp"

#include <future>
#include <list>

namespace mthread {

    /**
         @brief Класс для ожидания завершения группы асинхронных операций.
         @details Позволяет добавлять std::future объекты и ожидать завершения всех добавленных операций.
         Класс является потокобезопасным и реализует паттерн "только для перемещения" (MoveOnly).
    */
    class THREADS_EXPORT WaitList : public pattern::MoveOnlyBase {
    public:
        void add(std::future<void>&& future);   ///< Добавляет future-объект в список ожидания.
        void wait();                            ///< Ожидает завершения всех добавленных операций.

    private:
        std::mutex m_mutex;
        std::list<std::future<void>> m_list;    ///< Список future-объектов для ожидания
    };

} // namespace mthread
