///@file Connection.h
#pragma once

#include <list>
#include <mutex>

namespace mthread {

    /**
        @brief Отслеживает соединение между инициатором (InitiatorType) и владельцем.
        @details Позволяет выполнять команды напрямую между инициатором и владельцем,
                 минуя очередь команд. Обеспечивает безопасное взаимодействие между потоками.

        @tparam InitiatorType Тип инициатора соединения.
    */
    template <class InitiatorType>
    class Connection {
    public:

        /**
            @brief RAII-обертка для управления областью видимости соединения.
            @details При создании открывает область видимости соединения,
                     при разрушении - автоматически закрывает.
        */
        class Scope {
        public:
            Scope(Connection& stack, InitiatorType* initiator = nullptr)
                : m_stack(stack)
            {
                m_stack.open_scope(initiator);
            }
            ~Scope() { m_stack.close_scope(); }
        private:
            Connection& m_stack;        ///< Ссылка на родительский объект Connection.
        };

        Scope scope(InitiatorType* initiator = nullptr);        ///< Создает новую область видимости соединения
        bool in_scope();                ///< Проверяет, находится ли текущий поток в области видимости соединения
        void open_scope(InitiatorType* initiator = nullptr);    ///< Открывает новую область видимости соединения
        void close_scope();             ///< Закрывает текущую область видимости соединения

        bool make_friend(std::thread::id thread);   ///< Добавляет поток в список "дружественных" потоков
        void lose_friend(std::thread::id thread);   ///< Удаляет поток из списка "дружественных" потоков

        InitiatorType* get_initiator();             ///< Возвращает текущего инициатора соединения

    private:
        std::recursive_mutex m_mutex;               
        std::list<std::thread::id> m_friends;       ///< Список "дружественных" потоков
        InitiatorType* m_initiator = nullptr;       ///< Текущий инициатор соединения
        size_t m_depth = 0;                         ///< Глубина вложенности областей видимости
    };

} // namespace mthread

#include "Connection.tpp"
