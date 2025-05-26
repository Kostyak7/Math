#pragma once

#include <list>
#include <mutex>

namespace fem {

    /**
     * @brief Отслеживает соединение InitiatorType c владельцем.
     * Позволяет выполнять команды InitiatorType напрямую, а не через очередь.
     * Позволяет выполнять команды владельца напрямую, а не через очередь.
     */
    template <class InitiatorType>
    class Connection {
    public:
        class Scope {
        public:
            Scope(Connection& stack, InitiatorType* initiator = nullptr)
                : m_stack(stack)
            {
                m_stack.open_scope(initiator);
            }
            ~Scope() { m_stack.close_scope(); }
        private:
            Connection& m_stack;
        };

        Scope scope(InitiatorType* initiator = nullptr);
        bool in_scope();
        void open_scope(InitiatorType* initiator = nullptr);
        void close_scope();

        bool make_friend(std::thread::id thread);
        void lose_friend(std::thread::id thread);

        InitiatorType* get_initiator();

    private:
        std::recursive_mutex m_mutex;
        std::list<std::thread::id> m_friends;
        InitiatorType* m_initiator = nullptr;
        size_t m_depth = 0;
    };

} // namespace fem

#include "Connection.tpp"
