template <class InitiatorType>
mthread::Connection<InitiatorType>::Scope mthread::Connection<InitiatorType>::scope(InitiatorType* initiator) { return Scope(*this, initiator); }

template <class InitiatorType>
bool mthread::Connection<InitiatorType>::in_scope() {
    std::unique_lock lock(m_mutex);
    return m_depth != 0 &&
        std::any_of(m_friends.begin(), m_friends.end(),
            [this](auto& key) {return key == std::this_thread::get_id(); }
        );
}

template <class InitiatorType>
void mthread::Connection<InitiatorType>::open_scope(InitiatorType* initiator) {
    std::unique_lock lock(m_mutex);
    if (m_depth == 0) {
        m_friends.push_back(std::this_thread::get_id());
        m_initiator = initiator;
    }
    ++m_depth;
}

template <class InitiatorType>
void mthread::Connection<InitiatorType>::close_scope() {
    std::unique_lock lock(m_mutex);
    --m_depth;
    if (m_depth == 0) {
        m_initiator = nullptr;
        m_friends.clear();
    }
}

template <class InitiatorType>
bool mthread::Connection<InitiatorType>::make_friend(std::thread::id thread) {
    std::unique_lock lock(m_mutex);
    if (in_scope()) {
        m_friends.push_back(thread);
        return true;
    }
    return false;
}

template <class InitiatorType>
void mthread::Connection<InitiatorType>::lose_friend(std::thread::id thread) {
    std::unique_lock lock(m_mutex);
    if (in_scope())
        m_friends.remove(thread);
}

template <class InitiatorType>
InitiatorType* mthread::Connection<InitiatorType>::get_initiator() {
    return m_initiator;
}
