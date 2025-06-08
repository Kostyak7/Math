#pragma once

namespace pattern {

    class ItemWithID {
    public:
        ItemWithID(const char* id) : m_id(id) {}
        const char* get_id() const noexcept { return m_id; }
    protected:
        const char* m_id;
    };

} // namespace pattern
