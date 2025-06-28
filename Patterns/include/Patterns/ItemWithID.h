#pragma once

namespace pattern {

    class ItemWithID {
    public:
        ItemWithID(const char* id);
        const char* get_id() const noexcept;
        
    protected:
        const char* m_id;
    };

} // namespace pattern
