#include <Patterns/ItemWithID.h>

pattern::ItemWithID::ItemWithID(const char* id) 
    : m_id(id) 
{
}

const char* pattern::ItemWithID::get_id() const noexcept { 
    return m_id; 
}
