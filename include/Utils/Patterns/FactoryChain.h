#pragma once

#include <memory>
#include <cstring>
#include <vector>

namespace patterns {

    /**
     @brief ��������� ������� ������.
     @details � ������ �������� ���������� ���������� ������ Factory<IItem, Item>, ������� ��������� ���� �� ����� � ������ ������������� �� ��������. \n
     ���������� �� ���� ��������, ������� ������ � ������� ������ �������.
    */
    template <typename IItem>
    class BaseFactory {
    public:
        using item_type = IItem;
        static std::unique_ptr<IItem> Create(const char* id) noexcept;
        static std::vector<const char*> GetIDs() noexcept;
    protected:
        explicit BaseFactory(const char* id) noexcept : m_id(id), m_next(s_head) { s_head = this; }
        virtual std::unique_ptr<IItem> create() const = 0;
    private:
        static BaseFactory* Find(const char* id) noexcept;
    protected:
        const char* m_id;
    private:
        BaseFactory* m_next;
        inline static BaseFactory* s_head = nullptr;
    };

    /**
     @brief ��������� ��������� �����.
     @returns ���������� ��������� �� ������� �������, ��������������� id.
     */
    template<typename IItem>
    std::unique_ptr<IItem> BaseFactory<IItem>::Create(const char* id)  noexcept {
        BaseFactory<IItem>* ptr_factory = BaseFactory<IItem>::Find(id);
        if (ptr_factory)
            return ptr_factory->create();
        return {};
    }

    /**
     @brief ����� ������ ������� �� ��������������.
     @details ����������� �� ���� ����������� ������ � ���� �������, �� ���������� ��������� �� ������ �������, ����� ���������� nullptr.
     */
    template <typename IItem>
    BaseFactory<IItem>* BaseFactory<IItem>::Find(const char* id) noexcept {
        BaseFactory<IItem>* ptr_factory = s_head;
        while (ptr_factory) {
            if (!std::strcmp(ptr_factory->m_id, id))
                break;
            ptr_factory = ptr_factory->m_next;
        }
        return ptr_factory;
    }

    /**
     @brief ���������� json � id ���� ��������� ������.
     @details ����������� �� ���� ����������� ������ � ��������� id � ������.
     */
    template<typename IItem>
    std::vector<const char*> BaseFactory<IItem>::GetIDs() noexcept {
        std::vector<const char*> ids;
        for (BaseFactory* item_id = s_head; item_id != nullptr; item_id = item_id->m_next)
            ids.push_back(item_id->m_id);
        return ids;
    }

    /// ���������� �������, ����������� ���������� ��������� �����.
    template <typename IItem, typename Item>
    class Factory : public BaseFactory<IItem> {
    public:
        explicit Factory(const char* item_id) : BaseFactory<IItem>(item_id) {}
    protected:
        std::unique_ptr<IItem> create() const override {
            return std::make_unique<Item>(this->m_id);
        }
    };

} // namespace patterns
