#pragma once

#include <functional>
#include <memory>
#include <unordered_map>
#include <string>

namespace pattern {

    template <typename BaseT, typename KeyT = std::string>
    class Factory {
    public:
        using Creator = std::function<std::unique_ptr<BaseT>()>;

        bool bind(const KeyT& key, Creator creator) {
            return m_creators.emplace(key, std::move(creator)).second;
        }

        std::unique_ptr<BaseT> create(const KeyT& key) const {
            auto it = m_creators.find(key);
            return it != m_creators.end() ? it->second() : nullptr;
        }

    private:
        std::unordered_map<KeyT, Creator> m_creators;
    };

} // namespace pattern
