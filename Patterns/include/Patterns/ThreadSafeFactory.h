#pragma once

#include <functional>
#include <memory>
#include <mutex>
#include <string>

namespace pattern {

    template <typename BaseT, typename KeyT = std::string>
    class ThreadSafeFactory {
    public:
        using Creator = std::function<std::unique_ptr<BaseT>()>;

        bool bind(const KeyT& key, Creator creator) {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_creators.emplace(key, std::move(creator)).second;
        }

        std::unique_ptr<BaseT> create(const KeyT& key) const {
            std::lock_guard<std::mutex> lock(m_mutex);
            auto it = m_creators.find(key);
            return it != m_creators.end() ? it->second() : nullptr;
        }

    private:
        std::unordered_map<KeyT, Creator> m_creators;
        mutable std::mutex m_mutex;
    };

} // namespace pattern
