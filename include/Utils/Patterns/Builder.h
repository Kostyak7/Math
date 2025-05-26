#pragma once

#include <memory>

namespace pattern {

    template <typename T>
    class Builder {
    public:
        virtual ~Builder() = default;
        virtual std::unique_ptr<T> build() const = 0;
    };

} // namespace pattern
