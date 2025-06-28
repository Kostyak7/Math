#pragma once

namespace pattern {

    class MoveOnlyBase {
    public:
        MoveOnlyBase() = default;

        MoveOnlyBase(MoveOnlyBase&&) noexcept = default;
        MoveOnlyBase& operator=(MoveOnlyBase&&) noexcept = default;

        MoveOnlyBase(const MoveOnlyBase&) = delete;
        MoveOnlyBase& operator=(const MoveOnlyBase&) = delete;

        virtual ~MoveOnlyBase() = default;        
    };

} // namespace pattern
