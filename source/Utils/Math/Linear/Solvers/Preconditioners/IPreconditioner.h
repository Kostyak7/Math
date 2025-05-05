#pragma once 

#include <Utils/Math/Linear/AnyMatrix.h>

#include <memory>

namespace math::linal {

    class IPreconditioner {
    public:
        virtual ~IPreconditioner() = default;
        virtual void init(const AnyMatrix& matrix) = 0;

        FVector apply(const FVector& x) const {
            return m_impl->apply(x);
        }

        class IImpl {
        public:
            virtual ~IImpl() = default;
            virtual FVector apply(const FVector& x) const = 0;
        };

    protected:
        std::unique_ptr<IImpl> m_impl;
    }; 

} // namespace math::linal