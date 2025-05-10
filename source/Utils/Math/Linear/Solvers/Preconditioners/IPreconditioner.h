#pragma once 

#include <Utils/Math/Linear/AnyMatrix.h>

#include <memory>

namespace math::linal {

    class IPreconditioner {
    public:
        virtual ~IPreconditioner() = default;
        virtual void init(const AnyMatrix& matrix) = 0;

        DVector apply(const DVector& x) const { // по дефолту это только левое предобуславливание
            return m_impl->apply(x);
        }

        class IImpl {
        public:
            virtual ~IImpl() = default;
            virtual DVector apply(const DVector& x) const = 0;
        };

    protected:
        std::unique_ptr<IImpl> m_impl;
    }; 

} // namespace math::linal