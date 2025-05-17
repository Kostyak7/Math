#pragma once 

#include "IPreconditioner.h"

#include <memory>

namespace math::linal {

    class JacobiPreconditioner : public IPreconditioner {
    public:
        void init(const AnyMatrixConstRef& matrix) override;
        class Impl;
    };

} // namespace math::linal