#pragma once 

#include "IPreconditioner.h"

#include <memory>

namespace math::linal {

    class ILUPreconditioner : public IPreconditioner {
    public:
        void init(const AnyMatrix& matrix) override;
        class Impl;
    };

} // namespace math::linal