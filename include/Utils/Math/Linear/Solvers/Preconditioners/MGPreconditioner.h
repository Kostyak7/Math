#pragma once 

#include "IPreconditioner.h"

#include <Utils/Math/Linear/BandMatrix.h>

#include <memory>

namespace math::linal {

    class MATH_EXPORT MGPreconditioner : public IPreconditioner {
    public:
        enum class MGType {
            Algebraic,
            Geometric,
        };

        enum class SmootherType {
            JACOBI,
            GAUSS_SEIDEL,
            ILU0,
        };

        enum class CycleType {
            V_CYCLE,
            W_CYCLE,
            F_CYCLE,
        };

        struct Params {
            MGType mg_type = MGType::Algebraic;
            size_t min_coarse_size = 50;    // ����������� ������ ������ �����
            size_t max_levels = 10;         // ����. ����� �������
            size_t n_pre_smooth = 2;        // ����� ����-�����������
            size_t n_post_smooth = 2;       // ����� ����-�����������
            double strong_threshold = 0.25; // ����� ��� ������� ������
            SmootherType smoother_type = SmootherType::GAUSS_SEIDEL;
            CycleType cycle_type = CycleType::V_CYCLE;
        };

        explicit MGPreconditioner(const Params& params = {});
        void init(const AnyMatrixConstRef& matrix) override;

        class Impl;

    private:
        Params m_params;
    };

} // namespace math::linal