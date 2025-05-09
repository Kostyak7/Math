#pragma once 

#include "IPreconditioner.h"

#include <Utils/Math/Linear/BandMatrix.h>

#include <memory>

namespace math::linal {

    class MGPreconditioner : public IPreconditioner {
    public:
        enum class SmootherType {
            JACOBI,
            GAUSS_SEIDEL,
            ILU0
        };

        enum class CycleType {
            V_CYCLE,
            W_CYCLE,
            F_CYCLE
        };

        struct Params {
            size_t min_coarse_size = 50;    // ����������� ������ ������ �����
            size_t max_levels = 10;         // ����. ����� �������
            size_t n_pre_smooth = 2;        // ����� ����-�����������
            size_t n_post_smooth = 2;       // ����� ����-�����������
            double strong_threshold = 0.25; // ����� ��� ������� ������
            SmootherType smoother_type = SmootherType::JACOBI;
            CycleType cycle_type = CycleType::F_CYCLE;
        };

        explicit MGPreconditioner(const Params& params = {});
        void init(const AnyMatrix& matrix) override;

        class Impl;

    private:
        Params m_params;
    };

} // namespace math::linal