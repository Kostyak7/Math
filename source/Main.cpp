#include "Tests/Math/LinearTests.h"


int main() {

    //tests::base_tests<math::linal::BandMatrix>();

    tests::one_random_test<math::linal::BandMatrix>();

    //tests::diff_sizes_test<math::linal::BandMatrix>();
    
    //tests::diff_isl_test();

    //tests::diff_sparsity_test();

    return 0;
}
