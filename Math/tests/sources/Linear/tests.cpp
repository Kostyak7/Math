#include "Linear.h"

using namespace tests;
using namespace math::linal;

template <class TMatrix>
void base_tests() {
	auto slaes = get_base_SLAEs<TMatrix>();
	SolverTester<TMatrix>::Params params = { true, true, true, false, false, false };
	std::vector<std::shared_ptr<SolverTester<TMatrix>>> testers{
		//create_solver_tester(slaes, SolverFactory::SolverName::Gauss,    SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::CG,       SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::MINRES,   SolverFactory::PredictionerName::None, params), // +
		// 
		//// ���� �������� ������ ����� ������������������� 
		//// BiCGSTAB, �G, MINRES ������������� � ����� � Jacobi ��� ��� �� ��������
		//// ��� ������������� � BiCGSTAB (MINRES?) ����� ����������� ������������ ������������������
		create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::JACOBI, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::JACOBI, params), // +
		// 
		//// CG ����� ������������ � ������������������, �� ������� ��������������
		//// MINRES ������� ������������ ������������� �������������������
		create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::ILU, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::CG,	   SolverFactory::PredictionerName::ILU, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::ILU, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::ILU, params), // +
		// 
		//create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::MG, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::CG,       SolverFactory::PredictionerName::MG, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::MG, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::MG, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::MINRES,   SolverFactory::PredictionerName::MG, params), // +
	};
	TestRunner<TMatrix>::Params runner_params{ false, TestRunner<TMatrix>::SaveType::None };
	TestRunner<TMatrix>::run("base_tests", testers, runner_params);
}
template <class TMatrix>
void one_random_test() {
	TestGenerator<TMatrix> test_generator;
	auto slae = test_generator.get_SLAE(200);
	slae.solution = slae.b;
	slae.b = slae.A * slae.solution;
	std::vector<SLAE<TMatrix>> slaes = { std::move(slae) };
	SolverTester<TMatrix>::Params params = { true, true, true, false, false, false };
	std::vector<std::shared_ptr<SolverTester<TMatrix>>> testers{
		create_solver_tester(slaes, SolverFactory::SolverName::Gauss,    SolverFactory::PredictionerName::None, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::None, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::CG,       SolverFactory::PredictionerName::None, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::None, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::None, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::MINRES,   SolverFactory::PredictionerName::None, params), // +
		 
		create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::JACOBI, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::JACOBI, params), // +
		 
		create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::ILU, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::CG,	     SolverFactory::PredictionerName::ILU, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::ILU, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::ILU, params), // +
		 
		create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::CG,       SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::MINRES,   SolverFactory::PredictionerName::MG, params), // +
	};
	TestRunner<TMatrix>::Params runner_params{ false, TestRunner<TMatrix>::SaveType::OnlyTime };
	TestRunner<TMatrix>::run("one_random_test", testers, runner_params);
}
template <class TMatrix>
void diff_sizes_test() {
	std::vector<size_t> sizes = { 100, 200, 300, 400, 500, 1000, 1500, 2000, 3000, 4000, 5000 };
	TestGenerator<TMatrix> test_generator;
	auto slaes = test_generator.get_SLAEs(sizes);
	SolverTester<TMatrix>::Params params = { false, false, false, false, false, false };
	std::vector<std::shared_ptr<SolverTester<TMatrix>>> testers{
		//create_solver_tester(slaes, SolverFactory::SolverName::Gauss,    SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::CG,       SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::MINRES,   SolverFactory::PredictionerName::None, params), // +
		// 
		//create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::JACOBI, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::JACOBI, params), // +
		// 
		//create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::ILU, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::CG,	     SolverFactory::PredictionerName::ILU, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::ILU, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::ILU, params), // +
		// 
		create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::CG,       SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::MINRES,   SolverFactory::PredictionerName::MG, params), // +
	};
	TestRunner<TMatrix>::Params runner_params{ true, TestRunner<TMatrix>::SaveType::OnlyTime };
	TestRunner<TMatrix>::run("diff_sizes_test", testers, runner_params);
}
void diff_isl_test() {
	using TMatrix = ::BandMatrix;
	const size_t matrix_size = 1000;
	std::vector<size_t> isl_sizes = { 10, 25, 50, 75, 100, 150, 300, 500 };
	std::vector<SLAE<TMatrix>> slaes;
	TestGenerator<TMatrix> test_generator;
	for (size_t i = 0; i < isl_sizes.size(); ++i) {
		slaes.emplace_back(test_generator.get_SLAE(matrix_size, isl_sizes[i]));
	}
	SolverTester<TMatrix>::Params params = { false, false, false, false, false, false };
	std::vector<std::shared_ptr<SolverTester<TMatrix>>> testers{
		//create_solver_tester(slaes, SolverFactory::SolverName::Gauss,    SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::CG,       SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::None, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::MINRES,   SolverFactory::PredictionerName::None, params), // +
		//
		//create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::JACOBI, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::JACOBI, params), // +
		//
		//create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::ILU, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::CG,	     SolverFactory::PredictionerName::ILU, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::ILU, params), // +
		//create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::ILU, params), // +
		//
		create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::CG,       SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::MINRES,   SolverFactory::PredictionerName::MG, params), // +
	};
	TestRunner<TMatrix>::Params runner_params{ true, TestRunner<TMatrix>::SaveType::OnlyTime };
	TestRunner<TMatrix>::run("diff_isl_test", testers, runner_params);
}
void diff_sparsity_test() {
	using TMatrix = ::BandMatrix;
	const size_t matrix_size = 1000;
	const size_t isl = matrix_size / 5;
	std::vector<SLAE<TMatrix>> slaes;
	TestGenerator<TMatrix> test_generator;
	for (int i = 0; i < 5; ++i) {
		double sparce = 0.1 + 0.2 * i;
		slaes.emplace_back(test_generator.get_SLAE(matrix_size, isl, sparce));
	}
	SolverTester<TMatrix>::Params params = { false, false, false, false, false, false };
	std::vector<std::shared_ptr<SolverTester<TMatrix>>> testers{
		create_solver_tester(slaes, SolverFactory::SolverName::Gauss,    SolverFactory::PredictionerName::None, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::None, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::CG,       SolverFactory::PredictionerName::None, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::None, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::None, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::MINRES,   SolverFactory::PredictionerName::None, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::JACOBI, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::JACOBI, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::ILU, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::CG,	     SolverFactory::PredictionerName::ILU, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::ILU, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::ILU, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::BiCGSTAB, SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::CG,       SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::FGMRES,   SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::GMRES,    SolverFactory::PredictionerName::MG, params), // +
		create_solver_tester(slaes, SolverFactory::SolverName::MINRES,   SolverFactory::PredictionerName::MG, params), // +
	};
	TestRunner<TMatrix>::Params runner_params{ true, TestRunner<TMatrix>::SaveType::OnlyTime };
	TestRunner<TMatrix>::run("diff_sparsity_test", testers, runner_params);
}

TEST(LINEAR, BASE_SLAES) {
    std::thread([]() {
    base_tests<::BandMatrix>();  
    }).join();
}

TEST(LINEAR, ONE_RANDOM_SLAE) {
    std::thread([]() {
	one_random_test<::BandMatrix>();

	//diff_sizes_test<::BandMatrix>();

	//diff_isl_test();

	//diff_sparsity_test();    
    }).join();
}

TEST(LINEAR, DIFF_SIZES_SLAES) {
    std::thread([]() {
	//diff_sizes_test<::BandMatrix>();
    }).join();
}

TEST(LINEAR, BAND_MATRIX_DIFF_ISL) {
    std::thread([]() {
	//diff_isl_test();
    }).join();
}

TEST(LINEAR, BAND_MATRIX_DIFF_SPARSITY) {
    std::thread([]() {
	//diff_sparsity_test();    
    }).join();
}
