#pragma once 

#include "../../Common.h"

#include <Math/Linear/Solvers/BiCGSTABSystemSolver.h>
#include <Math/Linear/Solvers/CGSystemSolver.h>
#include <Math/Linear/Solvers/GaussSystemSolver.h>
#include <Math/Linear/Solvers/GMRESSystemSolver.h>
#include <Math/Linear/Solvers/FGMRESSystemSolver.h>
#include <Math/Linear/Solvers/MINRESSystemSolver.h>
#include <Math/Linear/Solvers/SORSystemSolver.h>

#include <Math/Linear/Solvers/Preconditioners/ILUPreconditioner.h>
#include <Math/Linear/Solvers/Preconditioners/JacobiPreconditioner.h>
#include <Math/Linear/Solvers/Preconditioners/MGPreconditioner.h>

namespace tests {

	void print_vector(const math::linal::DVector& vec);

	template <class TMatrix>
	void print_matrix(const TMatrix& matrix);

	template <class TMatrix>
	struct SLAE {
		TMatrix A;
		math::linal::DVector solution;
		math::linal::DVector b;
	};

	template <class TMatrix>
	class SolverTester {
	public:
		struct TestResult {
			size_t index = 0;
			bool success = false;

			math::linal::DVector x;
			std::chrono::microseconds duration = std::chrono::microseconds{ 0 };
			size_t operation = 0;
			size_t memory = 0;
		};

		struct Params {
			bool print_all_results = true;
			bool print_every_result = false;
			bool print_slae = false;

			bool print_time = true;
			bool print_operation_count = false;
			bool print_memory = false;
		};

	public:
		SolverTester(const std::vector<SLAE<TMatrix>>& slaes,
			std::unique_ptr<math::linal::ILinearSystemSolver>&& solver,
			const std::string& solver_name,
			const Params& params = {});
		std::string get_solver_name() const;
		std::vector<TestResult> run();

	private:
		TestResult test(size_t index);
		void print_test_result(const TestResult& res) const;
		bool only_one_print() const;

	private:
		const std::vector<SLAE<TMatrix>>& m_slaes;
		const std::unique_ptr<math::linal::ILinearSystemSolver> m_solver;
		const std::string m_solver_name;

		const Params m_params;
	};

	class ITestGenerator {
	public:
		ITestGenerator();

	protected:
		double get_random_ranged_double(double lower_bound = -1000., double upper_bound = 1000., bool can_be_zero = true);
		double get_random_double();
		math::linal::DVector get_random_vector(size_t n, const double& sparsity);

	protected:
		std::random_device rd;
		std::mt19937 gen;
		std::uniform_real_distribution<double> sparsity_dist;
	};

	class SolverFactory {
	public:
		enum class SolverName {
			BiCGSTAB,
			CG,
			FGMRES,
			Gauss,
			GMRES,
			MINRES,
			SOR,
		};

		enum class PredictionerName {
			None,
			JACOBI,
			ILU,
			MG,
		};

	public:
		static std::unique_ptr<math::linal::ILinearSystemSolver> get_solver(SolverName name, PredictionerName pred);
		template <class ...Args> static std::unique_ptr<math::linal::IPreconditioner> get_predictioner(PredictionerName name, Args&&... args);
		static std::string get_full_solver_name(SolverName name, PredictionerName pred);
		static std::string get_solver_name(SolverName name);
		static std::string get_predictioner_name(PredictionerName name);

	private:
		template <class SolverT, class ...Args>
		static std::unique_ptr<math::linal::ILinearSystemSolver> make_solver(std::unique_ptr<math::linal::IPreconditioner> preconditioner, Args&&... args);

	private:
		static std::unordered_map<SolverName, std::string> s_solver_names;
		static std::unordered_map<PredictionerName, std::string> s_predictioner_names;
	};

	template <class TMatrix>
	std::shared_ptr<SolverTester<TMatrix>> create_solver_tester(const std::vector<SLAE<TMatrix>>& slaes, SolverFactory::SolverName name, SolverFactory::PredictionerName pred, typename SolverTester<TMatrix>::Params const& params = {});

	template <class TMatrix>
	class TestGenerator : public ITestGenerator {
	public:
		std::vector<SLAE<TMatrix>> get_SLAEs(const std::vector<size_t>& sizes, double sparsity, bool async) = delete;
		SLAE<TMatrix> get_SLAE(size_t n, size_t isl, double sparsity) = delete;
		SLAE<TMatrix> get_SLAE(size_t n, double sparsity) = delete;
	};

	template <class TMAtrix>
	std::vector<SLAE<TMAtrix>> get_base_SLAEs();

	template <class TMatrix>
	class TestRunner {
	public:
		enum class SaveType {
			None,
			All,
			OnlyTime,
		};

		struct Params {
			bool async = false;
			SaveType save_type = SaveType::None;
		};

		static void run(const std::string& name, const std::vector<std::shared_ptr<SolverTester<TMatrix>>>& testers, const Params& params = {});

	private:
		static auto async_run(const std::string& name, const std::vector<std::shared_ptr<SolverTester<TMatrix>>>& testers, const Params& params);
		static void save_all_tests(const std::string& name, const std::unordered_map<std::string, std::vector<typename SolverTester<TMatrix>::TestResult>>& results);
		static void save_only_time(const std::string& name, const std::unordered_map<std::string, std::vector<typename SolverTester<TMatrix>::TestResult>>& results);
	};

} // namespace tests

#include "Linear.tpp"
