#pragma once

#include <iostream>
#include <future>
#include <fstream>


template <class TMatrix>
void tests::print_matrix(const TMatrix& matrix) {
	std::cout << "[";
	for (size_t i = 0; i < matrix.get_height() - 1; ++i) {
		print_vector(matrix[i]);
		std::cout << ",\n";
	}
	if (matrix.get_height())
		print_vector(matrix[matrix.get_height() - 1]);
	std::cout << "]";
}


template <class TMatrix>
tests::SolverTester<TMatrix>::SolverTester(const std::vector<SLAE<TMatrix>>& slaes, std::unique_ptr<math::linal::ILinearSystemSolver>&& solver, const std::string& solver_name, const Params& params)
	: m_slaes(slaes)
	, m_solver(std::move(solver))
	, m_solver_name(solver_name)
	, m_params(params)
{
}

template <class TMatrix>
std::string tests::SolverTester<TMatrix>::get_solver_name() const {
	return m_solver_name;
}

template <class TMatrix>
std::vector<typename tests::SolverTester<TMatrix>::TestResult> tests::SolverTester<TMatrix>::run() {
	std::vector<TestResult> results;
	if (m_params.print_all_results) {
		std::cout << "\n________________________________________________________________________________________\n";
		std::cout << "Solver name: " << m_solver_name << "\n\n";
		if (only_one_print()) {
			if (m_params.print_time) {
				std::cout << "time: ";
			}
			if (m_params.print_operation_count) {
				std::cout << "operations: ";
			}
			if (m_params.print_memory) {
				std::cout << "memory: ";
			}
		}
	}
	for (size_t i = 0; i < m_slaes.size(); ++i) {
		results.emplace_back(test(i));
		if (m_params.print_every_result) {
			std::cout << "\n######################## Test - " << i << " ########################\n";
			print_test_result(results.back());
		}
		else if (m_params.print_all_results && only_one_print()) {
			if (m_params.print_time) {
				std::cout << results.back().duration.count() / 1000.f << " ";
			}
			if (m_params.print_operation_count) {
				std::cout << results.back().operation << " ";
			}
			if (m_params.print_memory) {
				std::cout << results.back().memory << " ";
			}

		}
	}
	if (m_params.print_all_results && !only_one_print()) {
		if (m_params.print_time) {
			std::cout << "\ntime: ";
			for (const auto& res : results) {
				std::cout << res.duration.count() / 1000.f << ' ';
			}

		}
		if (m_params.print_operation_count) {
			std::cout << "\noperations: ";
			for (const auto& res : results) {
				std::cout << res.operation << ' ';
			}
		}
		if (m_params.print_memory) {
			std::cout << "\nmemory: ";
			for (const auto& res : results) {
				std::cout << res.memory << ' ';
			}
		}
		std::cout << "\n\n";
	}
	return results;
}

template <class TMatrix>
typename tests::SolverTester<TMatrix>::TestResult tests::SolverTester<TMatrix>::test(size_t index) {
	TestResult res;
	res.index = index;
	try {
		auto start = std::chrono::high_resolution_clock::now();
		{
			res.x = m_solver->solve(m_slaes[index].A, m_slaes[index].b);
		}
		/*auto itc = dynamic_cast<math::linal::IKrylovTypeLinearSystemSolver*>(m_solver.get())->get_last_solution_stats().iteration_count;
		auto acc = dynamic_cast<math::linal::IKrylovTypeLinearSystemSolver*>(m_solver.get())->get_last_solution_stats().accuracy;*/
		res.duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	}
	catch (std::exception& ex) {
		std::cout << ex.what() << '\n';
		return res;
	}
	catch (...) {
		std::cout << "Unknown exception\n";
		return res;
	}
	res.success = true;
	return res;
}

template <class TMatrix>
void tests::SolverTester<TMatrix>::print_test_result(const TestResult& res) const {
	if (m_params.print_slae) {
		//std::cout << "A = \n";
		//print_matrix(m_slaes[res.index].A);

		//std::cout << "\nb = ";
		//print_vector(m_slaes[res.index].b);

		if (m_slaes[res.index].solution.size()) {
			std::cout << "\n\ns = ";
			print_vector(m_slaes[res.index].solution);
		}

		std::cout << "\nx = ";
		print_vector(res.x);

		if (m_slaes[res.index].solution.size()) {
			std::cout << "\n\nIs equal: ";
			print_bool(res.x == m_slaes[res.index].solution);
		}
	}

	std::cout << "\n\nsuccess: ";
	print_bool(res.success);
	std::cout << "\n";

	if (m_params.print_time) {
		std::cout << "\ntime: " << res.duration.count() / 1000.f << " ms";
	}
	if (m_params.print_operation_count) {
		std::cout << "\noperations: " << res.operation;
	}
	if (m_params.print_memory) {
		std::cout << "\n\nmemory: " << res.memory << " kB";
	}
}

template <class TMatrix>
bool tests::SolverTester<TMatrix>::only_one_print() const {
	return static_cast<int>(m_params.print_time) + static_cast<int>(m_params.print_operation_count) + static_cast<int>(m_params.print_memory) == 1;
}

template <class ...Args>
static std::unique_ptr<math::linal::IPreconditioner> tests::SolverFactory::get_predictioner(PredictionerName name, Args&&... args) {
	switch (name) {
	case SolverFactory::PredictionerName::JACOBI:
		return std::make_unique<math::linal::JacobiPreconditioner>(std::forward<Args>(args)...);
	case SolverFactory::PredictionerName::ILU:
		return std::make_unique<math::linal::ILUPreconditioner>(std::forward<Args>(args)...);
	case SolverFactory::PredictionerName::MG:
		return std::make_unique<math::linal::MGPreconditioner>(std::forward<Args>(args)...);
	case SolverFactory::PredictionerName::None:
	default:
		return nullptr;
	}
}

template <class SolverT, class ...Args>
static std::unique_ptr<math::linal::ILinearSystemSolver> tests::SolverFactory::make_solver(std::unique_ptr<math::linal::IPreconditioner> preconditioner, Args&&... args) {
	if (preconditioner) {
		return std::make_unique<SolverT>(std::forward<Args>(args)..., std::move(preconditioner));
	}
	return std::make_unique<SolverT>(std::forward<Args>(args)...);
}

template <class TMatrix>
std::shared_ptr<tests::SolverTester<TMatrix>> tests::create_solver_tester(const std::vector<SLAE<TMatrix>>& slaes, SolverFactory::SolverName name, SolverFactory::PredictionerName pred, typename SolverTester<TMatrix>::Params const& params) {
	static SolverFactory factory;
	return std::make_shared<SolverTester<TMatrix>>(slaes, factory.get_solver(name, pred), factory.get_full_solver_name(name, pred), params);
}

namespace tests {
	
	template <>
	class TestGenerator<math::linal::BandMatrix> : public ITestGenerator {
	public:
		std::vector<SLAE<math::linal::BandMatrix>> get_SLAEs(const std::vector<size_t>& sizes, double sparsity = 0.5, bool async = false) {
			std::vector<SLAE<math::linal::BandMatrix>> slaes;
			if (async) {
				std::vector<std::future<SLAE<math::linal::BandMatrix>>> futures;

				for (const auto& size : sizes) {
					futures.emplace_back(std::async(std::launch::async, [this, size, sparsity]() {
						return get_SLAE(size, 0, sparsity);
						}));
				}

				std::vector<SLAE<math::linal::BandMatrix>> slaes;
				slaes.reserve(sizes.size());
				for (auto& fut : futures) {
					slaes.emplace_back(std::move(fut.get()));
				}
			}
			else {
				for (const auto& size : sizes) {
					slaes.emplace_back(get_SLAE(size, 0, sparsity));
				}
			}
			return slaes;
		}

		SLAE<math::linal::BandMatrix> get_SLAE(size_t n, size_t isl = 0, double sparsity = 0.5) {
			if (isl == 0) {
				isl = n / 3;
			}
			SLAE<math::linal::BandMatrix> slae;
			for (int i = 0; i < n; ++i) {
				slae.A.data().resize(n);
				for (auto& row : slae.A.data()) {
					row = get_random_vector(isl, sparsity);
				}
			}
			sparsity = 0.0;
			slae.b = get_random_vector(n, sparsity);
			return slae;
		}
	};

	template <>
	class TestGenerator<math::linal::DenseMatrix> : public ITestGenerator {
	public:
		std::vector<SLAE<math::linal::DenseMatrix>> get_SLAEs(const std::vector<size_t>& sizes, double sparsity = 0.5, bool async = false) {
			std::vector<SLAE<math::linal::DenseMatrix>> slaes;
			for (const auto& size : sizes) {
				slaes.emplace_back(get_SLAE(size, sparsity));
			}
			return slaes;
		}

		SLAE<math::linal::DenseMatrix> get_SLAE(size_t n, double sparsity = 0.5) {
			SLAE<math::linal::DenseMatrix> slae;
			for (int i = 0; i < n; ++i) {
				slae.A.data().resize(n);
				for (auto& row : slae.A.data()) {
					row = get_random_vector(n, sparsity);
				}
			}
			sparsity = 0.0;
			slae.b = get_random_vector(n, sparsity);
			return slae;
		}
	};

	template <>
	class TestGenerator<math::linal::SparseMatrix> : public ITestGenerator {
	public:
		std::vector<SLAE<math::linal::SparseMatrix>> get_SLAEs(const std::vector<size_t>& sizes, double sparsity = 0.5, bool async = false) {
			std::vector<SLAE<math::linal::SparseMatrix>> slaes;
			for (const auto& size : sizes) {
				slaes.emplace_back(get_SLAE(size, sparsity));
			}
			return slaes;
		}

		SLAE<math::linal::SparseMatrix> get_SLAE(size_t n, double sparsity = 0.5) {
			SLAE<math::linal::SparseMatrix> slae;
			// ...
			sparsity = 0.0;
			slae.b = get_random_vector(n, sparsity);
			return slae;
		}
	};

} // namespace tests

template <class TMatrix>
void tests::TestRunner<TMatrix>::run(const std::string& name, const std::vector<std::shared_ptr<SolverTester<TMatrix>>>& testers, const Params& params) {
	std::unordered_map<std::string, std::vector<typename SolverTester<TMatrix>::TestResult>> results;
	if (params.async) {
		results = async_run(name, testers, params);
	}
	else {
		for (const auto& tester : testers) {
			std::string name = tester->get_solver_name();
			std::vector<typename SolverTester<TMatrix>::TestResult> res;
			try {
				res = tester->run();
			}
			catch (...) {
				std::cout << "Exception in: " << name << "\n";
			}
			results.emplace(name, std::move(res));
		}
	}

	switch (params.save_type) {
	case SaveType::All:
		save_all_tests(name, results);
		break;
	case SaveType::OnlyTime:
		save_only_time(name, results);
		break;
	default:
		break;
	}
}

template <class TMatrix>
auto tests::TestRunner<TMatrix>::async_run(const std::string& name, const std::vector<std::shared_ptr<SolverTester<TMatrix>>>& testers, const Params& params) {
	std::unordered_map<std::string, std::future<std::vector<typename SolverTester<TMatrix>::TestResult>>> future_results;
	for (const auto& tester : testers) {
		std::string name = tester->get_solver_name();
		future_results.emplace(
			name,
			std::async(std::launch::async, [tester]() { return tester->run(); })
		);
	}

	std::unordered_map<std::string, std::vector<typename SolverTester<TMatrix>::TestResult>> results;
	for (auto& [name, fut] : future_results) {
		results.emplace(name, std::move(fut.get()));
	}

	return results;
}

template <class TMatrix>
void tests::TestRunner<TMatrix>::save_all_tests(const std::string& name, const std::unordered_map<std::string, std::vector<typename SolverTester<TMatrix>::TestResult>>& results) {
	std::string filename = name + ".txt";
	std::ofstream file(filename);
	if (!file.is_open()) {
		throw std::runtime_error("Failed to open file: " + filename);
	}

	for (const auto& [solver_name, test_results] : results) {
		file << "Solver: " << solver_name << "\n";
		file << "time: \n";

		for (const auto& result : test_results) {
			file << "Test - " << result.index
				<< ": success=" << result.success
				<< ", time=" << result.duration.count() << "µs"
				<< ", operations=" << result.operation
				<< ", memory=" << result.memory << "\n";
		}
		file << "\n\n";
	}

	file.close();
}

template <class TMatrix>
void tests::TestRunner<TMatrix>::save_only_time(const std::string& name, const std::unordered_map<std::string, std::vector<typename SolverTester<TMatrix>::TestResult>>& results) {
	std::string filename = name + ".txt";
	std::ofstream file(filename);
	if (!file.is_open()) {
		throw std::runtime_error("Failed to open file: " + filename);
	}

	for (const auto& [solver_name, test_results] : results) {
		file << "Solver: " << solver_name << "\n";
		file << "time: \n";

		for (const auto& result : test_results) {
			file << result.duration.count() / 1000.f << ", \t";
		}
		file << "\n\n";
	}

	file.close();
}
