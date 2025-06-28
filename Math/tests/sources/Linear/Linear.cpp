#include "Linear.h"

#include <iostream>

using namespace math;

void tests::print_vector(const linal::DVector& vec) {
	std::cout << "[ ";
	if (!vec.empty()) {
		for (size_t i = 0; i < vec.size() - 1; ++i) {
			std::cout << vec[i] << ", ";
		}
		std::cout << vec.back();
	}
	std::cout << " ]";
}

std::unique_ptr<linal::ILinearSystemSolver> tests::SolverFactory::get_solver(SolverName name, PredictionerName pred) {
	using namespace linal;
	using ItParams = IIterativeLinearSystemSolver::IterativeSolvingParams;
	using SolverParams = ILinearSystemSolver::Params;
	auto predictioner = get_predictioner(pred);
	switch (name) {
	case SolverName::BiCGSTAB:
		return make_solver<BiCGSTABLinearSystemSolver>(std::move(predictioner), 30, ItParams{ TOLERANCE, 7000 }, SolverParams{ });
	case SolverName::CG:
		return make_solver<CGLinearSystemSolver>(std::move(predictioner), 30, ItParams{ TOLERANCE, 7000 }, SolverParams{ });
	case SolverName::FGMRES:
		return make_solver<FGMRESLinearSystemSolver>(std::move(predictioner), 30, ItParams{ TOLERANCE, 7000 }, SolverParams{ });
	case SolverName::GMRES:
		return make_solver<GMRESLinearSystemSolver>(std::move(predictioner), 30, ItParams{ TOLERANCE, 7000 }, SolverParams{ });
	case SolverName::MINRES:
		return make_solver<MINRESLinearSystemSolver>(std::move(predictioner), 30, ItParams{ TOLERANCE, 7000 }, SolverParams{ });
	case SolverName::SOR:
		return make_solver<SORLinearSystemSolver>(std::move(predictioner), 1.0, 0, ItParams{ TOLERANCE, 7000 }, SolverParams{ });
	case SolverName::Gauss:
	default:
		return std::make_unique<GaussLinearSystemSolver>();
	}
}

std::string tests::SolverFactory::get_full_solver_name(SolverName name, PredictionerName pred) {
	auto solver_name = get_solver_name(name);
	if (pred != PredictionerName::None) {
		solver_name += "+" + get_predictioner_name(pred);
	}
	return solver_name;
}

std::string tests::SolverFactory::get_solver_name(SolverName name) {
	if (auto it = s_solver_names.find(name); it != s_solver_names.end()) {
		return it->second;
	}
	return {};
}

std::string tests::SolverFactory::get_predictioner_name(PredictionerName name) {
	if (auto it = s_predictioner_names.find(name); it != s_predictioner_names.end()) {
		return it->second;
	}
	return get_predictioner_name(PredictionerName::None);
}

std::unordered_map<tests::SolverFactory::SolverName, std::string> tests::SolverFactory::s_solver_names = {
	{SolverName::BiCGSTAB, "BiCGSTAB"},
	{SolverName::CG, "CG"},
	{SolverName::FGMRES, "FGMRES"},
	{SolverName::Gauss, "Gauss"},
	{SolverName::GMRES, "GMRES"},
	{SolverName::MINRES, "MINRES"},
	{SolverName::SOR, "SOR"},
};

std::unordered_map<tests::SolverFactory::PredictionerName, std::string> tests::SolverFactory::s_predictioner_names = {
	{PredictionerName::None, "None"},
	{PredictionerName::JACOBI, "JACOBI"},
	{PredictionerName::ILU, "ILU"},
	{PredictionerName::MG, "MG"},
};

tests::ITestGenerator::ITestGenerator()
	: gen(rd())
	, sparsity_dist(0.0, 1.0)
{
}

double tests::ITestGenerator::get_random_ranged_double(double lower_bound, double upper_bound, bool can_be_zero) {
	std::uniform_real_distribution<double> dist(lower_bound, upper_bound);

	double result;
	do {
		result = dist(gen);
	} while (!can_be_zero && dcmp(result) == 0);

	return result;
}

double tests::ITestGenerator::get_random_double() {
	static std::uniform_real_distribution<double> dist(-1000., 1000.);
	return dist(gen);
}

linal::DVector tests::ITestGenerator::get_random_vector(size_t n, const double& sparsity) {
	linal::DVector vector(n, 0.0);
	for (auto& el : vector) {
		if (sparsity > 0 || sparsity_dist(gen) > sparsity) {
			el = get_random_double();
		}
	}
	vector.at(0) = get_random_ranged_double(-1000., 1000., false);
	return vector;
}

template <>
std::vector<tests::SLAE<math::linal::BandMatrix>> tests::get_base_SLAEs() {
	using BandSLAE = SLAE<math::linal::BandMatrix>;
	std::vector<BandSLAE> res;
	{ // 0
		BandSLAE axb;
		axb.A = { {1}, {1}, {1}, {1}, {1} };
		axb.b = { 1, 2, 3, 4, 5 };
		axb.solution = { 1, 2, 3, 4, 5 };
		res.emplace_back(std::move(axb));
	}
	{ // 1
		BandSLAE axb;
		axb.A = { {5, 1}, {5, 1}, {5, 1}, {5, 1}, {5, 0} };
		axb.b = { 1, 2, 3, 4, 5 };
		axb.solution = { 0.143182, 0.284091, 0.436364, 0.534091, 0.893182 };
		res.emplace_back(std::move(axb));
	}
	{ // 2
		BandSLAE axb;
		axb.A = {
			{2, 1, 1, 1, 1},
			{3, 1, 1, 1, 0},
			{4, 1, 1, 0, 0},
			{5, 1, 1, 0, 0},
			{6, 0, 0, 0, 0} };
		axb.b = { 6, 7, 8, 9, 10 };
		axb.solution = { 1, 1, 1, 1, 1 };
		res.emplace_back(std::move(axb));
	}
	{ // 3
		BandSLAE axb;
		axb.A = {
			{3, 2, 1},
			{3, 2, 1},
			{3, 2, 1},
			{3, 2, 0},
			{3, 0, 0} };
		axb.b = { 6, 8, 10, 8, 6 };
		axb.solution = { 1, 0.5, 2, 0.5, 1 };
		res.emplace_back(std::move(axb));
	}
	{ // 4
		BandSLAE axb;
		axb.A = {
			{6, -1, 1, -1, 1},
			{6, -1, 1, -1, 0},
			{6, -1, 1, 0, 0},
			{6, -1, 0, 0, 0},
			{6,  0, 0, 0, 0} };
		axb.b = { 6, -6, 6, -6, 6 };
		axb.solution = { 0.6, -0.6,  0.6, -0.6,  0.6 };
		res.emplace_back(std::move(axb));
	}
	return res;
}

template <>
std::vector<tests::SLAE<math::linal::DenseMatrix>> tests::get_base_SLAEs() {
	using FSLAE = SLAE<math::linal::DenseMatrix>;
	std::vector<FSLAE> res;
	{ // 0
		FSLAE axb;
		axb.A = {
			{1, 0, 0, 0, 0},
			{0, 1, 0, 0, 0},
			{0, 0, 1, 0, 0},
			{0, 0, 0, 1, 0},
			{0, 0, 0, 0, 1} };
		axb.b = { 1, 2, 3, 4, 5 };
		axb.solution = { 1, 2, 3, 4, 5 };
		res.emplace_back(std::move(axb));
	}
	{ // 1
		FSLAE axb;
		axb.A = {
			{5, 1, 0, 0, 0},
			{1, 5, 1, 0, 0},
			{0, 1, 5, 1, 0},
			{0, 0, 1, 5, 1},
			{0, 0, 0, 1, 5} };
		axb.b = { 1, 2, 3, 4, 5 };
		axb.solution = { 0.143182, 0.284091, 0.436364, 0.534091, 0.893182 };
		res.emplace_back(std::move(axb));
	}
	{ // 2
		FSLAE axb;
		axb.A = {
			{2, 1, 1, 1, 1},
			{1, 3, 1, 1, 1},
			{1, 1, 4, 1, 1},
			{1, 1, 1, 5, 1},
			{1, 1, 1, 1, 6} };
		axb.b = { 6, 7, 8, 9, 10 };
		axb.solution = { 1, 1, 1, 1, 1 };
		res.emplace_back(std::move(axb));
	}
	{ // 3
		FSLAE axb;
		axb.A = {
			{3, 2, 1, 0, 0},
			{2, 3, 2, 1, 0},
			{1, 2, 3, 2, 1},
			{0, 1, 2, 3, 2},
			{0, 0, 1, 2, 3} };
		axb.b = { 6, 8, 10, 8, 6 };
		axb.solution = { 1, 0.5, 2, 0.5, 1 };
		res.emplace_back(std::move(axb));
	}
	{ // 4
		FSLAE axb;
		axb.A = {
			{6, -1, 1, -1, 1},
			{-1, 6, -1, 1, -1},
			{1, -1, 6, -1, 1},
			{-1, 1, -1, 6, -1},
			{1, -1, 1, -1, 6} };
		axb.b = { 6, -6, 6, -6, 6 };
		axb.solution = { 0.6, -0.6,  0.6, -0.6,  0.6 };
		res.emplace_back(std::move(axb));
	}
	{ // 5
		FSLAE axb;
		axb.A = {
			{1, 0, 0, 0, 0},
			{0, 1, 0, 0, 0},
			{0, 0, 1, 0, 0},
			{0, 0, 0, 1, 0},
			{0, 0, 0, 0, 1} };
		axb.b = { 1, 2, 3, 4, 5 };
		axb.solution = { 1, 2, 3, 4, 5 };
		//res.emplace_back(std::move(axb));
	}
	return res;
}

template <>
std::vector<tests::SLAE<math::linal::SparseMatrix>> tests::get_base_SLAEs() {
	using SSLAE = SLAE<math::linal::SparseMatrix>;
	std::vector<SSLAE> res;
	{ // 0
		SSLAE axb;
		axb.A = {};
		axb.b = { 1, 2, 3, 4, 5 };
		axb.solution = { 1, 2, 3, 4, 5 };
		//res.emplace_back(std::move(axb));
	}
	return res;
}