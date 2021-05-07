#pragma once

#include <Eigen/Dense>
#include <list>

namespace odepso
{
	//Common class to hold common enumerations or states
	class Common
	{
	public:

		//Types of solvers
		enum class SOLVER_TYPES
		{
			EULER = 10,
			RK2 = 20,
			RK4 = 30,
			IMPLICT_EULER = 40,
			CRANK_NICOLSON = 50
		};

		//Solver Comparison
		static const bool solverComp(const SOLVER_TYPES& solverOne, const SOLVER_TYPES& solverTwo);

		//Smart Adder for Vectors
		static Eigen::VectorXd smartAdd(const std::list<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>>& vectorArguments);
	};
}