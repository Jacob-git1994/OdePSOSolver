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

		//Smart Adder for Vectors
		static Eigen::VectorXd smartAdd(const std::list<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>>& vectorArguments);

		//Compare our solvers
		struct CompSolvers
		{
			inline bool operator()(SOLVER_TYPES solver1, SOLVER_TYPES solver2) const
			{
				return static_cast<unsigned int>(solver1) < static_cast<unsigned int>(solver2);
			}
		};

		static CompSolvers compSolvers;
	};
}