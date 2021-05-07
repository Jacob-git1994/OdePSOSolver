#include "Common.h"

/// <summary>
/// Comparison between two solvers
/// </summary>
/// <param name="solverOne"></param>
/// <param name="solverTwo"></param>
/// <returns></returns>
const bool odepso::Common::solverComp(const SOLVER_TYPES& solverOne, const SOLVER_TYPES& solverTwo)
{
	return static_cast<unsigned int>(solverOne) < static_cast<unsigned int>(solverTwo);
}

/// <summary>
/// Smart way to add vectors to correct for errors
/// </summary>
/// <param name="vectorArguments"></param>
/// <returns></returns>
Eigen::VectorXd odepso::Common::smartAdd(const std::list<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>>& vectorArguments)
{
	//Sum Vector Accumulator
	Eigen::VectorXd sum(vectorArguments.begin()->size());
	sum.fill(0.0);

	//Compensation Vector
	Eigen::VectorXd c(vectorArguments.begin()->size());
	c.fill(0.0);

	//Accumate all the vectors
	for (std::list<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>>::const_iterator vectorItr = vectorArguments.cbegin();
		vectorItr != vectorArguments.cend(); ++vectorItr)
	{
		//Correct the sum
		Eigen::VectorXd y = *(vectorItr) - c;

		//low order didgets are lost when accumulating
		Eigen::VectorXd t = sum + y;

		//Cancle the high-order part of y
		c = (t - sum) - y;

		//Update the sum
		sum = t;
	}

	//Return our new sum
	return sum;
}
