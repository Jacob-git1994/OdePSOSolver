#include "RungeKutta2.h"

/// <summary>
/// Copy constructor
/// </summary>
/// <param name="rungeKutta2In"></param>
odepso::RungeKutta2::RungeKutta2(const RungeKutta2& rungeKutta2In) : SolverIF(rungeKutta2In)
{
	//Copy Solver Vectors
	k1 = rungeKutta2In.k1;
	k2 = rungeKutta2In.k2;
}

/// <summary>
/// Assign constructor
/// </summary>
/// <param name="rungeKutta2In"></param>
/// <returns></returns>
odepso::RungeKutta2& odepso::RungeKutta2::operator=(const RungeKutta2& rungeKutta2In)
{
	//Assign to base class
	SolverIF::operator=(rungeKutta2In);

	//Copy the solver vectors
	k1 = rungeKutta2In.k1;
	k2 = rungeKutta2In.k2;

	//Return this
	return *this;
}

/// <summary>
/// Initalize the solver vectors
/// </summary>
/// <param name="vecSizeIn"></param>
void odepso::RungeKutta2::initalizeSolverVectors(const size_t& vecSizeIn)
{
	//Build our list of solver vectors
	std::list<Eigen::VectorXd*> solverHelperVectors = { &k1, &k2 };

	//Initalize our vectors
	initalizeVectors(solverHelperVectors, vecSizeIn);
}

/// <summary>
/// Propagrate the derivative forward in time n times.
/// </summary>
/// <param name="initalConditionsIn"></param>
/// <param name="updateStateIn"></param>
/// <param name="dtIn"></param>
/// <param name="numberOfStepsIn"></param>
/// <param name="currentTimeIn"></param>
void odepso::RungeKutta2::updateTimeStep(
	const Eigen::VectorXd& initalConditionsIn,
	Eigen::VectorXd& updateStateIn,
	const double& dtIn,
	const unsigned int& numberOfStepsIn,
	const double& currentTimeIn)
{
	//Initalize our vectors
	updateStateIn = initalConditionsIn;

	//Propagate forward number of steps in time
	for (unsigned int i = 0; i < numberOfStepsIn; ++i)
	{
		//Get the derivative at multiple steps
		problem->getFunctionVector(updateStateIn, k1, currentTimeIn);
		problem->getFunctionVector(odepso::Common::smartAdd({ updateStateIn, .5 * dtIn * k1 }), k2, currentTimeIn + .5 * dtIn);

		//Update the state
		updateStateIn += dtIn * k2;
	}
}

/// <summary>
/// Propagrate the derivative forward in time by dt one time
/// </summary>
/// <param name="initalConditionsIn"></param>
/// <param name="updateStateIn"></param>
/// <param name="dtIn"></param>
/// <param name="currentTimeIn"></param>
void odepso::RungeKutta2::updateTimeStep(
	const Eigen::VectorXd& initalConditionsIn,
	Eigen::VectorXd& updateStateIn,
	const double& dtIn,
	const double& currentTimeIn)
{
	//Initalize our vectors
	updateStateIn = initalConditionsIn;

	//Get the derivative at multiple steps
	problem->getFunctionVector(updateStateIn, k1, currentTimeIn);
	problem->getFunctionVector(updateStateIn + .5 * dtIn * k1, k2, currentTimeIn + .5 * dtIn);

	//Update the state
	updateStateIn += dtIn * (k1 + k2);
}

/// <summary>
/// Get the theoretical convergence
/// </summary>
/// <returns></returns>
const double odepso::RungeKutta2::getTheoreticalConvergence() const
{
	return 3.0;
}

/// <summary>
/// Clone our current method
/// </summary>
/// <returns></returns>
odepso::SolverIF* odepso::RungeKutta2::deepClone() const
{
	//Create a new clone
	odepso::SolverIF* tempSolverPtr = new RungeKutta2(*this);

	//Make sure this is valid
	if (!tempSolverPtr)
	{
		throw std::bad_alloc();
	}
	else
	{
		return tempSolverPtr;
	}
}
