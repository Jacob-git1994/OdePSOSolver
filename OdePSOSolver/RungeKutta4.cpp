#include "RungeKutta4.h"

/// <summary>
/// Copy constructor. Copies the vectors and problem
/// </summary>
/// <param name="rungeKutta4In"></param>
odepso::RungeKutta4::RungeKutta4(const RungeKutta4& rungeKutta4In) : odepso::SolverIF(rungeKutta4In)
{
	//Copy over the solver vectors
	k1 = rungeKutta4In.k1;
	k2 = rungeKutta4In.k2;
	k3 = rungeKutta4In.k3;
	k4 = rungeKutta4In.k4;
}

/// <summary>
/// Assign constructor. Copy the problem to the base class and the vectors to rk4 class.
/// </summary>
/// <param name="rungeKutta4In"></param>
/// <returns></returns>
odepso::RungeKutta4& odepso::RungeKutta4::operator=(const RungeKutta4& rungeKutta4In)
{
	//Copy the problem
	odepso::SolverIF::operator=(rungeKutta4In);

	//Copy the solver vectors
	k1 = rungeKutta4In.k1;
	k2 = rungeKutta4In.k2;
	k3 = rungeKutta4In.k3;
	k4 = rungeKutta4In.k4;

	//Return this
	return *this;
}

/// <summary>
/// Initalize the vectors sizes
/// </summary>
/// <param name="vecSizeIn"></param>
void odepso::RungeKutta4::initalizeSolverVectors(const size_t& vecSizeIn)
{
	//Get pointers to out solver vectors
	std::list<Eigen::VectorXd*> solverVectors = { &k1,&k2,&k3,&k4 };

	//Allocate
	initalizeVectors(solverVectors, vecSizeIn);
}

/// <summary>
/// Propagrate the derivative forward in time
/// </summary>
/// <param name="initalConditionsIn"></param>
/// <param name="updateStateIn"></param>
/// <param name="dtIn"></param>
/// <param name="numberOfStepsIn"></param>
/// <param name="currentTimeIn"></param>
void odepso::RungeKutta4::updateTimeStep(
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
		problem->getFunctionVector(updateStateIn + .5 * dtIn * k1, k2, currentTimeIn + .5 * dtIn);
		problem->getFunctionVector(updateStateIn + .5 * dtIn * k2, k3, currentTimeIn + .5 * dtIn);
		problem->getFunctionVector(updateStateIn + dtIn * k3, k4, currentTimeIn + dtIn);

		//Update the state
		updateStateIn += (dtIn / 6.0) * (k1 + 2. * k2 + 2. * k3 + k4);
	}
}

/// <summary>
/// Propagate the derivative forward ONCE in time. 
/// </summary>
/// <param name="initalConditionsIn"></param>
/// <param name="updateStateIn"></param>
/// <param name="dtIn"></param>
/// <param name="currentTimeIn"></param>
void odepso::RungeKutta4::updateTimeStep(
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
	problem->getFunctionVector(updateStateIn + .5 * dtIn * k2, k3, currentTimeIn + .5 * dtIn);
	problem->getFunctionVector(updateStateIn + dtIn * k3, k4, currentTimeIn + dtIn);

	//Update the state
	updateStateIn += (dtIn / 6.0) * odepso::Common::smartAdd({ k1, 2.0 * k2, 2.0 * k3, k4 }); //(k1 + 2. * k2 + 2. * k3 + k4);
}

const double odepso::RungeKutta4::getTheoreticalConvergence() const
{
	return 4.0;
}

/// <summary>
/// Clone our current method
/// </summary>
/// <returns></returns>
odepso::SolverIF* odepso::RungeKutta4::deepClone() const
{
	//Create a new clone
	odepso::SolverIF* tempSolverPtr = new RungeKutta4(*this);

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
