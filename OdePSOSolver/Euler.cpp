#include "Euler.h"

/// <summary>
/// Copy Constructor. Copies the base class and our current solver vectors
/// </summary>
/// <param name="eulerIFIn"></param>
odepso::Euler::Euler(const Euler& eulerIn) : SolverIF(eulerIn)
{
	//Copy the vectors
	k1 = eulerIn.k1;
}

/// <summary>
/// Assign Constructor. This will copy over the vectors
/// </summary>
/// <param name="eulerIn"></param>
/// <returns></returns>
odepso::Euler& odepso::Euler::operator=(const Euler& eulerIn)
{
	//Assign constructor to the base
	odepso::SolverIF::operator=(eulerIn);

	//This will copy over the vectors
	k1 = eulerIn.k1;

	//Return this
	return *this;
}

/// <summary>
/// Initalize the vectors sizes
/// </summary>
/// <param name="vecSizeIn"></param>
void odepso::Euler::initalizeSolverVectors(const size_t& vecSizeIn)
{
	//Get pointers to out solver vectors
	std::list<Eigen::VectorXd*> solverVectors = { &k1};

	//Allocate
	initalizeVectors(solverVectors, vecSizeIn);
}

/// <summary>
/// Propagate the derivative forward in time.
/// </summary>
/// <param name="initalConditionsIn"></param>
/// <param name="updateStateIn"></param>
/// <param name="dtIn"></param>
/// <param name="numberOfStepsIn"></param>
/// <param name="currentTimeIn"></param>
void odepso::Euler::updateTimeStep(
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
		//Get the derivative
		problem->getFunctionVector(updateStateIn, k1, currentTimeIn);

		//Update the state
		updateStateIn += dtIn * k1;
	}
}

/// <summary>
/// This method will update the time step ONCE by dt
/// </summary>
/// <param name="initalConditionsIn"></param>
/// <param name="updateStateIn"></param>
/// <param name="dtIn"></param>
/// <param name="currentTimeIn"></param>
void odepso::Euler::updateTimeStep(
	const Eigen::VectorXd& initalConditionsIn,
	Eigen::VectorXd& updateStateIn,
	const double& dtIn,
	const double& currentTimeIn)
{
	//Initalize our vectors
	updateStateIn = initalConditionsIn;

	//Get the derivative
	problem->getFunctionVector(updateStateIn, k1, currentTimeIn);

	//Update the state
	updateStateIn += dtIn * k1;
}

/// <summary>
/// Return the theortical convergence.
/// </summary>
/// <returns></returns>
const double odepso::Euler::getTheoreticalConvergence() const
{
	return 2.;
}

/// <summary>
/// Clone our current method
/// </summary>
/// <returns></returns>
odepso::SolverIF* odepso::Euler::deepClone() const
{
	//Create a new clone
	odepso::SolverIF* tempSolverPtr = new Euler(*this);

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
