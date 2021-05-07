#include "SolverIF.h"
/// <summary>
/// Copy the infomration over
/// </summary>
/// <param name="SolverIFIn"></param>
odepso::SolverIF::SolverIF(const SolverIF& solverIFIn)
{
	problem = solverIFIn.problem;
}
/// <summary>
/// Assigment constructor
/// </summary>
/// <param name="solverIFIn"></param>
/// <returns></returns>
odepso::SolverIF& odepso::SolverIF::operator=(const SolverIF& solverIFIn)
{
	//Copy the problem over
	problem = solverIFIn.problem;

	//Return this
	return *this;
}

/// <summary>
/// Set our problem for this solver
/// </summary>
/// <param name="problemIn"></param>
void odepso::SolverIF::initalizeProblem(std::shared_ptr<ProblemWrapperIF> problemIn)
{
	//"copy" the shared pointer
	problem = problemIn;
}

/// <summary>
/// Initalize our solverHelperVector regardless of what method is being used.
/// </summary>
/// <param name="solverHelperVectorsIn"></param>
/// <param name="vecSizeIn"></param>
void odepso::SolverIF::initalizeVectors(std::list<Eigen::VectorXd*>& solverHelperVectorsIn, const size_t& vecSizeIn)
{
	//Make sure we input something non-zero
	if (vecSizeIn < 1)
	{
		throw std::invalid_argument("Invalid Size");
	}

	//Go through all the vectors
	for (std::list<Eigen::VectorXd*>::iterator listElementPtr = solverHelperVectorsIn.begin();
		listElementPtr != solverHelperVectorsIn.end(); ++listElementPtr)
	{
		//Get our solver vector
		Eigen::VectorXd& currentSolverVector = *(*listElementPtr);

		//Check if we need to resize the vector
		if (currentSolverVector.size() != vecSizeIn || static_cast<bool>(currentSolverVector.size()))
		{
			//Check if we have enough memory to resize/build
			try
			{
				//Resize the vectors
				currentSolverVector.resize(vecSizeIn);
			}
			catch (std::exception& e)
			{
				std::cerr << e.what();
				exit(1);
			}
		}
	}
}

