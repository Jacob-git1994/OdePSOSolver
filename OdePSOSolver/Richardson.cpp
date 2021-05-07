#include "Richardson.h"

/// <summary>
/// Build our tables based off table size
/// </summary>
void odepso::Richardson::buildTables()
{
	//Ensure we don't try to allocate bad data
	if (verifyParameters())
	{
		throw std::logic_error("Invalid Parameters");
	}
	else
	{
		//Resize the rows
		richardsonTables.resize(tableSize);

		//Resize result vector and initalize everything to 0
		bestResult.resize(vectorSize);
		bestResult.fill(0.0);

		//Resize the cols
		for (std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>>>::iterator rowItr =
			richardsonTables.begin(); rowItr != richardsonTables.end(); ++rowItr)
		{
			//Resize
			rowItr->resize(tableSize);

			//Update the vectors with their new state and initalize them
			for (std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>>::iterator colItr =
				rowItr->begin(); colItr != rowItr->end(); ++colItr)
			{
				//Resize the vector and fill
				colItr->resize(vectorSize);
				colItr->fill(0.0);
			}
		}
	}
}

/// <summary>
/// Verify our parameters are valid before we start anything
/// </summary>
/// <returns></returns>
const bool odepso::Richardson::verifyParameters() const
{
	return (tableSize < 2) || (vectorSize < 1) || (reductionFactor < 2) || (dt <= 0.) || methodTruncationOrder < 1.5;
}

/// <summary>
/// Contructor
/// </summary>
odepso::Richardson::Richardson() : 
	reductionFactor(0.0), 
	error(0.0), 
	c(0.0), 
	dt(0.0),
	methodTruncationOrder(0.0),
	tableSize(0),
	vectorSize(0)
{
	//Nothing else to do here
}

/// <summary>
/// Copy constructor
/// </summary>
/// <param name="richardsonIn"></param>
odepso::Richardson::Richardson(const Richardson& richardsonIn) : 
	reductionFactor(richardsonIn.reductionFactor), 
	error(richardsonIn.error), 
	c(richardsonIn.c), 
	dt(richardsonIn.dt),
	methodTruncationOrder(richardsonIn.methodTruncationOrder),
	tableSize(richardsonIn.tableSize),
	vectorSize(richardsonIn.vectorSize),
	richardsonTables(richardsonIn.richardsonTables),
	bestResult(richardsonIn.bestResult)
{
	//Nothing else to do here
}

//Assign Constructor
odepso::Richardson& odepso::Richardson::operator=(const Richardson& richardsonIn)
{
	//Copy over the tables
	reductionFactor = richardsonIn.reductionFactor;
	error = richardsonIn.error;
	c = richardsonIn.c;
	dt = richardsonIn.dt;
	methodTruncationOrder = richardsonIn.methodTruncationOrder;
	tableSize = richardsonIn.tableSize;
	vectorSize = richardsonIn.vectorSize;
	richardsonTables = richardsonIn.richardsonTables;
	bestResult = richardsonIn.bestResult;
	
	//Return this
	return *this;
}

/// <summary>
/// Reset our tables and reduction factor
/// </summary>
/// <param name="tableSizeIn"></param>
/// <param name="vectorSizeIn"></param>
/// <param name="reductionFactorIn"></param>
/// <param name="truncationOrderIn"></param>
/// <param name="dtIn"></param>
void odepso::Richardson::reset(const size_t& tableSizeIn, const size_t& vectorSizeIn, const long double& reductionFactorIn, const long double& dtIn, const long double& truncationOrderIn)
{
	//Set our new parameters
	reductionFactor = reductionFactorIn;
	tableSize = tableSizeIn;
	vectorSize = vectorSizeIn;
	dt = dtIn;
	methodTruncationOrder = truncationOrderIn;

	//Build our tables. Fail tables if memory can't be alllocated
	try
	{
		buildTables();
	}
	catch (std::exception& e)
	{
		std::cerr << e.what();
		exit(1);
	}
}

/// <summary>
/// Add a copy vector to our table
/// </summary>
/// <param name="rowColIn"></param>
/// <param name="solutionIn"></param>
/// <returns></returns>
odepso::Richardson& odepso::Richardson::append(const size_t& rowColIn, const Eigen::VectorXd& solutionIn)
{
	//Copy solution into table
	richardsonTables.at(rowColIn).at(0) = solutionIn;

	//Return this
	return *this;
}

/// <summary>
/// Add a copy vector to any position in our tables
/// </summary>
/// <param name="rowIn"></param>
/// <param name="colIn"></param>
/// <param name="solutionIn"></param>
/// <returns></returns>
odepso::Richardson& odepso::Richardson::append(const size_t& rowIn, const size_t colIn, const Eigen::VectorXd& solutionIn)
{
	//Copy solution into table
	richardsonTables.at(rowIn).at(colIn) = solutionIn;

	//Return this
	return *this;
}

/// <summary>
/// Move the solution vector into our table
/// </summary>
/// <param name="rowColIn"></param>
/// <param name="solutionIn"></param>
/// <returns></returns>
odepso::Richardson& odepso::Richardson::append(const size_t& rowColIn, Eigen::VectorXd&& solutionIn)
{
	//Copy solution into table
	richardsonTables.at(rowColIn).at(0) = std::move(solutionIn);

	//Return this
	return *this;
}

/// <summary>
/// Move the solution anywhere into our table
/// </summary>
/// <param name="rowIn"></param>
/// <param name="colIn"></param>
/// <param name="solutionIn"></param>
/// <returns></returns>
odepso::Richardson& odepso::Richardson::append(const size_t& rowIn, const size_t colIn, Eigen::VectorXd&& solutionIn)
{
	//Copy solution into table
	richardsonTables.at(rowIn).at(colIn) = std::move(solutionIn);

	//Return this
	return *this;
}

/// <summary>
/// Get the reduction factor
/// </summary>
/// <returns></returns>
const double& odepso::Richardson::getReductionFactor() const
{
	return reductionFactor;
}

/// <summary>
/// Get the error
/// </summary>
/// <returns></returns>
const double& odepso::Richardson::getError() const
{
	return error;
}

/// <summary>
/// Get the convergence estimate
/// </summary>
/// <returns></returns>
const double& odepso::Richardson::getC() const
{
	return c;
}

/// <summary>
/// Getter to the best result
/// </summary>
/// <returns></returns>
const Eigen::VectorXd& odepso::Richardson::getBestResult() const
{
	return bestResult;
}

/// <summary>
/// Getter to the table size
/// </summary>
/// <returns></returns>
const size_t& odepso::Richardson::getTableSize() const
{
	return tableSize;
}

/// <summary>
/// Getter to the trunction order
/// </summary>
/// <returns></returns>
const double& odepso::Richardson::getTruncationOrder() const
{
	return methodTruncationOrder;
}

/// <summary>
/// Calculate the all the parameters used for this extrapolation
/// </summary>
/// <returns></returns>
void odepso::Richardson::calculateRichardsonExtrapolation()
{
	//Our extrapolation coeifficents
	double alphaOne = 0.0;
	double alphaTwo = 0.0;

	//Check if everything is valid before running
	if (verifyParameters())
	{
		throw std::logic_error("Parameters not initalize correctly");
	}
	else
	{
		for (size_t i = 0; i < richardsonTables.size() - 1; ++i)
		{
			for (size_t j = 0; j <= i; ++j)
			{
				//Set our coefficents. methodTruncError - theoretical error. i is the reduction done for that row. j is assumed that the next method has an order of theretical error + 1
				alphaOne = std::pow(reductionFactor, (methodTruncationOrder * (static_cast<long double>(i) + 1.0)) + static_cast<long double>(j));
				alphaTwo = std::pow(reductionFactor, (methodTruncationOrder * (static_cast<long double>(i))) + static_cast<long double>(j));

				//Get a new extrapolated state
				Eigen::VectorXd extrapResult = odepso::Common::smartAdd({ alphaOne * richardsonTables.at(i + 1).at(j),-alphaTwo * richardsonTables.at(i).at(j) }) / (alphaOne - alphaTwo);

				//Add the new result into our table
				append(i + 1, j + 1, std::move(extrapResult));
			}
		}

		//Best results are at the bottom right corner
		bestResult = richardsonTables.back().back();

		//Calculate the error
		error = (richardsonTables.at(tableSize - 1).at(tableSize - 1) -
			richardsonTables.at(tableSize - 2).at(tableSize - 2)).cwiseAbs().maxCoeff();

		//Calculate the convergence
		c = std::abs(std::log(error) / std::log(dt));
	}
}
