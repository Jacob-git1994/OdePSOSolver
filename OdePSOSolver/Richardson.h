#pragma once

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "Common.h"

namespace odepso
{
	//Richardson class to calculate the error and extrapolate to get a better answer
	class Richardson
	{
	private:

		//How much each solution was reduced by
		double reductionFactor;

		//Error found
		double error;

		//Convergence calculated
		double c;

		//Step Size
		double dt;

		//Method theoretical truncation order
		double methodTruncationOrder;

		//How many tables are required
		size_t tableSize;

		//Vector Size
		size_t vectorSize;

		//Richardson tables to store the results
		std::vector<std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>>> richardsonTables;

		//Best result
		Eigen::VectorXd bestResult;

		//Build our table vector
		void buildTables();

		//Verify our parameters
		const bool verifyParameters() const;

	public:

		//Constructor
		Richardson();

		//Copy Constructor
		Richardson(const Richardson& richardsonIn);

		//Assign constructor
		Richardson& operator=(const Richardson& richardsonIn);

		//Default deconstructor
		~Richardson() = default;

		//Reset the parameters for this table
		void reset(const size_t& tableSizeIn, const size_t& vectorSizeIn, const long double& reductionFactorIn, const long double& dtIn, const long double& truncationOrderIn);

		//Append to the tables
		Richardson& append(const size_t& rowColIn, const Eigen::VectorXd& solutionIn);

		//Append to the tables
		Richardson& append(const size_t& rowIn, const size_t colIn, const Eigen::VectorXd& solutionIn);

		//Append to the tables
		Richardson& append(const size_t& rowColIn, Eigen::VectorXd&& solutionIn);

		//Append to the tables
		Richardson& append(const size_t& rowIn, const size_t colIn, Eigen::VectorXd&& solutionIn);

		//ReductionFactor Getter
		const double& getReductionFactor() const;

		//Error Getter
		const double& getError() const;

		//Convergence Getter
		const double& getC() const;

		//Best result getter
		const Eigen::VectorXd& getBestResult() const;

		//Getter to table size
		const size_t& getTableSize() const;

		const double& getTruncationOrder() const;

		//Calculate richardson extrapolation on results
		void calculateRichardsonExtrapolation();
	};
}