#pragma once

#include <Eigen/Dense>
#include <utility>

#include "Common.h"

namespace odepso
{
	class ProblemWrapperIF
	{
	public:

		//Get the y' = f(X)
		virtual void getFunctionVector(const Eigen::VectorXd& initalConditionIn, Eigen::VectorXd& outStateIn, const double& currentTime) = 0;

	};
}
