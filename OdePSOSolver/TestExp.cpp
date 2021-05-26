#include "TestExp.h"

/// <summary>
/// Test Problem
/// </summary>
/// <param name="initalConditionIn"></param>
/// <param name="outStateIn"></param>
/// <param name="currentTime"></param>
void TestExp::getFunctionVector(const Eigen::VectorXd& initalConditionIn, Eigen::VectorXd& outStateIn, const double& currentTime)
{
	outStateIn(0) = -1*initalConditionIn(0);// (initalConditionIn(0) * initalConditionIn(0));
}
