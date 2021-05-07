#pragma once
#include "ProblemWrapperIF.h"

using namespace odepso;

class TestExp :
    public ProblemWrapperIF
{
    virtual void getFunctionVector(const Eigen::VectorXd& initalConditionIn, Eigen::VectorXd& outStateIn, const double& currentTime) override;

};

