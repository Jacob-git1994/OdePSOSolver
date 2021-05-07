#include "OdeSolverParameters.h"

using namespace odepso;

const bool OdeSolverParameters::isUniqueAndOrdered(const double& lowerIn, const double& upperIn) const
{
    return (lowerIn > std::numeric_limits<double>::epsilon()) &&
        (upperIn > std::numeric_limits<double>::epsilon()) &&
        (std::fabs(lowerIn - upperIn) > std::numeric_limits<double>::epsilon()) &&
        (upperIn - lowerIn > 0.0);
}

/// <summary>
/// Initialize the parameters
/// </summary>
odepso::OdeSolverParameters::OdeSolverParameters() :
    currentError(0.0),
    totalError(0.0),
    minAllowedError(0.0),
    maxAllowedError(0.0),
    currentRunTime(0.0),
    totalRunTime(0.0),
    richardsonLevels(0),
    minRichardsonTables(0),
    maxRichardsonTables(0),
    c(0.0),
    dt(0.0),
    minDt(0.0),
    maxDt(0.0),
    allowedEuler(false),
    allowedRK2(false),
    allowedRK4(false),
    isSystemStiff(false),
    particleParams(odepso::ParticleParameters()),
    currentTime(0.0),
    reductionFactor(0.0)
{
    //Nothing else to do here
}

/// <summary>
/// Initalize our solver parameters
/// </summary>
/// <param name="currentErrorIn"></param>
/// <param name="totalErrorIn"></param>
/// <param name="minAllowedErrorIn"></param>
/// <param name="maxAllowedErrorIn"></param>
/// <param name="currentRunTimeIn"></param>
/// <param name="totalRunTimeIn"></param>
/// <param name="maxRunTimeIn"></param>
/// <param name="richardsonLevelsIn"></param>
/// <param name="minRichardsonTablesIn"></param>
/// <param name="maxRichardsonTablesIn"></param>
/// <param name="cIn"></param>
/// <param name="dtIn"></param>
/// <param name="minDtIn"></param>
/// <param name="maxDtIn"></param>
/// <param name="currentCostIn"></param>
/// <param name="totalCostIn"></param>
OdeSolverParameters::OdeSolverParameters(
    const double& minAllowedErrorIn, 
    const double& maxAllowedErrorIn, 
    const size_t& minRichardsonTablesIn,
    const size_t& maxRichardsonTablesIn,
    const double& minDtIn, 
    const double& maxDtIn, 
    const ParticleParameters& particleParametersIn,
    const double& reductionFactorIn) : 
    currentError(0.0),
    totalError(0.0),
    minAllowedError(minAllowedErrorIn),
    maxAllowedError(maxAllowedErrorIn),
    currentRunTime(0.0),
    totalRunTime(0.0),
    richardsonLevels(0),
    minRichardsonTables(minRichardsonTablesIn),
    maxRichardsonTables(maxRichardsonTablesIn),
    c(0.0),
    dt(0.0),
    minDt(minDtIn),
    maxDt(maxDtIn),
    allowedEuler(true),
    allowedRK2(true),
    allowedRK4(true),
    isSystemStiff(false),
    particleParams(particleParametersIn),
    currentTime(0.0),
    reductionFactor(reductionFactorIn)
{
    //Nothing else to do here
}

/// <summary>
/// Get the current error
/// </summary>
/// <returns></returns>
const double& OdeSolverParameters::getCurrentError() const
{
    return currentError;
}

/// <summary>
/// Set the current error
/// </summary>
/// <param name="currentErrorIn"></param>
void OdeSolverParameters::setCurrentError(const double& currentErrorIn)
{
    this->currentError = currentErrorIn;
}

/// <summary>
/// Get the total error
/// </summary>
/// <returns></returns>
const double& OdeSolverParameters::getTotalError() const
{
    return totalError;
}

/// <summary>
/// Set the total error
/// </summary>
/// <param name="totalErrorIn"></param>
void OdeSolverParameters::setTotalError(const double& totalErrorIn)
{
    this->totalError = totalErrorIn;
}

/// <summary>
/// Get min allowed error
/// </summary>
/// <returns></returns>
const double& OdeSolverParameters::getMinAllowedError() const
{
    return minAllowedError;
}

/// <summary>
/// Set min allowed error
/// </summary>
/// <param name="minAllowedErrorIn"></param>
void OdeSolverParameters::setMinAllowedError(const double& minAllowedErrorIn)
{
    this->minAllowedError = minAllowedErrorIn;
}

/// <summary>
/// Get max allowed error
/// </summary>
/// <returns></returns>
const double& OdeSolverParameters::getMaxAllowedError() const
{
    return maxAllowedError;
}

/// <summary>
/// Set max allowed error
/// </summary>
/// <param name="maxAllowedErrorIn"></param>
void OdeSolverParameters::setMaxAllowedError(const double& maxAllowedErrorIn)
{
    this->maxAllowedError = maxAllowedErrorIn;
}

/// <summary>
/// Get the current runtime
/// </summary>
/// <returns></returns>
const double& OdeSolverParameters::getCurrentRunTime() const
{
    return currentRunTime;
}

/// <summary>
/// Set the current runtime
/// </summary>
/// <param name="currentRunTimeIn"></param>
void OdeSolverParameters::setCurrentRunTime(const double& currentRunTimeIn)
{
    this->currentRunTime = currentRunTimeIn;
}

/// <summary>
/// Get the total runtime
/// </summary>
/// <returns></returns>
const double& OdeSolverParameters::getTotalRunTime() const
{
    return totalRunTime;
}

/// <summary>
/// Set the total runtime
/// </summary>
/// <param name="totalRunTimeIn"></param>
void OdeSolverParameters::setTotalRunTime(const double& totalRunTimeIn)
{
    this->totalRunTime = totalRunTimeIn;
}

/// <summary>
/// Get the richardson levels
/// </summary>
/// <returns></returns>
const size_t& OdeSolverParameters::getRichardsonLevels() const
{
    return richardsonLevels;
}

/// <summary>
/// Set the richardson levels
/// </summary>
/// <param name="richardsonLevelsIn"></param>
void OdeSolverParameters::setRichardsonLevels(const size_t& richardsonLevelsIn)
{
    this->richardsonLevels = richardsonLevelsIn;
}

/// <summary>
/// Get the mininum richardson tables
/// </summary>
/// <returns></returns>
const size_t& OdeSolverParameters::getMinRichardsonTables() const
{
    return minRichardsonTables;
}

/// <summary>
/// Get the mininum richardson tables
/// </summary>
/// <param name="minRichardsonTablesIn"></param>
void OdeSolverParameters::setMinRichardsonTables(const size_t& minRichardsonTablesIn)
{
    this->minRichardsonTables = minRichardsonTablesIn;
}

/// <summary>
/// Get the max number of richardson tables
/// </summary>
/// <returns></returns>
const size_t& OdeSolverParameters::getMaxRichardsonTables() const
{
    return maxRichardsonTables;
}

/// <summary>
/// Set the max number of richardson tables
/// </summary>
/// <param name="maxRichardsonTablesIn"></param>
void OdeSolverParameters::setMaxRichardsonTables(const size_t& maxRichardsonTablesIn)
{
    this->maxRichardsonTables = maxRichardsonTablesIn;
}

/// <summary>
/// Get the convergence 
/// </summary>
/// <returns></returns>
const double& OdeSolverParameters::getC() const
{
    return c;
}

/// <summary>
/// Set the convergence
/// </summary>
/// <param name="cIn"></param>
void OdeSolverParameters::setC(const double& cIn)
{
    this->c = cIn;
}

/// <summary>
/// Get delta time
/// </summary>
/// <returns></returns>
const double& OdeSolverParameters::getDt() const
{
    return dt;
}

/// <summary>
/// Set delta time
/// </summary>
/// <param name="dtIn"></param>
void OdeSolverParameters::setDt(const double& dtIn)
{
    this->dt = dtIn;
}

/// <summary>
/// Get min dt allowed
/// </summary>
/// <returns></returns>
const double& OdeSolverParameters::getMinDt() const
{
    return minDt;
}

/// <summary>
/// Set min dt allowed
/// </summary>
/// <param name="minDtIn"></param>
void OdeSolverParameters::setMinDt(const double& minDtIn)
{
    this->minDt = minDtIn;
}

/// <summary>
/// Get the max allowed dt
/// </summary>
/// <returns></returns>
const double& OdeSolverParameters::getMaxDt() const
{
    return maxDt;
}

/// <summary>
/// Set the max allowed dt
/// </summary>
/// <param name="maxDtIn"></param>
void OdeSolverParameters::setMaxDt(const double& maxDtIn)
{
    this->maxDt = maxDtIn;
}

/// <summary>
/// Verifiy if our parameters are valid
/// </summary>
/// <returns></returns>
const bool OdeSolverParameters::paramsValid() const
{
    //Initalize our set to true
    bool tempGoodness = true;

    //Check if error bounds are valid
    tempGoodness &= isUniqueAndOrdered(minAllowedError, maxAllowedError);

    //Check richardson levels
    tempGoodness &= (minRichardsonTables >= 2) && (maxRichardsonTables > minRichardsonTables);

    //Particle parameters are valid
    tempGoodness &= (particleParams.getNumOfParticles() > 0) && ((particleParams.getLearningRate() >= 0.0) && (particleParams.getLearningRate() <= 1.0) && particleParams.getMaxPSOIterations() > 0 && particleParams.getPSOTol() > 0);

    //Check dt parameters and return
    return tempGoodness &= isUniqueAndOrdered(minDt, maxDt);
}

/// <summary>
/// Get if euler is allowed
/// </summary>
/// <returns></returns>
const bool& OdeSolverParameters::getAllowedEuler() const
{
    return allowedEuler;
}

/// <summary>
/// Set is euler allowed
/// </summary>
/// <param name="allowedEulerIn"></param>
void OdeSolverParameters::setAllowedEuler(const bool& allowedEulerIn)
{
    allowedEuler = allowedEulerIn;
}

/// <summary>
/// Get if rk2 is allowed
/// </summary>
/// <returns></returns>
const bool& OdeSolverParameters::getAllowedRK2() const
{
    return allowedRK2;
}

/// <summary>
/// Set rk2 is allowed
/// </summary>
/// <param name="allowedRK2In"></param>
void OdeSolverParameters::setAllowedRK2(const bool& allowedRK2In)
{
    allowedRK2 = allowedRK2In;
}

/// <summary>
/// Get if rk4 is allowed
/// </summary>
/// <returns></returns>
const bool& OdeSolverParameters::getAllowedRK4() const
{
    return allowedRK4;
}

/// <summary>
/// Set if rk4 is allowed
/// </summary>
/// <param name="allowedRK4In"></param>
void OdeSolverParameters::setAllowedRK4(const bool& allowedRK4In)
{
    allowedRK4 = allowedRK4In;
}

/// <summary>
/// Get if system is stiff
/// </summary>
/// <returns></returns>
const bool& OdeSolverParameters::getIsSystemStiff() const
{
    return isSystemStiff;
}

/// <summary>
/// Set if system is stiff
/// </summary>
/// <param name="isSystemStiffIn"></param>
void OdeSolverParameters::setIsSystemStiff(const bool& isSystemStiffIn)
{
    isSystemStiff = isSystemStiffIn;
}

/// <summary>
/// Get our particle parameters const
/// </summary>
/// <returns></returns>
const ParticleParameters& odepso::OdeSolverParameters::getParticleParameters() const
{
    return particleParams;
}

/// <summary>
/// Get the handle to the particle params for modification
/// </summary>
/// <returns></returns>
ParticleParameters& odepso::OdeSolverParameters::getParticleParameters()
{
    return particleParams;
}

/// <summary>
/// Set our particle parameters
/// </summary>
/// <param name="particleParametersIn"></param>
void odepso::OdeSolverParameters::setParticleParameters(const ParticleParameters& particleParametersIn)
{
    particleParams = particleParametersIn;
}

/// <summary>
/// Get the current time
/// </summary>
/// <returns></returns>
const double& odepso::OdeSolverParameters::getCurrentTime() const
{
    return currentTime;
}

/// <summary>
/// Set the current time
/// </summary>
/// <param name="currentTimeIn"></param>
void odepso::OdeSolverParameters::setCurrentTime(const double& currentTimeIn)
{
    currentTime = currentTimeIn;
}

/// <summary>
/// Get the reduction factor
/// </summary>
/// <returns></returns>
const double& odepso::OdeSolverParameters::getReductionFactor() const
{
    return reductionFactor;
}

/// <summary>
/// Set the reduction factor
/// </summary>
/// <param name="reductionFactorIn"></param>
void odepso::OdeSolverParameters::setReductionFactor(const double& reductionFactorIn)
{
    reductionFactor = reductionFactorIn;
}



