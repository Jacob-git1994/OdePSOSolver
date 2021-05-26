#include "Particle.h"

odepso::Particle::Particle() :
	particleParameters(odepso::OdeSolverParameters()),
	particleTables(odepso::Richardson()),
	errorCost(0.0),
	timeCost(0.0),
	errorVector({}),
	dtVelocity(0.0),
	richardsonVelocity(0.0),
	bestResult(Eigen::VectorXd::Zero(0))
{
	//Nothing else to do here
}

odepso::Particle::Particle(
	const odepso::OdeSolverParameters& odeSolverParametersIn,
	const odepso::Richardson& richardsonIn) :
	particleParameters(odeSolverParametersIn),
	particleTables(richardsonIn),
	errorCost(0.0),
	timeCost(0.0),
	errorVector({}),
	dtVelocity(0.0),
	richardsonVelocity(0.0),
	bestResult(Eigen::VectorXd::Zero(0))
{
	//Nothing else to do here
}

odepso::Particle::Particle(const odepso::Particle& particleIn) :
	particleParameters(particleIn.particleParameters),
	particleTables(particleIn.particleTables),
	errorCost(particleIn.errorCost),
	timeCost(particleIn.timeCost),
	errorVector(particleIn.errorVector),
	dtVelocity(particleIn.dtVelocity),
	richardsonVelocity(particleIn.richardsonVelocity),
	bestResult(particleIn.bestResult)
{
	//Nothing else to do here
}

odepso::Particle& odepso::Particle::operator=(const odepso::Particle& particleIn)
{
	particleParameters = particleIn.particleParameters;
	particleTables = particleIn.particleTables;
	errorCost = particleIn.errorCost;
	timeCost = particleIn.timeCost;
	errorVector = particleIn.errorVector;
	dtVelocity = particleIn.dtVelocity;
	richardsonVelocity = particleIn.richardsonVelocity;
	bestResult = particleIn.bestResult;

	return *this;
}

void odepso::Particle::initalize(std::minstd_rand& randEngineIn)
{
	//Random Generators
	std::uniform_real_distribution<double> dtRandomGen(particleParameters.getMinDt(), particleParameters.getMaxDt());
	std::uniform_int_distribution<size_t> tableRandomGen(particleParameters.getMinRichardsonTables(), particleParameters.getMaxRichardsonTables());
	std::uniform_real_distribution<double> dtV(-std::fabs(particleParameters.getMaxDt() - particleParameters.getMinDt()), std::fabs(particleParameters.getMaxDt() - particleParameters.getMinDt()));
	std::uniform_real_distribution<double> dtR(-std::fabs(static_cast<double>(particleParameters.getMaxRichardsonTables()) - static_cast<double>(particleParameters.getMinRichardsonTables())), std::fabs(static_cast<double>(particleParameters.getMaxRichardsonTables()) - static_cast<double>(particleParameters.getMinRichardsonTables())));

	//Initialize our parameteres
	particleParameters.setDt(dtRandomGen(randEngineIn));
	particleParameters.setRichardsonLevels(tableRandomGen(randEngineIn));
	dtVelocity = dtV(randEngineIn);
	richardsonVelocity = dtR(randEngineIn);

	//Make sure everything is valid
	if (!particleParameters.paramsValid())
	{
		throw std::invalid_argument("Bad Solver Parameters");
	}
}

void odepso::Particle::updateTimeStep(std::unique_ptr<SolverIF>& solverTypeIn, OdeSolverParameters& odeSolverParametersIn, const Eigen::VectorXd& currentStateIn, const double& finalTimeIn)
{
	//Gather parameters
	const double& currentTime = particleParameters.getCurrentTime();
	const double currentDt = (particleParameters.getDt() + currentTime > finalTimeIn ? finalTimeIn - currentTime : particleParameters.getDt()); //If we need to clamp
	const size_t& currentTableSize = particleParameters.getRichardsonLevels();

	//Initalize our dt and richardson tables
	double tempDt = 0.0;
	unsigned int numberOfSteps = 0;

	//Run Time
	double runTime = 0.0;

	//Initalize our clock
	std::chrono::time_point<std::chrono::high_resolution_clock> tNow  = std::chrono::time_point<std::chrono::high_resolution_clock>::time_point::min();
	std::chrono::time_point<std::chrono::high_resolution_clock> tEnd = std::chrono::time_point<std::chrono::high_resolution_clock>::time_point::min();

	//Temp vector to hold result
	Eigen::VectorXd tempState = currentStateIn;

	//Update the richardson tables
	particleTables.reset(
		currentTableSize,
		currentStateIn.size(), 
		particleParameters.getReductionFactor(), 
		currentDt, 
		solverTypeIn->getTheoreticalConvergence());

	//Record current time
	tNow = std::chrono::high_resolution_clock::now();

	//Go through all types of table sizes
	for (size_t i = 0; i < particleTables.getTableSize(); ++i)
	{
		//Update new time step and number of steps to statify that time step
		tempDt = currentDt / std::pow(particleParameters.getReductionFactor(), static_cast<double>(i));
		numberOfSteps = static_cast<unsigned int>(std::pow(particleParameters.getReductionFactor(), static_cast<double>(i)));

		//Propagate forward in time
		solverTypeIn->updateTimeStep(currentStateIn, tempState, tempDt, numberOfSteps, currentTime);

		//Append to our tables
		particleTables.append(i, tempState);
	}

	//Extrapolate our results
	particleTables.calculateRichardsonExtrapolation();

	//End Timing 
	tEnd = std::chrono::high_resolution_clock::now();

	//Calculate the run time in seconds
	runTime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(tEnd - tNow).count());

	//Update all the parameters
	particleParameters.setC(particleTables.getC());
	particleParameters.setCurrentError(particleTables.getError());
	particleParameters.setTotalError(particleTables.getError() + particleParameters.getTotalError());
	particleParameters.setCurrentRunTime(runTime);
	particleParameters.setTotalRunTime(runTime + particleParameters.getTotalRunTime());
	particleParameters.setDt(particleParameters.getDt()); //Set the the known dt not the fixed one
	errorCost = particleTables.getError();
	timeCost = runTime;
	bestResult = particleTables.getBestResult();
	errorVector.push_back(particleParameters.getTotalError());

	//Update the results to the next time
	particleParameters.setCurrentTime(currentDt + particleParameters.getCurrentTime());

	//Set the parameters to the new set of parameters found here
	odeSolverParametersIn = particleParameters;
}

void odepso::Particle::updateTimeStep(std::unique_ptr<SolverIF>& solverTypeIn, const Eigen::VectorXd& currentStateIn, const double& finalTimeIn)
{
	Particle::updateTimeStep(solverTypeIn, particleParameters, currentStateIn, finalTimeIn);
}

const Eigen::VectorXd& odepso::Particle::getBestResult() const
{
	return bestResult;
}

const odepso::OdeSolverParameters& odepso::Particle::getParameters() const
{
	return particleParameters;
}

const double& odepso::Particle::getErrorCost() const
{
	return errorCost;
}

const double& odepso::Particle::getTimeCost() const
{
	return timeCost;
}

const double odepso::Particle::estimateGlobalError(const double& finalTimeIn) const
{
	const size_t& vecSize = errorVector.size();

	double estimatedGlobalError = 0.;

	if (vecSize == 0)
	{
		estimatedGlobalError = particleParameters.getCurrentError();
	}
	else if (vecSize == 1)
	{
		estimatedGlobalError = particleParameters.getTotalError() + ((finalTimeIn - particleParameters.getCurrentTime()) / particleParameters.getDt()) * errorVector.at(0);
	}
	else if (vecSize == 2)
	{
		estimatedGlobalError = particleParameters.getTotalError() + (finalTimeIn - particleParameters.getCurrentTime()) * (errorVector.at(vecSize - 1) - errorVector.at(vecSize - 2)) / particleParameters.getDt();
	}
	else
	{
		estimatedGlobalError = particleParameters.getTotalError() +
			std::fabs((finalTimeIn - particleParameters.getCurrentTime()) * (errorVector.at(vecSize - 1) - errorVector.at(vecSize - 3)) / 2. * particleParameters.getDt()) +
			std::fabs(.5 * std::pow((finalTimeIn - particleParameters.getCurrentTime()), 2.0) * (errorVector.at(vecSize - 1) - 2.0 * errorVector.at(vecSize - 2) + errorVector.at(vecSize - 3)) / (particleParameters.getDt() * particleParameters.getDt()));
	}
	return estimatedGlobalError;
}

void odepso::Particle::setParticleParameters(const OdeSolverParameters& odeSolverParametersIn)
{
	particleParameters = odeSolverParametersIn;
}

const double& odepso::Particle::getDtVelocity() const
{
	return dtVelocity;
}

void odepso::Particle::setDtVelocity(const double& dtVelocityIn)
{
	dtVelocity = dtVelocityIn;
}

const double& odepso::Particle::getRichardsonVelocity() const
{
	return richardsonVelocity;
}

void odepso::Particle::setRichardsonVelocity(const double& richardsonVelocityIn)
{
	richardsonVelocity = richardsonVelocityIn;
}

void odepso::Particle::clearParticleHistory(const OdeSolverParameters& paramsIn)
{
	particleParameters.setC(0.0);
	particleParameters.setCurrentError(0.0);
	particleParameters.setCurrentRunTime(0.0);
	particleParameters.setCurrentTime(paramsIn.getCurrentTime());
	particleParameters.setTotalError(paramsIn.getTotalError());
	particleParameters.setTotalRunTime(paramsIn.getTotalRunTime());
}

void odepso::Particle::setBestResult(const Eigen::VectorXd& bestStateIn)
{
	bestResult = bestStateIn;
}
