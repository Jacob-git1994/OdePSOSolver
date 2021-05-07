#include "ParticleProcessor.h"

void odepso::ParticleProcessor::runParticle(std::unique_ptr<SolverIF>& solverTypeIn, odepso::Particle& particleIn, const Eigen::VectorXd& currentStateIn)
{
	//Get a temp state so the next pass through knows the updated result
	Eigen::VectorXd currentStateTemp = currentStateIn;

	//Propograte through time
	while (particleIn.getParameters().getCurrentTime() < wayPointEnd)
	{
		//Propagrate solution forward in time
		particleIn.updateTimeStep(solverTypeIn, currentStateTemp, wayPointEnd);

		//Update the state vector
		currentStateTemp = particleIn.getBestResult();
	}
}

bool odepso::ParticleProcessor::compParticle(const Particle& lhs, const Particle& rhs) const
{
	//The errors are comparable
	if (globalErrorInRange(lhs) && globalErrorInRange(rhs))
	{
		//Best error return lhs.estimateGlobalError(wayPointEnd) < rhs.estimateGlobalError(wayPointEnd);
		//Get fastest running in range
		//return lhs.getParameters().getDt() > rhs.getParameters().getDt();
		return lhs.getParameters().getTotalRunTime() < rhs.getParameters().getTotalRunTime();
	}
	else if (globalErrorInRange(lhs) && !globalErrorInRange(rhs))
	{
		return true;
	}
	else if (!globalErrorInRange(lhs) && globalErrorInRange(rhs))
	{
		return false;
	}
	else //Pick the particle closest to the set
	{
		if (lhs.estimateGlobalError(finalEndTime) > bestParameters.getMaxAllowedError() && rhs.estimateGlobalError(finalEndTime) > bestParameters.getMaxAllowedError())
		{
			if (std::fabs(lhs.estimateGlobalError(finalEndTime) - bestParameters.getMaxAllowedError()) < std::fabs(rhs.estimateGlobalError(finalEndTime) - bestParameters.getMaxAllowedError()))
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		else if (lhs.estimateGlobalError(finalEndTime) > bestParameters.getMaxAllowedError() && rhs.estimateGlobalError(finalEndTime) < bestParameters.getMinAllowedError())
		{
			if (std::fabs(lhs.estimateGlobalError(finalEndTime) - bestParameters.getMaxAllowedError()) < std::fabs(rhs.estimateGlobalError(finalEndTime) - bestParameters.getMinAllowedError()))
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		else if (lhs.estimateGlobalError(finalEndTime) < bestParameters.getMinAllowedError() && rhs.estimateGlobalError(finalEndTime) > bestParameters.getMaxAllowedError())
		{
			if (std::fabs(lhs.estimateGlobalError(finalEndTime) - bestParameters.getMinAllowedError()) < std::fabs(rhs.estimateGlobalError(finalEndTime) - bestParameters.getMaxAllowedError()))
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			if (std::fabs(lhs.estimateGlobalError(finalEndTime) - bestParameters.getMinAllowedError()) < std::fabs(rhs.estimateGlobalError(finalEndTime) - bestParameters.getMinAllowedError()))
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		//Default case
		return false;
	}
}

bool odepso::ParticleProcessor::globalErrorInRange(const Particle& particleIn) const
{
	return particleIn.estimateGlobalError(finalEndTime) <= bestParameters.getMaxAllowedError() && particleIn.estimateGlobalError(finalEndTime) >= bestParameters.getMinAllowedError();
	//return particleIn.getParameters().getTotalError() <= bestParameters.getMaxAllowedError() && particleIn.getParameters().getTotalError() >= bestParameters.getMinAllowedError();
}

bool odepso::ParticleProcessor::isLessThanMaxError(const Particle& particleIn) const
{
	return (particleIn.estimateGlobalError(finalEndTime) <= bestParameters.getMaxAllowedError()) && (particleIn.estimateGlobalError(finalEndTime) > std::numeric_limits<double>::epsilon());
}

odepso::Particle odepso::ParticleProcessor::generateAndSetParticle(std::unique_ptr<SolverIF>& solverTypeIn, const Particle& particleIn, Particle& bestLocalParticleIn, std::minstd_rand& randIn, const Eigen::VectorXd& currentStateIn)
{
	//Gather important parameters from this particle
	const double& dt = particleIn.getParameters().getDt();
	const double richLevels = static_cast<double>(particleIn.getParameters().getRichardsonLevels());
	const double& dtVelocity = particleIn.getDtVelocity();
	const double& richVelocity = particleIn.getRichardsonVelocity();
	const double& omega = particleIn.getParameters().getParticleParameters().getOmega();
	const double& phiOne = particleIn.getParameters().getParticleParameters().getPhi1();
	const double& phiTwo = particleIn.getParameters().getParticleParameters().getPhi2();
	const double& learningRate = particleIn.getParameters().getParticleParameters().getLearningRate();

	//Gather important information about this particle's history
	const double& bestKnownDt = bestLocalParticleIn.getParameters().getDt();
	const double bestKnownRich = static_cast<double>(bestLocalParticleIn.getParameters().getRichardsonLevels());

	//Gather Information about the global best particle
	const double& globalBestDt = bestParticle.getParameters().getDt();
	const double globalBestRich = static_cast<double>(bestParticle.getParameters().getRichardsonLevels());

	//Uniform random generators
	std::uniform_real_distribution<double> uniformRandomGen;

	//Generate our random numbers
	const double rpDt = uniformRandomGen(randIn);
	const double rgDt = uniformRandomGen(randIn);
	const double rpRich = uniformRandomGen(randIn);
	const double rgRich = uniformRandomGen(randIn);

	//Update the particles velocity
	const double newDtVelocity = omega * dtVelocity + phiOne * rpDt * (bestKnownDt - dt) + phiTwo * rgDt * (globalBestDt - dt);
	const double newRichVelocity = omega * richVelocity + phiOne * rpRich * (bestKnownRich - richLevels) + phiTwo * rgRich * (globalBestRich - richLevels);

	//Update the new position of the particle
	const double newDt = dt + learningRate * newDtVelocity;
	double integerComp = 0;
	size_t newRich = 0;
	bool ceilRich = (std::modf(richLevels + learningRate * newRichVelocity, &integerComp) > uniformRandomGen(randIn));
	if (ceilRich)
	{
		newRich = static_cast<size_t>(integerComp + 1);
	}
	else
	{
		newRich = static_cast<size_t>(integerComp);
	}

	//Create a new set of Parameters
	odepso::OdeSolverParameters newParams = bestParameters;

	//Set the new variables for our parameters
	newParams.setDt(newDt);
	newParams.setRichardsonLevels(newRich);

	//Ensure the parameters are valid
	ensureAllParametersAreValid(newParams);

	//Temp Particle based on current set of parameters
	odepso::Particle tempParticle(newParams, odepso::Richardson());

	//Update the velocity
	tempParticle.setDtVelocity(newDtVelocity);
	tempParticle.setRichardsonVelocity(newRichVelocity);

	//Run the temp particle
	runParticle(solverTypeIn, tempParticle, currentStateIn);

	//Check if this particle is better than its history
	if (compParticle(tempParticle, bestLocalParticleIn))
	{
		bestLocalParticleIn = tempParticle;

		//Update global particle is this is way better
		if (compParticle(bestLocalParticleIn, bestParticle))
		{
			bestParticle = bestLocalParticleIn;
		}
	}

	return tempParticle;
}

void odepso::ParticleProcessor::ensureAllParametersAreValid(OdeSolverParameters& paramsIn)
{
	//TODO Check if all parameters are valid
	//std::cout << paramsIn.getDt() << std::endl;
	//std::cout << paramsIn.getRichardsonLevels() << std::endl;
	//Dt out of bounds
	if ((paramsIn.getDt() < paramsIn.getMinDt()))
	{
		paramsIn.setDt(paramsIn.getMinDt());
	}
	else if (paramsIn.getDt() > paramsIn.getMaxDt())
	{
		paramsIn.setDt(paramsIn.getMaxDt());
	}

	//Richard levels out of bounds
	if(paramsIn.getRichardsonLevels() > paramsIn.getMaxRichardsonTables())
	{
		paramsIn.setRichardsonLevels(paramsIn.getMaxRichardsonTables());
	}
	else if (paramsIn.getRichardsonLevels() < paramsIn.getMinRichardsonTables())
	{
		paramsIn.setRichardsonLevels(paramsIn.getMinRichardsonTables());
	}

	//std::cout << std::endl << paramsIn.getDt() << std::endl;
	//std::cout << paramsIn.getRichardsonLevels() << std::endl;
}

void odepso::ParticleProcessor::calcStatistics()
{
	//Clear previous results
	expectedDt = 0.0;
	expectedRichardson = 0.0;
	varDt = 0.0;
	varRich = 0.0;

	//Get expected dt
	std::for_each(particles.begin(), particles.end(), [&](const odepso::Particle& particleIn)
		{
			expectedDt += particleIn.getParameters().getDt();
		});

	//Update expected dt
	expectedDt *= 1. / static_cast<double>(bestParameters.getParticleParameters().getNumOfParticles());

	//Get expected richardson
	std::for_each(particles.begin(), particles.end(), [&](const odepso::Particle& particleIn)
		{
			expectedRichardson += static_cast<double>(particleIn.getParameters().getRichardsonLevels());
		});

	//Update expected richardson
	expectedRichardson *= 1. / static_cast<double>(bestParameters.getParticleParameters().getNumOfParticles());

	//Calculate dt variance
	std::for_each(particles.begin(), particles.end(), [&](const odepso::Particle& particleIn)
		{
			varDt += (particleIn.getParameters().getDt() - expectedDt) * (particleIn.getParameters().getDt() - expectedDt);
		});

	//Calculate richardson variance
	std::for_each(particles.begin(), particles.end(), [&](const odepso::Particle& particleIn)
		{
			varRich += (static_cast<double>(particleIn.getParameters().getRichardsonLevels()) - expectedRichardson) * (static_cast<double>(particleIn.getParameters().getRichardsonLevels()) - expectedRichardson);
		});

	//Update the variance calculations
	varDt *= 1. / static_cast<double>(bestParameters.getParticleParameters().getNumOfParticles() - 1);
	varRich *= 1. / static_cast<double>(bestParameters.getParticleParameters().getNumOfParticles() - 1);
}

odepso::ParticleProcessor::ParticleProcessor() :
	bestParameters(OdeSolverParameters()),
	particles({}),
	bestLocalParticles({}),
	wayPointBegin(0.0),
	wayPointEnd(0.0),
	bestParticle(odepso::Particle()),
	expectedDt(0.0),
	expectedRichardson(0.0),
	varDt(std::numeric_limits<double>::max()),
	varRich(std::numeric_limits<double>::max()),
	finalEndTime(0.0)
{
	//Nothing to do here
}

void odepso::ParticleProcessor::setParamsAndWayPoint(const OdeSolverParameters& odeSolverParametersIn, const double& wayPointBeginIn, const double& wayPointEndIn, const double& finalEndTimeIn)
{
	if (wayPointEnd < wayPointBegin)
	{
		throw std::logic_error("Way points inproperly ordered");
	}

	//Copy infomation
	bestParameters = odeSolverParametersIn;
	wayPointBegin = wayPointBeginIn;
	wayPointEnd = wayPointEndIn;
	finalEndTime = finalEndTimeIn;

	//Clear all our lists
	particles.clear();
	bestLocalParticles.clear();

	//Clear best particle
	bestParticle = odepso::Particle();
}

const odepso::OdeSolverParameters& odepso::ParticleProcessor::getParams() const
{
	return bestParameters;
}

const odepso::Particle& odepso::ParticleProcessor::getBestParticle() const
{
	return bestParticle;
}

void odepso::ParticleProcessor::PSO(std::unique_ptr<SolverIF>& solverTypeIn, std::minstd_rand& randIn, const Eigen::VectorXd& currentStateIn)
{
	//Get number of particles
	const size_t& numParticles = bestParameters.getParticleParameters().getNumOfParticles();
	const size_t& MAX_ITER = bestParameters.getParticleParameters().getMaxPSOIterations();
	const double& CONVERGENCE_TOL = bestParameters.getParticleParameters().getPSOTol();
	size_t currentIterationCount = 0;

	//Build our particles
	for (size_t i = 0; i < numParticles; ++i)
	{
		//Create our partile
		odepso::Particle tempParticle(bestParameters, odepso::Richardson());

		//Initialize the state and velocity
		tempParticle.initalize(randIn);

		//Iterate this particle through the way points
		runParticle(solverTypeIn, tempParticle, currentStateIn);

		//Add particle to particle vector
		particles.push_back(tempParticle);
		bestLocalParticles.push_back(tempParticle);

		//Check if this value is better than the global particle. Also adds if this is the first pass through
		if (compParticle(tempParticle, bestParticle) || (i == 0))
		{
			bestParticle = tempParticle;
		}
	}

	//PSO
	while ((currentIterationCount++ < MAX_ITER) && (((varDt > CONVERGENCE_TOL) || varRich > CONVERGENCE_TOL) || (!globalErrorInRange(bestParticle))))
	{
		//Get iterators to be best particle and current list of particles
		std::vector<odepso::Particle>::iterator partItr = particles.begin();
		std::vector<odepso::Particle>::iterator localBestItr = bestLocalParticles.begin();
		for (; ((partItr != particles.end()) && (localBestItr != bestLocalParticles.end())); )
		{
			//Get particles
			odepso::Particle& currentParticle = *partItr;
			odepso::Particle& thisBestParticle = *localBestItr;

			//Clear this particles information
			currentParticle.clearParticleHistory(bestParameters);

			//Create our new particle and compare it. It will update this "Best" of this particle if it is better. It will update the global particle if it is better.
			currentParticle = generateAndSetParticle(solverTypeIn, currentParticle, thisBestParticle, randIn, currentStateIn);

			//Update to the next particle
			++partItr;
			++localBestItr;

			//std::cout << (partItr - 1)->getParameters().getRichardsonLevels() << "\n";// << (partItr - 1)->getParameters().getDt() << "\n";
			//std::cout << (partItr - 1)->getParameters().getTotalRunTime() << "\t" << (partItr - 1)->estimateGlobalError(finalEndTime) << "\n";
		}

		//Calcualte statistics from the previous pass through
		calcStatistics();

		std::cout << bestParticle.getParameters().getDt() << "\t" << bestParticle.getParameters().getRichardsonLevels() << "\t" << bestParticle.getParameters().getTotalRunTime() << "\n";
		std::cout << varDt << "\t" << varRich << "\n";
	}

	//Update the best parameters for this best particle
	bestParameters = bestParticle.getParameters();
}

const double& odepso::ParticleProcessor::getExpectedDt() const
{
	return expectedDt;
}

const double& odepso::ParticleProcessor::getExpectedRich() const
{
	return expectedRichardson;
}

const double& odepso::ParticleProcessor::getVarDt() const
{
	return varDt;
}

const double& odepso::ParticleProcessor::getVarRich() const
{
	return varRich;
}



