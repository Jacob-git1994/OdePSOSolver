#include "MethodWrapper.h"

odepso::MethodWrapper::MethodWrapper() :
	methodsMap(std::map<Common::SOLVER_TYPES, std::unique_ptr<SolverIF>, Common::CompSolvers>()),
	paramsMap(std::map<Common::SOLVER_TYPES, OdeSolverParameters, Common::CompSolvers>()),
	resultsMap(std::map<Common::SOLVER_TYPES, std::vector<std::pair<double, Eigen::VectorXd>>, Common::CompSolvers>()),
	psoProcessorMap(std::map<Common::SOLVER_TYPES, ParticleProcessor, Common::CompSolvers>()),
	randomStruct(RandomStruct())
{
	//Nothing else to do here
}

bool odepso::MethodWrapper::initialize(const OdeSolverParameters& paramsIn)
{
	//Make sure everything is valid
	if (!paramsIn.paramsValid())
	{
		return false;
	}

	//If valid clear our lists
	methodsMap.clear();
	paramsMap.clear();
	psoProcessorMap.clear();
	resultsMap.clear();
	randomStruct = RandomStruct();

	try
	{
		//Build Euler
		if (paramsIn.getAllowedEuler() && !paramsIn.getIsSystemStiff())
		{
			methodsMap.emplace(Common::SOLVER_TYPES::EULER, std::move(std::unique_ptr<odepso::SolverIF>(new Euler())));
		}

		//Build RK2
		if (paramsIn.getAllowedRK2() && !paramsIn.getIsSystemStiff())
		{
			methodsMap.emplace(Common::SOLVER_TYPES::RK2, std::move(std::unique_ptr<odepso::SolverIF>(new RungeKutta2())));
		}

		//Build RK4
		if (paramsIn.getAllowedRK4() && !paramsIn.getIsSystemStiff())
		{
			methodsMap.emplace(Common::SOLVER_TYPES::RK4, std::move(std::unique_ptr<odepso::SolverIF>(new RungeKutta4())));
		}

		//Build our remaining maps
		for (std::map<Common::SOLVER_TYPES, std::unique_ptr<SolverIF>, Common::CompSolvers>::const_iterator methodItr = methodsMap.cbegin();
			methodItr != methodsMap.cend(); ++methodItr)
		{
			//Build parameter map
			paramsMap.emplace(methodItr->first, paramsIn);

			//Build result map
			resultsMap.emplace(methodItr->first, std::move(std::vector<std::pair<double, Eigen::VectorXd>>()));

			//Build our particle processor map
			psoProcessorMap.emplace(methodItr->first, std::move(odepso::ParticleProcessor()));
		}

		//Set up our random generators
		randomStruct.seed = paramsIn.getParticleParameters().getSeed();
		randomStruct.randomEngine = std::minstd_rand(static_cast<unsigned int>(randomStruct.seed));
	}
	catch (std::exception& e)
	{
		std::cerr << e.what();
		return false;
	}

	//If everything passed return true
	return true;
}
