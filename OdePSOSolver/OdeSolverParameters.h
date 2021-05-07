#pragma once

#include "ParticleParameters.h"

#include <cmath>
#include <limits>

namespace odepso
{
	class OdeSolverParameters
	{
	private:

		//Current truncation error
		double currentError;

		//Total accumulated error
		double totalError;

		//Mininum allowed error
		double minAllowedError;

		//Max allowed error
		double maxAllowedError;

		//Current runtime for the step
		double currentRunTime;

		//Total runtime for the step
		double totalRunTime;

		//Number of richardson levels
		size_t richardsonLevels;

		//Mininum allowed tables
		size_t minRichardsonTables;

		//Maximum allowed tables
		size_t maxRichardsonTables;

		//Reduction factor
		double reductionFactor;

		//Our convergence
		double c;

		//Current time step
		double dt;

		//Mininum delta time
		double minDt;

		//Maximim delta time
		double maxDt;

		//Ensure double ranges are unique and valid
		const bool isUniqueAndOrdered(const double& lowerIn, const double& upperIn) const;

		//Use this Euler's method
		bool allowedEuler;

		//Use Runge Kutta two method
		bool allowedRK2;

		//Use Runge Kutta four method
		bool allowedRK4;

		//Is this system stiff
		bool isSystemStiff;

		//The particle parameters
		ParticleParameters particleParams;

		//The current Time
		double currentTime;

	public:

		//Default operator
		OdeSolverParameters();

		//Constructor to set our parameters
		OdeSolverParameters(
			const double& minAllowedErrorIn,
			const double& maxAllowedErrorIn,
			const size_t& minRichardsonTablesIn,
			const size_t& maxRichardsonTablesIn,
			const double& minDtIn,
			const double& maxDtIn,
			const ParticleParameters& particleParametersIn,
			const double& reductionFactor);

		//Default Copy Constructor
		OdeSolverParameters(const OdeSolverParameters& OdeSolverParametersIn) = default;

		//Default Assign Constructor
		OdeSolverParameters& operator=(const OdeSolverParameters & OdeSolverParametersIn) = default;

		//Default Destructor
		virtual ~OdeSolverParameters() = default;

		//Get the current error
		const double& getCurrentError() const;

		//Set the current error
        void setCurrentError(const double& currentErrorIn);

		//Get the total error
        const double& getTotalError() const;

		//Set the total error
        void setTotalError(const double& totalErrorIn);

		//Get the min allowed error
		const double& getMinAllowedError() const;

		//Set the min allowed error
        void setMinAllowedError(const double& minAllowedErrorIn);

		//Get the max allowed error
		const double& getMaxAllowedError() const;

		//Set the max allowed error
        void setMaxAllowedError(const double& maxAllowedErrorIn);

		//Get the current run time
		const double& getCurrentRunTime() const;

		//Set the current runtime
        void setCurrentRunTime(const double& currentRunTimeIn);

		//Get the total run time
		const double& getTotalRunTime() const;

		//Set the total run time
        void setTotalRunTime(const double& totalRunTimeIn);

		//Get the number of current richardson levels
        const size_t& getRichardsonLevels() const;

		//Set the number of richardson levels
        void setRichardsonLevels(const size_t& richardsonLevelsIn);

		//Get the min number of richardson levels
        const size_t& getMinRichardsonTables() const;

		//Set the min number of richardson levels
        void setMinRichardsonTables(const size_t& minRichardsonTablesIn);

		//Get the max number of richardson levels
        const size_t& getMaxRichardsonTables() const;

		//Set the max number of richardson levels
        void setMaxRichardsonTables(const size_t& maxRichardsonTablesIn);

		//Get the convergence of the tables
		const double& getC() const;

		//Set the convergence of the tables
        void setC(const double& cIn);

		//Get the current time step
		const double& getDt() const;

		//Set the current time step
        void setDt(const double& dtIn);

		//Get the mininum dt allowed
		const double& getMinDt() const;

		//Set the mininum dt allowed
        void setMinDt(const double& minDtIn);

		//Get the maximum dt allowed
		const double& getMaxDt() const;

		//Set the maxinum dt allowed
        void setMaxDt(const double& maxDtIn);

		//Ensure the parameters are valid
		const bool paramsValid() const;

		//Get if we can use Euler
		const bool& getAllowedEuler() const;

		//Set Euler
		void setAllowedEuler(const bool& allowedEulerIn);

		//Get if we can use rk2
		const bool& getAllowedRK2() const;

		//Set if we can use rk2
		void setAllowedRK2(const bool& allowedRK2In);

		//Get if we can use rk4
		const bool& getAllowedRK4() const;

		//Set if we can use rk4
		void setAllowedRK4(const bool& allowedRK4In);

		//Get if our system is stiff
		const bool& getIsSystemStiff() const;

		//Set if our system is stiff
		void setIsSystemStiff(const bool& isSystemStiffIn);

		//Get our particle parameters
		const ParticleParameters& getParticleParameters() const;

		//Get our particle params handle non-const
		ParticleParameters& getParticleParameters();

		//Set our partical parameters
		void setParticleParameters(const ParticleParameters& particleParametersIn);

		//Get our current time
		const double& getCurrentTime() const;

		//Set our current Time
		void setCurrentTime(const double& currentTimeIn);

		//Get our reduction factor
		const double& getReductionFactor() const;

		//Set reduction factor
		void setReductionFactor(const double& reductionFactorIn);
	};
}