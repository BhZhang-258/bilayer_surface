#ifndef ELASTICBENDINGFORCE_H
#define ELASTICBENDINGFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"
#include "GeometryUtils.h"
#include "ConstitutiveModel.h"

class elasticBendingForce
{
public:
	elasticBendingForce(elasticPlate &m_plate, timeStepper &m_stepper);
	~elasticBendingForce();

	void computeFb();
	void computeJb();
    void setFirstJacobian();

private:
	elasticPlate *plate;
    timeStepper *stepper;

    VectorXd derivative_old;
    MatrixXd hessian_old;

};

#endif
