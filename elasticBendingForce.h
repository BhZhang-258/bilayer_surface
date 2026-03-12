#ifndef ELASTICBENDINGFORCE_H
#define ELASTICBENDINGFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"
#include "GeometryUtils.h"

class elasticBendingForce
{
public:
	elasticBendingForce(elasticPlate &m_plate, timeStepper &m_stepper);
	~elasticBendingForce();

	void computeFb();
    void computeFb_old(int kkk);
	void computeJb();
    void setFirstJacobian();

private:
	elasticPlate *plate;
    timeStepper *stepper;

    VectorXd derivative_old;
    MatrixXd hessian_old;

};

#endif
