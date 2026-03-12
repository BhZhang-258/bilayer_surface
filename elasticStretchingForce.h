#ifndef ELASTICSTRETCHINGFORCE_H
#define ELASTICSTRETCHINGFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"
#include "GeometryUtils.h"

using namespace Geometry;

class elasticStretchingForce
{
public:
	elasticStretchingForce(elasticPlate &m_plate, timeStepper &m_stepper);
	~elasticStretchingForce();
	void computeFs();
    void computeFs_old(int i);

	void computeJs();
    void setFirstJacobian();

private:
	elasticPlate *plate;
    timeStepper *stepper;

    Matrix3d I3;
    Matrix3d Z3;

    VectorXd derivative_1_old;
    VectorXd derivative_2_old;
    MatrixXd hessian_1_old;
    MatrixXd hessian_2_old;


};

#endif
