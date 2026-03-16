#ifndef ELASTICSTRETCHINGFORCE_H
#define ELASTICSTRETCHINGFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"
#include "GeometryUtils.h"
#include "ConstitutiveModel.h"

using namespace Geometry;
using namespace Constitutive;

class elasticStretchingForce
{
public:
	elasticStretchingForce(elasticPlate &m_plate, timeStepper &m_stepper);
	~elasticStretchingForce();
	void computeFs();

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
