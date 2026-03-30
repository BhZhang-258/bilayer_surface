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
	void computeFb_layer(int idx , int l);

	void computeJb();
    void setFirstJacobian();

private:
	elasticPlate *plate;
    timeStepper *stepper;

    VectorXd derivative_old;
    MatrixXd hessian_old;

	Matrix2d A;
	Matrix2d B;
	Matrix<double, 4, 18> aderiv;
    std::vector<Matrix<double, 18, 18>> ahess;

	Matrix<double, 4, 18 > bderiv;
    std::vector<Matrix<double, 18 , 18 >> bhess;


};

#endif
