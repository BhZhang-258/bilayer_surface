#include "elasticStretchingForce.h"
#include <iostream>

elasticStretchingForce::elasticStretchingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
	stepper = &m_stepper;

	I3 << 1,0,0,
          0,1,0,
          0,0,1;

	Z3 << 0,0,0,
          0,0,0,
          0,0,0;
}

elasticStretchingForce::~elasticStretchingForce()
{
	;
}

void elasticStretchingForce::computeFs()
{
	for (int i = 0; i < plate->v_triangularElement.size(); i++)
	{
		

		VectorXi arrayIndex = plate->v_triangularElement[i].arrayIndex;

		Matrix2d A;

		Matrix<double, 4, 9> fFFderiv;
		std::vector <Eigen::Matrix<double, 9, 9> > fFFhess;

		A = plate->getFFF(i, &fFFderiv, &fFFhess);

	// for h 1
    	double thickness_1 = plate->thickness_1;
    	double coefficient_1 = thickness_1 / 4;

    	Matrix2d abar_1 = plate->v_triangularElement[i].abar_1;
    	Matrix2d abarInv_1 = plate->v_triangularElement[i].abarinv_1;

    	double dA_1 = 0.5 * sqrt( abar_1.determinant() );

    	double alpha_1, beta_1;
    	alpha_1 = plate->alpha_1;
    	beta_1 = plate->beta_1;

		
		Matrix4d C_1 = Matrix<double, 4, 4>::Zero();
		// C_1 << alpha_1 + 2 * beta_1,  alpha_1, 0,
		// 	   alpha_1, alpha_1 + 2 * beta_1,  0,
		// 		0,0,beta_1;
		C_1(0,0) = alpha_1 + 2 * beta_1;
		C_1(0,3) = alpha_1;
		C_1(3,0) = alpha_1;
		C_1(3,3) = alpha_1 + 2 * beta_1;
		C_1(1,2) = 2 * beta_1;
		C_1(2,1) = 2 * beta_1;

		
		Matrix<double, 9, 1> derivative_1;
		Matrix<double, 9, 9> hessian_1;
		

		Matrix2d M1 =  A - abar_1;
		Matrix4d S =  Matrix4d::Identity();
		Vector4d eps = Vector4d( M1(0,0), M1(0,1), M1(1,0), M1(1,1) );
		Matrix4d T_1 = Geometry::vecLeftMultiplyOperator(abarInv_1);
		
		// stvk energy/grad/hessian  for epsilon inner product <eps,eps> with C_1 as the stiffness matrix
		double StVK_11 = Constitutive::stvk<9>(C_1, T_1, 
			eps, eps,
			fFFderiv,fFFderiv, 
			fFFhess, fFFhess, 
			&derivative_1, &hessian_1);			

		derivative_1 = coefficient_1 * dA_1 * derivative_1;
		hessian_1 = coefficient_1 * dA_1 * hessian_1;

	// for h 2
		double thickness_2 = plate->thickness_2;
        double coefficient_2 = thickness_2 / 4;

        Matrix2d abar_2 = plate->v_triangularElement[i].abar_2;
        Matrix2d abarInv_2 = plate->v_triangularElement[i].abarinv_2;

        double dA_2 = 0.5 * sqrt( abar_2.determinant() );

        double alpha_2, beta_2;
        alpha_2 = plate->alpha_2;
        beta_2 = plate->beta_2;

		Matrix<double, 4, 4> C_2 = Matrix<double, 4, 4>::Zero();
		C_2(0,0) = alpha_2 + 2 * beta_2;
		C_2(0,3) = alpha_2;
		C_2(3,0) = alpha_2;
		C_2(3,3) = alpha_2 + 2 * beta_2;
		C_2(1,2) = 2 * beta_2;
		C_2(2,1) = 2 * beta_2;

		Matrix<double, 9, 1> derivative_2;
		Matrix<double, 9, 9> hessian_2;
		Matrix2d M2 =  A - abar_2;
		Vector4d eps2 = Vector4d( M2(0,0), M2(0,1), M2(1,0), M2(1,1) );
		Matrix4d T_2 = Geometry::vecLeftMultiplyOperator(abarInv_2);
		
		double StVK_2 = Constitutive::stvk<9>(C_2, T_2, 
			eps2, eps2,
			fFFderiv,fFFderiv, 
			fFFhess, fFFhess, 
			&derivative_2, &hessian_2);			
		
		derivative_2 = coefficient_2 * dA_2 * derivative_2;
		hessian_2 = coefficient_2 * dA_2 * hessian_2;
		

    	
		
		for (int j = 0; j < 9; j++)
        {
            stepper->addForce(arrayIndex(j), derivative_1(j) + derivative_2(j) );
        }
		

		for (int j = 0; j < 9; j++)
        {
            for (int k = 0; k < 9; k++)
            {
                stepper->addJacobian(arrayIndex(j), arrayIndex(k), hessian_1(j, k) + hessian_2(j, k) );
            }
        }
		
	}

	
}


void elasticStretchingForce::computeJs()
{
}

void elasticStretchingForce::setFirstJacobian()
{
	for (int i = 0; i < plate->v_triangularElement.size(); i++)
	{
		VectorXi arrayIndex = plate->v_triangularElement[i].arrayIndex;

		for (int j = 0; j < 9; j++)
        {
            for (int k = 0; k < 9; k++)
            {
                stepper->addJacobian(arrayIndex(j), arrayIndex(k), 1);
            }
        }
	}
}

