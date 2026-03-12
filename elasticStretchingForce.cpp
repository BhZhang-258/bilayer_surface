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


		VectorXd gradA1, gradA2, gradA3, gradA4;
		MatrixXd hessA1, hessA2, hessA3, hessA4;


    	gradA1 = fFFderiv.row(0).transpose();
		gradA2 = fFFderiv.row(1).transpose();
		gradA3 = fFFderiv.row(2).transpose();
		gradA4 = fFFderiv.row(3).transpose();

    	hessA1 = fFFhess[0];
    	hessA2 = fFFhess[1];
    	hessA3 = fFFhess[2];
    	hessA4 = fFFhess[3];



    	double thickness_1 = plate->thickness_1;
    	double coefficient_1 = thickness_1 / 4;

    	Matrix2d abar_1 = plate->v_triangularElement[i].abar_1;
    	Matrix2d abarInv_1 = plate->v_triangularElement[i].abarinv_1;

    	double dA_1 = 0.5 * sqrt( abar_1.determinant() );

    	double alpha_1, beta_1;
    	alpha_1 = plate->alpha_1;
    	beta_1 = plate->beta_1;

    	Matrix2d M_1 = abarInv_1 * (A - abar_1);

    	//double StVK_1 = 0.5 * alpha_1 * ( M_1.trace() ) * ( M_1.trace() ) + beta_1 * (M_1 * M_1).trace();

	    Matrix2d tempMatrix_1;
    	tempMatrix_1 = alpha_1 * M_1.trace() * abarInv_1 + 2 * beta_1 * M_1 * abarInv_1;

    	VectorXd derivative_1;

    	derivative_1 = VectorXd::Zero(9);

    	derivative_1 = dA_1 * coefficient_1 * (gradA1 * tempMatrix_1(0,0) + gradA2 * tempMatrix_1(0,1) + gradA3 * tempMatrix_1(1,0) + gradA4 * tempMatrix_1(1,1) );  

    	VectorXd inner_1;

    	inner_1 = (gradA1 * abarInv_1(0,0) + gradA2 * abarInv_1(0,1) + gradA3 * abarInv_1(1,0) + gradA4 * abarInv_1(1,1) );

    	MatrixXd hessian_1;

    	hessian_1 = MatrixXd::Zero(9, 9);

    	hessian_1 = (coefficient_1 * dA_1) * alpha_1 * (inner_1 * inner_1.transpose());

    	Matrix2d MainV_1 = M_1 * abarInv_1;

    	double coeff_i_1;
    
    	coeff_i_1 = alpha_1 * M_1.trace() * abarInv_1(0, 0) + 2 * beta_1 * MainV_1(0, 0);
    	hessian_1 = hessian_1 + ( coefficient_1 * dA_1 ) * coeff_i_1 * hessA1;

    	coeff_i_1 = alpha_1 * M_1.trace() * abarInv_1(0, 1) + 2 * beta_1 * MainV_1(0, 1);
    	hessian_1 = hessian_1 + ( coefficient_1 * dA_1 ) * coeff_i_1 * hessA2;

    	coeff_i_1 = alpha_1 * M_1.trace() * abarInv_1(1, 0) + 2 * beta_1 * MainV_1(1, 0);
    	hessian_1 = hessian_1 + ( coefficient_1 * dA_1 ) * coeff_i_1 * hessA3;

    	coeff_i_1 = alpha_1 * M_1.trace() * abarInv_1(1, 1) + 2 * beta_1 * MainV_1(1, 1);
    	hessian_1 = hessian_1 + ( coefficient_1 * dA_1 ) * coeff_i_1 * hessA4;

    	VectorXd inner001 = abarInv_1(0,0) * gradA1 + abarInv_1(0,1) * gradA3; 
		VectorXd inner011 = abarInv_1(0,0) * gradA2 + abarInv_1(0,1) * gradA4;
		VectorXd inner101 = abarInv_1(1,0) * gradA1 + abarInv_1(1,1) * gradA3;
		VectorXd inner111 = abarInv_1(1,0) * gradA2 + abarInv_1(1,1) * gradA4;

		hessian_1 = hessian_1 + (coefficient_1 * dA_1) * 2 * beta_1 * (inner001 * inner001.transpose());
		hessian_1 = hessian_1 + (coefficient_1 * dA_1) * 2 * beta_1 * (inner011 * inner101.transpose() + inner101 * inner011.transpose());
		hessian_1 = hessian_1 + (coefficient_1 * dA_1) * 2 * beta_1 * (inner111 * inner111.transpose());



        double thickness_2 = plate->thickness_2;
        double coefficient_2 = thickness_2 / 4;

        Matrix2d abar_2 = plate->v_triangularElement[i].abar_2;
        Matrix2d abarInv_2 = plate->v_triangularElement[i].abarinv_2;

        double dA_2 = 0.5 * sqrt( abar_2.determinant() );

        double alpha_2, beta_2;
        alpha_2 = plate->alpha_2;
        beta_2 = plate->beta_2;

        Matrix2d M_2 = abarInv_2 * (A - abar_2);

        //double StVK_2 = 0.5 * alpha_2 * ( M_2.trace() ) * ( M_2.trace() ) + beta_2 * (M_2 * M_2).trace();

        Matrix2d tempMatrix_2;
        tempMatrix_2 = alpha_2 * M_2.trace() * abarInv_2 + 2 * beta_2 * M_2 * abarInv_2;

        VectorXd derivative_2;

        derivative_2 = VectorXd::Zero(9);

        derivative_2 = dA_2 * coefficient_2 * (gradA1 * tempMatrix_2(0,0) + gradA2 * tempMatrix_2(0,1) + gradA3 * tempMatrix_2(1,0) + gradA4 * tempMatrix_2(1,1) );  

        VectorXd inner_2;

        inner_2 = (gradA1 * abarInv_2(0,0) + gradA2 * abarInv_2(0,1) + gradA3 * abarInv_2(1,0) + gradA4 * abarInv_2(1,1) );

        MatrixXd hessian_2;

        hessian_2 = MatrixXd::Zero(9, 9);

        hessian_2 = (coefficient_2 * dA_2) * alpha_2 * (inner_2 * inner_2.transpose());

        Matrix2d MainV_2 = M_2 * abarInv_2;

        double coeff_i_2;
    
        coeff_i_2 = alpha_2 * M_2.trace() * abarInv_2(0, 0) + 2 * beta_2 * MainV_2(0, 0);
        hessian_2 = hessian_2 + ( coefficient_2 * dA_2 ) * coeff_i_2 * hessA1;

        coeff_i_2 = alpha_2 * M_2.trace() * abarInv_2(0, 1) + 2 * beta_2 * MainV_2(0, 1);
        hessian_2 = hessian_2 + ( coefficient_2 * dA_2 ) * coeff_i_2 * hessA2;

        coeff_i_2 = alpha_2 * M_2.trace() * abarInv_2(1, 0) + 2 * beta_2 * MainV_2(1, 0);
        hessian_2 = hessian_2 + ( coefficient_2 * dA_2 ) * coeff_i_2 * hessA3;

        coeff_i_2 = alpha_2 * M_2.trace() * abarInv_2(1, 1) + 2 * beta_2 * MainV_2(1, 1);
        hessian_2 = hessian_2 + ( coefficient_2 * dA_2 ) * coeff_i_2 * hessA4;

        VectorXd inner002 = abarInv_2(0,0) * gradA1 + abarInv_2(0,1) * gradA3; 
        VectorXd inner012 = abarInv_2(0,0) * gradA2 + abarInv_2(0,1) * gradA4;
        VectorXd inner102 = abarInv_2(1,0) * gradA1 + abarInv_2(1,1) * gradA3;
        VectorXd inner112 = abarInv_2(1,0) * gradA2 + abarInv_2(1,1) * gradA4;

        hessian_2 = hessian_2 + (coefficient_2 * dA_2) * 2 * beta_2 * (inner002 * inner002.transpose());
        hessian_2 = hessian_2 + (coefficient_2 * dA_2) * 2 * beta_2 * (inner012 * inner102.transpose() + inner102 * inner012.transpose());
        hessian_2 = hessian_2 + (coefficient_2 * dA_2) * 2 * beta_2 * (inner112 * inner112.transpose());

		
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

