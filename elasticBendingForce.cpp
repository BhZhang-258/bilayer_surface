#include "elasticBendingForce.h"

elasticBendingForce::elasticBendingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
    stepper = &m_stepper;
}

elasticBendingForce::~elasticBendingForce()
{
	;
}

void elasticBendingForce::computeFb()
{
    for (int kkk = 0; kkk < plate->triangular.size(); kkk++)
    {
        double thickness_1 = plate->thickness_1;
        Matrix2d abars_1 = plate->v_triangularElement[kkk].abar_1;
        Matrix2d bbars_1 = plate->v_triangularElement[kkk].bbar_1;
        double coeff_1 = pow(thickness_1, 3) / 24;
        Matrix2d abarinv_1 = plate->v_triangularElement[kkk].abarinv_1;
        double alpha_1 = plate->alpha_1;
        double beta_1 = plate->beta_1;
        double dA_1 = 0.5 * sqrt(abars_1.determinant());


        double thickness_2 = plate->thickness_2;
        Matrix2d abars_2 = plate->v_triangularElement[kkk].abar_2;
        Matrix2d bbars_2 = plate->v_triangularElement[kkk].bbar_2;
        double coeff_2 = pow(thickness_2, 3) / 24;
        Matrix2d abarinv_2 = plate->v_triangularElement[kkk].abarinv_2;
        double alpha_2 = plate->alpha_2;
        double beta_2 = plate->beta_2;
        double dA_2 = 0.5 * sqrt(abars_2.determinant());


        Matrix<double, 4, 9> aderiv_t;
        std::vector<Matrix<double, 9, 9>> ahess_t;

        Matrix2d a = plate->getFFF(kkk, &aderiv_t, &ahess_t);

        Matrix<double, 4, 18> aderiv = Matrix<double, 4, 18>::Zero();
        aderiv.block<4, 9>(0, 0) = aderiv_t;  
        
        std::vector<Matrix<double, 18, 18>> ahess(4, Matrix<double, 18, 18>::Zero());
        for (int i = 0; i < ahess_t.size(); i++)
        {
            ahess[i].block<9,9>(0,0) = ahess_t[i]; 
        }
        



        Matrix<double, 4, 18 > bderiv;
        std::vector<Matrix<double, 18 , 18 > > bhess;
        Matrix2d b = plate->getSFF( kkk, &bderiv, &bhess );


        double crossTermCoeff_1 = thickness_1 * thickness_1 / 8;
       
        double crossTermCoeff_2 = - thickness_2 * thickness_2 / 8;


        Matrix2d B_1 = (b - bbars_1);
        Matrix2d B_2 = (b - bbars_2);
        Matrix2d A_1 = (a - abars_1);
        Matrix2d A_2 = (a - abars_2);

        Vector4d epsB_1 = Vector4d( B_1(0,0), B_1(0,1), B_1(1,0), B_1(1,1) );
        Vector4d epsB_2 = Vector4d( B_2(0,0), B_2(0,1), B_2(1,0), B_2(1,1) );
        Vector4d epsA_1 = Vector4d( A_1(0,0), A_1(0,1), A_1(1,0), A_1(1,1) );
        Vector4d epsA_2 = Vector4d( A_2(0,0), A_2(0,1), A_2(1,0), A_2(1,1) );     

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

        Matrix<double, 4, 4> C_2 = Matrix<double, 4, 4>::Zero();
		C_2(0,0) = alpha_2 + 2 * beta_2;
		C_2(0,3) = alpha_2;
		C_2(3,0) = alpha_2;
		C_2(3,3) = alpha_2 + 2 * beta_2;
		C_2(1,2) = 2 * beta_2;
		C_2(2,1) = 2 * beta_2;

        Matrix4d T_1 = Geometry::vecLeftMultiplyOperator(abarinv_1);
        Matrix4d T_2 = Geometry::vecLeftMultiplyOperator(abarinv_2);

        Matrix<double, 18, 1> derivative;
        derivative.setZero(18, 1);

        Matrix<double, 18, 18> hessian;
        hessian.setZero(18, 18);

        Matrix<double, 18, 1> derivative_t;
        derivative.setZero(18, 1);
        Matrix<double, 18, 18> hessian_t;
        hessian.setZero(18, 18);

        double StVK_b1 = Constitutive::stvk<18>(C_1, T_1, 
			epsB_1, epsB_1,
			bderiv,bderiv, 
			bhess, bhess, 
			&derivative_t, &hessian_t);	
        		
        derivative += coeff_1 * dA_1 * derivative_t;
        hessian += coeff_1 * dA_1 * hessian_t;
        
        double StVK_c1 = Constitutive::stvk<18>(C_1, T_1, 
			epsA_1, epsB_1,
			aderiv,bderiv, 
			ahess, bhess, 
			&derivative_t, &hessian_t);
        derivative += crossTermCoeff_1 * dA_1 * derivative_t;
        hessian += crossTermCoeff_1 * dA_1 * hessian_t;

        double StVK_b2 = Constitutive::stvk<18>(C_2, T_2, 
			epsB_2, epsB_2,
			bderiv,bderiv, 
			bhess, bhess, 
			&derivative_t, &hessian_t);	
        		
        derivative += coeff_2 * dA_2 * derivative_t;
        hessian += coeff_2 * dA_2 * hessian_t;

        double StVK_c2 = Constitutive::stvk<18>(C_2, T_2, 
			epsA_2, epsB_2,
			aderiv,bderiv, 
			ahess, bhess, 
			&derivative_t, &hessian_t);
        derivative += crossTermCoeff_2 * dA_2 * derivative_t;
        hessian += crossTermCoeff_2 * dA_2 * hessian_t;            

        
        VectorXi arrayNum = plate->mesh.Fbenddof.row(kkk);

        for (int j = 0; j < 18; j++)
        {
            if ( arrayNum(j) >= 0 )
            {
                stepper->addForce(arrayNum(j), derivative(j) );
            }

        }

        for (int j = 0; j < 18; j++)
        {
            for (int k = 0; k < 18; k++)
            {
                if ( arrayNum(j) >= 0 && arrayNum(k) >= 0 )
                {
                    stepper->addJacobian(arrayNum(j), arrayNum(k), hessian(j,k) );
                }
                
            }
        }

    }


}
void elasticBendingForce::computeJb()
{
    ;
}

void elasticBendingForce::setFirstJacobian()
{
    for (int i = 0; i < plate->triangular.size(); i++)
    {
        VectorXi arrayNum = plate->mesh.Fbenddof.row(i);

        for (int j = 0; j < 18; j++)
        {
            for (int k = 0; k < 18; k++)
            {

                if ( arrayNum(j) >= 0 && arrayNum(k) >= 0 )
                {
                    stepper->addJacobian(arrayNum(j), arrayNum(k), 1);
                }
                
            }
        }
    }
}
