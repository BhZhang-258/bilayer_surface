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


        Matrix<double, 4, 9 > aderiv;
        std::vector<Matrix<double, 9 , 9 > > ahess;
        Matrix2d a = plate->getFFF( kkk, &aderiv, &ahess );


        Matrix<double, 4, 18 > bderiv;
        std::vector<Matrix<double, 18 , 18 > > bhess;
        Matrix2d b = plate->getSFF( kkk, &bderiv, &bhess );


        double crossTermCoeff_1 = thickness_1 * thickness_1 / 8;
        Matrix2d sigma_1 = abarinv_1 * (a - abars_1);

        double crossTermCoeff_2 = - thickness_2 * thickness_2 / 8;
        Matrix2d sigma_2 = abarinv_2 * (a - abars_2);


        Matrix2d M_1 = abarinv_1 * (b - bbars_1);
        Matrix2d M_2 = abarinv_2 * (b - bbars_2);

        
        VectorXd derivative;
        derivative.setZero(18, 1);

        Matrix2d temp_1 = 0.5 * alpha_1 * M_1.trace() * abarinv_1 + beta_1 * M_1 * abarinv_1;
        Matrix2d temp_2 = 0.5 * alpha_2 * M_2.trace() * abarinv_2 + beta_2 * M_2 * abarinv_2;

        derivative =              2 * coeff_1 * dA_1 * bderiv.transpose() * Map<Vector4d>(temp_1.data());
        derivative = derivative + 2 * coeff_2 * dA_2 * bderiv.transpose() * Map<Vector4d>(temp_2.data());


        VectorXd crossDeriv_1 = crossTermCoeff_1 * dA_1 * ( aderiv.transpose() * Map<Vector4d>(temp_1.data()) );
        VectorXd crossDeriv_2 = crossTermCoeff_2 * dA_2 * ( aderiv.transpose() * Map<Vector4d>(temp_2.data()) );

        derivative.segment(0,9) = derivative.segment(0,9) + crossDeriv_1 + crossDeriv_2;

        Matrix2d temp_3 = 0.5 * alpha_1 * sigma_1.trace() * abarinv_1 + beta_1 * sigma_1 * abarinv_1;
        Matrix2d temp_4 = 0.5 * alpha_2 * sigma_2.trace() * abarinv_2 + beta_2 * sigma_2 * abarinv_2;

        derivative = derivative + crossTermCoeff_1 * dA_1 * bderiv.transpose() * Map<Vector4d>(temp_3.data());
        derivative = derivative + crossTermCoeff_2 * dA_2 * bderiv.transpose() * Map<Vector4d>(temp_4.data());

        MatrixXd hessian;
        hessian.setZero(18, 18);

        Matrix<double, 1, 18 > inner_1 = bderiv.transpose() * Map<Vector4d>(abarinv_1.data());
        Matrix<double, 1, 18 > inner_2 = bderiv.transpose() * Map<Vector4d>(abarinv_2.data());

        Matrix<double, 1, 18 > inner_3;
        Matrix<double, 1, 18 > inner_4;
        inner_3.setZero(1, 18);
        inner_4.setZero(1, 18);
        inner_3.segment(0,9) = aderiv.transpose() * Map<Vector4d>(abarinv_1.data());
        inner_4.segment(0,9) = aderiv.transpose() * Map<Vector4d>(abarinv_2.data());

        
        hessian  = coeff_1 * dA_1 * alpha_1 * inner_1.transpose() * inner_1;
        hessian += coeff_2 * dA_2 * alpha_2 * inner_2.transpose() * inner_2;

        hessian += crossTermCoeff_1 * dA_1 * 0.5 * alpha_1 * ( (inner_1.transpose() * inner_3) + (inner_3.transpose() * inner_1) );
        hessian += crossTermCoeff_2 * dA_2 * 0.5 * alpha_2 * ( (inner_2.transpose() * inner_4) + (inner_4.transpose() * inner_2) );


        Matrix2d Mainv_1 = M_1 * abarinv_1;
        Matrix2d Mainv_2 = M_2 * abarinv_2;

        Matrix2d Sainv_1 = sigma_1 * abarinv_1;
        Matrix2d Sainv_2 = sigma_2 * abarinv_2;

        for (int i = 0; i < 4; i++) 
        {
            hessian += coeff_1 * dA_1 * ( 2 * temp_1(i) ) * bhess[i];
            hessian += coeff_2 * dA_2 * ( 2 * temp_2(i) ) * bhess[i];

            hessian.block(0,0,9,9) = hessian.block(0,0,9,9) + crossTermCoeff_1 * dA_1 * ( temp_1(i) ) * ahess[i];
            hessian.block(0,0,9,9) = hessian.block(0,0,9,9) + crossTermCoeff_2 * dA_2 * ( temp_2(i) ) * ahess[i];
            hessian += crossTermCoeff_1 * dA_1 * ( temp_3(i) ) * bhess[i];
            hessian += crossTermCoeff_2 * dA_2 * ( temp_4(i) ) * bhess[i];
        }

        
        Matrix<double, 1, 18> inner001 = abarinv_1(0, 0) * bderiv.row(0) + abarinv_1(0, 1) * bderiv.row(2);
        Matrix<double, 1, 18> inner011 = abarinv_1(0, 0) * bderiv.row(1) + abarinv_1(0, 1) * bderiv.row(3);
        Matrix<double, 1, 18> inner101 = abarinv_1(1, 0) * bderiv.row(0) + abarinv_1(1, 1) * bderiv.row(2);
        Matrix<double, 1, 18> inner111 = abarinv_1(1, 0) * bderiv.row(1) + abarinv_1(1, 1) * bderiv.row(3);
        hessian += 2 * coeff_1 * dA_1 * beta_1 *  inner001.transpose() * inner001;
        hessian += 2 * coeff_1 * dA_1 * beta_1 * (inner011.transpose() * inner101 + inner101.transpose() * inner011);
        hessian += 2 * coeff_1 * dA_1 * beta_1 *  inner111.transpose() * inner111;

        Matrix<double, 1, 18> inner002 = abarinv_2(0, 0) * bderiv.row(0) + abarinv_2(0, 1) * bderiv.row(2);
        Matrix<double, 1, 18> inner012 = abarinv_2(0, 0) * bderiv.row(1) + abarinv_2(0, 1) * bderiv.row(3);
        Matrix<double, 1, 18> inner102 = abarinv_2(1, 0) * bderiv.row(0) + abarinv_2(1, 1) * bderiv.row(2);
        Matrix<double, 1, 18> inner112 = abarinv_2(1, 0) * bderiv.row(1) + abarinv_2(1, 1) * bderiv.row(3);
        hessian += 2 * coeff_2 * dA_2 * beta_2 *  inner002.transpose() * inner002;
        hessian += 2 * coeff_2 * dA_2 * beta_2 * (inner012.transpose() * inner102 + inner102.transpose() * inner012);
        hessian += 2 * coeff_2 * dA_2 * beta_2 *  inner112.transpose() * inner112;

        Matrix<double, 1, 18> inner003;
        inner003.setZero(1, 18);
        inner003.segment(0,9) = abarinv_1(0, 0) * aderiv.row(0) + abarinv_1(0, 1) * aderiv.row(2);
        Matrix<double, 1, 18> inner013;
        inner013.setZero(1, 18);
        inner013.segment(0,9) = abarinv_1(0, 0) * aderiv.row(1) + abarinv_1(0, 1) * aderiv.row(3);
        Matrix<double, 1, 18> inner103;
        inner103.setZero(1, 18);
        inner103.segment(0,9) = abarinv_1(1, 0) * aderiv.row(0) + abarinv_1(1, 1) * aderiv.row(2);
        Matrix<double, 1, 18> inner113;
        inner113.setZero(1, 18);
        inner113.segment(0,9) = abarinv_1(1, 0) * aderiv.row(1) + abarinv_1(1, 1) * aderiv.row(3);

        hessian += crossTermCoeff_1 * dA_1 * beta_1 *  inner001.transpose() * inner003;
        hessian += crossTermCoeff_1 * dA_1 * beta_1 *  inner003.transpose() * inner001;
        hessian += crossTermCoeff_1 * dA_1 * beta_1 * (inner011.transpose() * inner103 + inner103.transpose() * inner011);
        hessian += crossTermCoeff_1 * dA_1 * beta_1 * (inner013.transpose() * inner101 + inner101.transpose() * inner013);
        hessian += crossTermCoeff_1 * dA_1 * beta_1 *  inner111.transpose() * inner113;
        hessian += crossTermCoeff_1 * dA_1 * beta_1 *  inner113.transpose() * inner111;

        Matrix<double, 1, 18> inner004;
        inner004.setZero(1, 18);
        inner004.segment(0,9) = abarinv_2(0, 0) * aderiv.row(0) + abarinv_2(0, 1) * aderiv.row(2);
        Matrix<double, 1, 18> inner014;
        inner014.setZero(1, 18);
        inner014.segment(0,9) = abarinv_2(0, 0) * aderiv.row(1) + abarinv_2(0, 1) * aderiv.row(3);
        Matrix<double, 1, 18> inner104;
        inner104.setZero(1, 18);
        inner104.segment(0,9) = abarinv_2(1, 0) * aderiv.row(0) + abarinv_2(1, 1) * aderiv.row(2);
        Matrix<double, 1, 18> inner114;
        inner114.setZero(1, 18);
        inner114.segment(0,9) = abarinv_2(1, 0) * aderiv.row(1) + abarinv_2(1, 1) * aderiv.row(3);

        hessian += crossTermCoeff_2 * dA_2 * beta_2 *  inner002.transpose() * inner004;
        hessian += crossTermCoeff_2 * dA_2 * beta_2 *  inner004.transpose() * inner002;
        hessian += crossTermCoeff_2 * dA_2 * beta_2 * (inner012.transpose() * inner104 + inner104.transpose() * inner012);
        hessian += crossTermCoeff_2 * dA_2 * beta_2 * (inner014.transpose() * inner102 + inner102.transpose() * inner014);
        hessian += crossTermCoeff_2 * dA_2 * beta_2 *  inner112.transpose() * inner114;
        hessian += crossTermCoeff_2 * dA_2 * beta_2 *  inner114.transpose() * inner112;


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
