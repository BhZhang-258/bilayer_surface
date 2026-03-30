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
        Matrix<double, 4, 9> aderiv_t;
        std::vector<Matrix<double, 9, 9>> ahess_t;

        A = plate->getFFF(kkk, &aderiv_t, &ahess_t);

        aderiv.setZero();
        aderiv.block<4, 9>(0, 0) = aderiv_t;  
        
        ahess.resize(4);
        for (int i = 0; i < ahess_t.size(); i++)
        {
            ahess[i].setZero();
            ahess[i].block<9,9>(0,0) = ahess_t[i]; 
        }
        
        B = plate->getSFF( kkk, &bderiv, &bhess );

        for (int l = 0; l < 2; l++)
		{
			computeFb_layer(kkk, l);
		}

        
    }


}

void elasticBendingForce::computeFb_layer(int idx, int l)
{
        VectorXi arrayNum = plate->mesh.Fbenddof.row(idx);
        Layer layer = plate->v_triangularElement[idx].layers[l];
        double thickness = layer.thickness;
       
        double coeff = pow(thickness, 3) / 12;
        double crossTermCoeff = 0;
       
        if (l == 0)
        {
            crossTermCoeff = thickness * thickness / 4;
        }
        else if (l == 1)
        {
            crossTermCoeff = - thickness * thickness / 4;
        }
        else 
        {
            cerr<< "warning, invalid layer index" <<endl;
            exit(1);
        }

		Matrix2d Fg = layer.Fg;

		Matrix2d G = Fg.transpose() * Fg;
		Matrix2d Ginv = G.inverse();
		
		
		Matrix2d M = layer.M;
		Matrix2d Minv = layer.Minv;
		
		Matrix2d abar_g = M.transpose() * G * M;
		Matrix2d abarinv = abar_g.inverse();

		Matrix4d C = layer.C;

		Matrix4d kronc = Geometry::kron( Minv.transpose(), M * abarinv);


    	double dA = 0.5 * sqrt( abar_g.determinant() );

        Matrix2d Epsilon_bend = (B - plate->v_triangularElement[idx].bbar);
        Matrix2d Epsilon_stretch = (A - abar_g);

        Vector4d epsB = Eigen::Map<Vector4d>(Epsilon_bend.data());
        Vector4d epsA = Eigen::Map<Vector4d>(Epsilon_stretch.data());


        Matrix<double, 18, 1> derivative;
        derivative.setZero(18, 1);

        Matrix<double, 18, 18> hessian;
        hessian.setZero(18, 18);

        Matrix<double, 18, 1> derivative_t;
        Matrix<double, 18, 18> hessian_t;

        double StVK_b1 = Constitutive::stvk<18>(C, kronc, 
			epsB, epsB,
			bderiv,bderiv, 
			bhess, bhess, 
			&derivative_t, &hessian_t);	
        		
        derivative += coeff * dA * derivative_t;
        hessian += coeff * dA * hessian_t;
        
        double StVK_c1 = Constitutive::stvk<18>(C, kronc, 
			epsA, epsB,
			aderiv,bderiv, 
			ahess, bhess, 
			&derivative_t, &hessian_t);
        derivative += crossTermCoeff * dA * derivative_t;
        hessian += crossTermCoeff * dA * hessian_t;


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
