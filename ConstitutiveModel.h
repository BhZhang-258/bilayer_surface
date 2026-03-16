#ifndef ConstitutiveModel_H
#define ConstitutiveModel_H

#include "eigenIncludes.h"
#include <optional>

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

using namespace Eigen;
namespace Constitutive
{

template<int DOF>
double stvk(
    const Eigen::Matrix4d &C,
    const Eigen::Matrix4d &T,
    const Eigen::Vector4d &A_,
    const Eigen::Vector4d &B_,
    const Eigen::Matrix<double,4,DOF> &gradA,
    const Eigen::Matrix<double,4,DOF> &gradB,
    const std::vector <Eigen::Matrix<double, DOF, DOF> > &hessA,
    const std::vector <Eigen::Matrix<double, DOF, DOF> > &hessB,
    Eigen::Matrix<double,DOF,1>* derivative,
    Eigen::Matrix<double,DOF,DOF>* hessian)
{

    Eigen::Vector4d eps1 = T * A_;
    Eigen::Vector4d eps2 = T * B_;

    double energy = 0.5 * eps2.transpose() * C * eps1;

    Eigen::MatrixXd B1 = T * gradA;
    Eigen::MatrixXd B2 = T * gradB;

    if (derivative)
    {
        *derivative =
            0.5 * B1.transpose() * C.transpose() * eps2 +
            0.5 * B2.transpose() * C.transpose() * eps1;
    }

    if (hessian)
    {
        Eigen::Matrix<double,DOF,DOF> Hmat =
            0.5 * B1.transpose() * C * B2 +
            0.5 * B2.transpose() * C * B1;

        Eigen::Matrix<double,DOF,DOF> Hgeo =
            Eigen::Matrix<double,DOF,DOF>::Zero();

        Eigen::RowVectorXd coeff1 = 0.5 * eps2.transpose() * C * T;
        Eigen::RowVectorXd coeff2 = 0.5 * eps1.transpose() * C * T;

        for(int i=0;i<4;i++)
        {
            Hgeo += coeff1(i) * hessA[i];
            Hgeo += coeff2(i) * hessB[i];
        }

        *hessian = Hmat + Hgeo;
    }

    return energy;
}





}

#endif