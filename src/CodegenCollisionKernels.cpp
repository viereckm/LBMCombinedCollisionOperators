#include "CodegenCollisionKernels.hpp"

namespace walberla {
namespace lbm {
namespace codegenKernels {

void D2Q9SRTKernel( int64_t const ctr_0, int64_t const ctr_1, double * RESTRICT _data_density, double * RESTRICT const _data_pdfs, double * RESTRICT const _data_pdfs_tmp, double * RESTRICT _data_velocity, int64_t const _stride_density_0, int64_t const _stride_density_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_tmp_0, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, double omega )
{
    const double xi_4 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 3*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 7*_stride_pdfs_2];
    const double xi_6 = -_data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 6*_stride_pdfs_2];
    const double vel0Term = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 4*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 8*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 6*_stride_pdfs_2];
    const double xi_5 = vel0Term - xi_4 - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 5*_stride_pdfs_2];
    const double vel1Term = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 5*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2];
    const double xi_7 = vel1Term - xi_6 - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 7*_stride_pdfs_2] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 2*_stride_pdfs_2] - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 8*_stride_pdfs_2];
    const double delta_rho = vel0Term + vel1Term + xi_4 + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 2*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1];
    const double rho = delta_rho + 1.0;
    const double xi_9 = rho*1.5;
    const double xi_2 = ((1.0) / (rho));
    const double u_0 = xi_2*xi_5;
    const double xi_8 = u_0*u_0;
    const double xi_14 = u_0*0.33333333333333331;
    const double xi_15 = xi_8*0.5;
    const double u_1 = xi_2*xi_7;
    const double xi_10 = u_1*u_1;
    const double xi_12 = u_1*0.33333333333333331;
    const double xi_13 = xi_10*0.5;
    const double momdensity_0 = xi_5;
    const double momdensity_1 = xi_7;
    const double xi_0 = momdensity_0*xi_2;
    const double xi_1 = momdensity_1*xi_2;
    const double u0Mu1 = u_0 - u_1;
    const double xi_17 = u0Mu1*0.083333333333333329;
    const double xi_18 = 0.125*(u0Mu1*u0Mu1);
    const double u0Pu1 = u_0 + u_1;
    const double xi_19 = u0Pu1*0.083333333333333329;
    const double xi_20 = 0.125*(u0Pu1*u0Pu1);
    const double f_eq_common = delta_rho - xi_10*xi_9 - xi_8*xi_9;
    const double xi_11 = f_eq_common*0.1111111111111111;
    const double xi_16 = f_eq_common*0.027777777777777776;
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1] = omega*(f_eq_common*0.44444444444444442 - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2] = omega*(rho*(xi_12 + xi_13) + xi_11 - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 2*_stride_pdfs_tmp_2] = omega*(rho*(-xi_12 + xi_13) + xi_11 - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 2*_stride_pdfs_2]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 2*_stride_pdfs_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 3*_stride_pdfs_tmp_2] = omega*(rho*(-xi_14 + xi_15) + xi_11 - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 3*_stride_pdfs_2]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 3*_stride_pdfs_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 4*_stride_pdfs_tmp_2] = omega*(rho*(xi_14 + xi_15) + xi_11 - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 4*_stride_pdfs_2]) + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 4*_stride_pdfs_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 5*_stride_pdfs_tmp_2] = omega*(rho*(-xi_17 + xi_18) + xi_16 - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 5*_stride_pdfs_2]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 5*_stride_pdfs_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 6*_stride_pdfs_tmp_2] = omega*(rho*(xi_19 + xi_20) + xi_16 + xi_6) + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 6*_stride_pdfs_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 7*_stride_pdfs_tmp_2] = omega*(rho*(-xi_19 + xi_20) + xi_16 - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 7*_stride_pdfs_2]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 7*_stride_pdfs_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 8*_stride_pdfs_tmp_2] = omega*(rho*(xi_17 + xi_18) + xi_16 - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 8*_stride_pdfs_2]) + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 8*_stride_pdfs_2];
    _data_density[_stride_density_0*ctr_0 + _stride_density_1*ctr_1] = rho;
    _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1] = xi_0;
    _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2] = xi_1;
}

void D2Q9TRTKernel( int64_t const ctr_0, int64_t const ctr_1, double * RESTRICT _data_density, double * RESTRICT const _data_pdfs, double * RESTRICT const _data_pdfs_tmp, double * RESTRICT _data_velocity, int64_t const _stride_density_0, int64_t const _stride_density_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_tmp_0, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, double omega )
{
    const double xi_2 = ((1.0) / (omega*-0.25 + 2.0));
    const double rr_0 = xi_2*(omega*-2.0 + 4.0);
   
    const double xi_5 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 3*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 7*_stride_pdfs_2];
    const double xi_11 = 0.5*_data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2];
    const double xi_12 = 0.5*_data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 2*_stride_pdfs_2];
    const double xi_17 = 0.5*_data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 3*_stride_pdfs_2];
    const double xi_18 = 0.5*_data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 4*_stride_pdfs_2];
    const double xi_21 = 0.5*_data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 5*_stride_pdfs_2];
    const double xi_22 = 0.5*_data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 8*_stride_pdfs_2];
    const double xi_27 = 0.5*_data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 6*_stride_pdfs_2];
    const double xi_28 = 0.5*_data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 7*_stride_pdfs_2];
    const double vel0Term = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 4*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 8*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 6*_stride_pdfs_2];
    const double xi_6 = vel0Term - xi_5 - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 5*_stride_pdfs_2];
    const double vel1Term = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 5*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2];
    const double xi_7 = vel1Term - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 7*_stride_pdfs_2] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 2*_stride_pdfs_2] - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 8*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 6*_stride_pdfs_2];
    const double delta_rho = vel0Term + vel1Term + xi_5 + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 2*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1];
    const double rho = delta_rho + 1.0;
    const double xi_9 = rho*1.5;
    const double xi_13 = rho*0.33333333333333331;
    const double xi_23 = rho*0.083333333333333329;
    const double xi_3 = ((1.0) / (rho));
    const double u_0 = xi_3*xi_6;
    const double xi_8 = u_0*u_0;
    const double xi_19 = u_0*xi_13 + xi_17 - xi_18;
    const double u_1 = xi_3*xi_7;
    const double xi_10 = u_1*u_1;
    const double xi_14 = u_1*xi_13 - xi_11 + xi_12;
    const double momdensity_0 = xi_6;
    const double momdensity_1 = xi_7;
    const double xi_0 = momdensity_0*xi_3;
    const double xi_1 = momdensity_1*xi_3;
    const double u0Mu1 = u_0 - u_1;
    const double xi_24 = u0Mu1*xi_23 + xi_21 - xi_22;
    const double u0Pu1 = u_0 + u_1;
    const double xi_29 = u0Pu1*xi_23 - xi_27 + xi_28;
    const double f_eq_common = delta_rho - xi_10*xi_9 - xi_8*xi_9;
    const double xi_15 = f_eq_common*-0.1111111111111111;
    const double xi_16 = omega*(rho*xi_10*0.5 - xi_11 - xi_12 - xi_15);
    const double xi_20 = omega*(rho*xi_8*0.5 - xi_15 - xi_17 - xi_18);
    const double xi_25 = f_eq_common*-0.027777777777777776;
    const double xi_26 = omega*(rho*0.125*(u0Mu1*u0Mu1) - xi_21 - xi_22 - xi_25);
    const double xi_30 = omega*(rho*0.125*(u0Pu1*u0Pu1) - xi_25 - xi_27 - xi_28);
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1] = omega*(f_eq_common*0.44444444444444442 - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2] = rr_0*xi_14 + xi_16 + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 2*_stride_pdfs_tmp_2] = -rr_0*xi_14 + xi_16 + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 2*_stride_pdfs_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 3*_stride_pdfs_tmp_2] = -rr_0*xi_19 + xi_20 + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 3*_stride_pdfs_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 4*_stride_pdfs_tmp_2] = rr_0*xi_19 + xi_20 + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 4*_stride_pdfs_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 5*_stride_pdfs_tmp_2] = -rr_0*xi_24 + xi_26 + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 5*_stride_pdfs_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 6*_stride_pdfs_tmp_2] = rr_0*xi_29 + xi_30 + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 6*_stride_pdfs_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 7*_stride_pdfs_tmp_2] = -rr_0*xi_29 + xi_30 + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 7*_stride_pdfs_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 8*_stride_pdfs_tmp_2] = rr_0*xi_24 + xi_26 + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 8*_stride_pdfs_2];
    _data_density[_stride_density_0*ctr_0 + _stride_density_1*ctr_1] = rho;
    _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1] = xi_0;
    _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2] = xi_1;
}

void D2Q9CumulantKernel( int64_t const ctr_0, int64_t const ctr_1, double * RESTRICT _data_density, double * RESTRICT const _data_pdfs, double * RESTRICT const _data_pdfs_tmp, double * RESTRICT _data_velocity, int64_t const _stride_density_0, int64_t const _stride_density_1, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_tmp_0, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, double omega )
{
    const double xi_6 = 0.33333333333333331;
    const double xi_26 = -xi_6;
    const double xi_7 = 0.1111111111111111;
    const int64_t xi_8 = 2;
    const double xi_9 = 0.50000000000000000;
    const double xi_10 = 0.33333333333333331;
    const double xi_11 = 0.25000000000000000;

    const double xi_12 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 8*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 6*_stride_pdfs_2];
    const double xi_13 = xi_12 + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 4*_stride_pdfs_2];
    const double xi_14 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 7*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 5*_stride_pdfs_2];
    const double xi_15 = xi_14 + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 3*_stride_pdfs_2];
    const double xi_16 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 2*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2];
    const double xi_17 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 7*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 2*_stride_pdfs_2];
    const double vel0Term = xi_13;
    const double vel1Term = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 5*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2];
    const double partial_m_m1_e_0 = xi_15;
    const double partial_m_0_e_0 = xi_16 + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1];
    const double partial_m_1_e_0 = xi_13;
    const double xi_18 = partial_m_1_e_0 + partial_m_m1_e_0;
    const double partial_m_m1_e_1 = -_data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 7*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 5*_stride_pdfs_2];
    const double partial_m_0_e_1 = -_data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 2*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2];
    const double partial_m_1_e_1 = -_data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 8*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 6*_stride_pdfs_2];
    const double partial_m_m1_e_2 = xi_14;
    const double partial_m_0_e_2 = xi_16;
    const double partial_m_1_e_2 = xi_12;
    const double delta_rho = vel0Term + vel1Term + xi_17 + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + 3*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1];
    const double m_00 = partial_m_0_e_0 + xi_18 + 1.0;
    const double xi_3 = ((1.0) / (m_00));
    const double m_10 = partial_m_1_e_0 - partial_m_m1_e_0;
    const double u_0 = m_10*xi_3;
    const double xi_22 = u_0*u_0;
    const double xi_23 = m_00*xi_22;
    const double m_01 = partial_m_0_e_1 + partial_m_1_e_1 + partial_m_m1_e_1;
    const double u_1 = m_01*xi_3;
    const double xi_19 = u_1*u_1;
    const double xi_20 = m_00*xi_19;
    const double xi_21 = m_00*u_1;
    const double m_20 = xi_18 + xi_6;
    const double m_02 = partial_m_0_e_2 + partial_m_1_e_2 + partial_m_m1_e_2 + xi_6;
    const double m_11 = partial_m_1_e_1 - partial_m_m1_e_1;
    const double chimera_kappa_0_e_2 = m_02 - xi_20;
    const double chimera_kappa_1_e_1 = m_11 - u_0*xi_21;
    const double kappa_20 = m_20 - xi_23;
    const double C_4 = -chimera_kappa_0_e_2 + kappa_20;
    const double C_post_3 = -chimera_kappa_1_e_1*omega + chimera_kappa_1_e_1;
    const double xi_25 = C_post_3*((double)(xi_8));
    const double C_post_4 = -C_4*omega + C_4;
    const double xi_24 = C_post_4*xi_9;
    const double c_post_20 = m_00*xi_10 + xi_24;
    const double c_post_02 = m_00*xi_10 - xi_24;
    const double kappa_post_22 = c_post_02*c_post_20*xi_3 + xi_3*((double)(xi_8))*(C_post_3*C_post_3);
    const double chimera_m_post_0_e_1 = xi_21;
    const double chimera_m_post_0_e_2 = c_post_02 + xi_20;
    const double chimera_m_post_2_e_1 = c_post_20*u_1;
    const double chimera_m_post_1_e_2 = u_1*xi_25;
    const double chimera_m_post_2_e_2 = c_post_20*xi_19 + kappa_post_22;
    const double m_post_10 = m_00*u_0;
    const double m_post_20 = c_post_20 + xi_23;
    const double m_post_11 = C_post_3 + chimera_m_post_0_e_1*u_0;
    const double m_post_21 = chimera_m_post_0_e_1*xi_22 + chimera_m_post_2_e_1 + u_0*xi_25;
    const double xi_28 = -m_post_21;
    const double m_post_12 = chimera_m_post_0_e_2*u_0 + chimera_m_post_1_e_2;
    const double m_post_22 = chimera_m_post_0_e_2*xi_22 + chimera_m_post_1_e_2*u_0*((double)(xi_8)) + chimera_m_post_2_e_2;
    const double sub_k_to_f_0 = m_00 - 1.0;
    const double sub_k_to_f_3 = m_post_20 + xi_26;
    const double sub_k_to_f_4 = chimera_m_post_0_e_2 + xi_26;
    const double sub_k_to_f_8 = m_post_22 - xi_7;
    const double xi_27 = -sub_k_to_f_8;
    const double sub_k_to_f_10 = xi_9*(sub_k_to_f_4 + xi_27);
    const double sub_k_to_f_11 = xi_9*(chimera_m_post_0_e_1 + xi_28);
    const double sub_k_to_f_12 = xi_9*(sub_k_to_f_3 + xi_27);
    const double sub_k_to_f_13 = xi_9*(-m_post_10 + m_post_12);
    const double sub_k_to_f_14 = xi_11*(-m_post_11 - xi_27);
    const double sub_k_to_f_15 = xi_11*(-m_post_12 - xi_28);
    const double sub_k_to_f_16 = xi_11*(m_post_11 + sub_k_to_f_8);
    const double sub_k_to_f_17 = xi_11*(m_post_12 + m_post_21);
    const double momdensity_0 = vel0Term - xi_15;
    const double momdensity_1 = vel1Term - xi_17 - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + 8*_stride_pdfs_2] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + 6*_stride_pdfs_2];
    const double xi_0 = delta_rho + 1.0;
    const double xi_4 = ((1.0) / (xi_0));
    const double xi_1 = momdensity_0*xi_4;
    const double xi_2 = momdensity_1*xi_4;
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1] = sub_k_to_f_0 - sub_k_to_f_3 - sub_k_to_f_4 + sub_k_to_f_8;
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2] = sub_k_to_f_10 + sub_k_to_f_11;
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 2*_stride_pdfs_tmp_2] = sub_k_to_f_10 - sub_k_to_f_11;
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 3*_stride_pdfs_tmp_2] = sub_k_to_f_12 + sub_k_to_f_13;
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 4*_stride_pdfs_tmp_2] = sub_k_to_f_12 - sub_k_to_f_13;
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 5*_stride_pdfs_tmp_2] = sub_k_to_f_14 + sub_k_to_f_15;
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 6*_stride_pdfs_tmp_2] = sub_k_to_f_16 + sub_k_to_f_17;
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 7*_stride_pdfs_tmp_2] = sub_k_to_f_16 - sub_k_to_f_17;
    _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + 8*_stride_pdfs_tmp_2] = sub_k_to_f_14 - sub_k_to_f_15;
    _data_density[_stride_density_0*ctr_0 + _stride_density_1*ctr_1] = xi_0;
    _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1] = xi_1;
    _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2] = xi_2;
}

} // !namespace codegenKernels
} // !namespace lbm
} // !namespace walberla