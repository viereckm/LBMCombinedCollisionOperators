#include "CodegenPdfsSetterKernels.hpp"

namespace walberla {
namespace lbm {
namespace codegenKernels {

void D2Q9SRTPdfsSetter( int64_t const ctr_0, int64_t const ctr_1, double * RESTRICT _data_pdfs, double * RESTRICT const _data_velocity, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, double rho_0 )
{
    const double rho = rho_0;
    const double delta_rho = rho - 1.0;

    const double u_0 = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1];
    const double u_1 = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2];
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1] = delta_rho*0.44444444444444442 + rho*-0.66666666666666663*(u_0*u_0) + rho*-0.66666666666666663*(u_1*u_1);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2] = delta_rho*0.1111111111111111 + rho*u_1*0.33333333333333331 + rho*-0.16666666666666666*(u_0*u_0) + rho*0.33333333333333331*(u_1*u_1);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 2*_stride_pdfs_2] = delta_rho*0.1111111111111111 + rho*u_1*-0.33333333333333331 + rho*-0.16666666666666666*(u_0*u_0) + rho*0.33333333333333331*(u_1*u_1);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 3*_stride_pdfs_2] = delta_rho*0.1111111111111111 + rho*u_0*-0.33333333333333331 + rho*-0.16666666666666666*(u_1*u_1) + rho*0.33333333333333331*(u_0*u_0);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 4*_stride_pdfs_2] = delta_rho*0.1111111111111111 + rho*u_0*0.33333333333333331 + rho*-0.16666666666666666*(u_1*u_1) + rho*0.33333333333333331*(u_0*u_0);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 5*_stride_pdfs_2] = delta_rho*0.027777777777777776 + rho*u_0*u_1*-0.25 + rho*u_0*-0.083333333333333329 + rho*u_1*0.083333333333333329 + rho*0.083333333333333329*(u_0*u_0) + rho*0.083333333333333329*(u_1*u_1);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 6*_stride_pdfs_2] = delta_rho*0.027777777777777776 + rho*u_0*u_1*0.25 + rho*u_0*0.083333333333333329 + rho*u_1*0.083333333333333329 + rho*0.083333333333333329*(u_0*u_0) + rho*0.083333333333333329*(u_1*u_1);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 7*_stride_pdfs_2] = delta_rho*0.027777777777777776 + rho*u_0*u_1*0.25 + rho*u_0*-0.083333333333333329 + rho*u_1*-0.083333333333333329 + rho*0.083333333333333329*(u_0*u_0) + rho*0.083333333333333329*(u_1*u_1);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 8*_stride_pdfs_2] = delta_rho*0.027777777777777776 + rho*u_0*u_1*-0.25 + rho*u_0*0.083333333333333329 + rho*u_1*-0.083333333333333329 + rho*0.083333333333333329*(u_0*u_0) + rho*0.083333333333333329*(u_1*u_1);
}

void D2Q9TRTPdfsSetter( int64_t const ctr_0, int64_t const ctr_1, double * RESTRICT _data_pdfs, double * RESTRICT const _data_velocity, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, double rho_0 )
{
    const double rho = rho_0;
    const double delta_rho = rho - 1.0;

    const double u_0 = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1];
    const double u_1 = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2];
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1] = delta_rho*0.44444444444444442 + rho*-0.66666666666666663*(u_0*u_0) + rho*-0.66666666666666663*(u_1*u_1);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2] = delta_rho*0.1111111111111111 + rho*u_1*0.33333333333333331 + rho*-0.16666666666666666*(u_0*u_0) + rho*0.33333333333333331*(u_1*u_1);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 2*_stride_pdfs_2] = delta_rho*0.1111111111111111 + rho*u_1*-0.33333333333333331 + rho*-0.16666666666666666*(u_0*u_0) + rho*0.33333333333333331*(u_1*u_1);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 3*_stride_pdfs_2] = delta_rho*0.1111111111111111 + rho*u_0*-0.33333333333333331 + rho*-0.16666666666666666*(u_1*u_1) + rho*0.33333333333333331*(u_0*u_0);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 4*_stride_pdfs_2] = delta_rho*0.1111111111111111 + rho*u_0*0.33333333333333331 + rho*-0.16666666666666666*(u_1*u_1) + rho*0.33333333333333331*(u_0*u_0);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 5*_stride_pdfs_2] = delta_rho*0.027777777777777776 + rho*u_0*u_1*-0.25 + rho*u_0*-0.083333333333333329 + rho*u_1*0.083333333333333329 + rho*0.083333333333333329*(u_0*u_0) + rho*0.083333333333333329*(u_1*u_1);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 6*_stride_pdfs_2] = delta_rho*0.027777777777777776 + rho*u_0*u_1*0.25 + rho*u_0*0.083333333333333329 + rho*u_1*0.083333333333333329 + rho*0.083333333333333329*(u_0*u_0) + rho*0.083333333333333329*(u_1*u_1);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 7*_stride_pdfs_2] = delta_rho*0.027777777777777776 + rho*u_0*u_1*0.25 + rho*u_0*-0.083333333333333329 + rho*u_1*-0.083333333333333329 + rho*0.083333333333333329*(u_0*u_0) + rho*0.083333333333333329*(u_1*u_1);
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 8*_stride_pdfs_2] = delta_rho*0.027777777777777776 + rho*u_0*u_1*-0.25 + rho*u_0*0.083333333333333329 + rho*u_1*-0.083333333333333329 + rho*0.083333333333333329*(u_0*u_0) + rho*0.083333333333333329*(u_1*u_1);
}

void D2Q9CumulantPdfsSetter( int64_t const ctr_0, int64_t const ctr_1, double * RESTRICT _data_pdfs, double * RESTRICT const _data_velocity, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, double rho_0 )
{
    const double rho = rho_0;

    const double u_0 = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1];
    const double u_1 = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2];
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1] = rho*-0.66666666666666663*(u_0*u_0) + rho*-0.66666666666666663*(u_1*u_1) + rho*0.44444444444444442 + rho*(u_0*u_0)*(u_1*u_1) - 0.44444444444444442;
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2] = rho*u_1*-0.5*(u_0*u_0) + rho*u_1*0.33333333333333331 + rho*-0.16666666666666666*(u_0*u_0) + rho*-0.5*(u_0*u_0)*(u_1*u_1) + rho*0.1111111111111111 + rho*0.33333333333333331*(u_1*u_1) - 0.1111111111111111;
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 2*_stride_pdfs_2] = rho*u_1*-0.33333333333333331 + rho*u_1*0.5*(u_0*u_0) + rho*-0.16666666666666666*(u_0*u_0) + rho*-0.5*(u_0*u_0)*(u_1*u_1) + rho*0.1111111111111111 + rho*0.33333333333333331*(u_1*u_1) - 0.1111111111111111;
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 3*_stride_pdfs_2] = rho*u_0*-0.33333333333333331 + rho*u_0*0.5*(u_1*u_1) + rho*-0.16666666666666666*(u_1*u_1) + rho*-0.5*(u_0*u_0)*(u_1*u_1) + rho*0.1111111111111111 + rho*0.33333333333333331*(u_0*u_0) - 0.1111111111111111;
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 4*_stride_pdfs_2] = rho*u_0*-0.5*(u_1*u_1) + rho*u_0*0.33333333333333331 + rho*-0.16666666666666666*(u_1*u_1) + rho*-0.5*(u_0*u_0)*(u_1*u_1) + rho*0.1111111111111111 + rho*0.33333333333333331*(u_0*u_0) - 0.1111111111111111;
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 5*_stride_pdfs_2] = rho*u_0*u_1*-0.25 + rho*u_0*-0.083333333333333329 + rho*u_0*-0.25*(u_1*u_1) + rho*u_1*0.083333333333333329 + rho*u_1*0.25*(u_0*u_0) + rho*0.027777777777777776 + rho*0.083333333333333329*(u_0*u_0) + rho*0.083333333333333329*(u_1*u_1) + rho*0.25*(u_0*u_0)*(u_1*u_1) - 0.027777777777777776;
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 6*_stride_pdfs_2] = rho*u_0*u_1*0.25 + rho*u_0*0.083333333333333329 + rho*u_0*0.25*(u_1*u_1) + rho*u_1*0.083333333333333329 + rho*u_1*0.25*(u_0*u_0) + rho*0.027777777777777776 + rho*0.083333333333333329*(u_0*u_0) + rho*0.083333333333333329*(u_1*u_1) + rho*0.25*(u_0*u_0)*(u_1*u_1) - 0.027777777777777776;
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 7*_stride_pdfs_2] = rho*u_0*u_1*0.25 + rho*u_0*-0.083333333333333329 + rho*u_0*-0.25*(u_1*u_1) + rho*u_1*-0.083333333333333329 + rho*u_1*-0.25*(u_0*u_0) + rho*0.027777777777777776 + rho*0.083333333333333329*(u_0*u_0) + rho*0.083333333333333329*(u_1*u_1) + rho*0.25*(u_0*u_0)*(u_1*u_1) - 0.027777777777777776;
    _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + 8*_stride_pdfs_2] = rho*u_0*u_1*-0.25 + rho*u_0*0.083333333333333329 + rho*u_0*0.25*(u_1*u_1) + rho*u_1*-0.083333333333333329 + rho*u_1*-0.25*(u_0*u_0) + rho*0.027777777777777776 + rho*0.083333333333333329*(u_0*u_0) + rho*0.083333333333333329*(u_1*u_1) + rho*0.25*(u_0*u_0)*(u_1*u_1) - 0.027777777777777776;
}

} // !namespace codegenKernels
} // !namespace lbm
} // !namespace walberla