import sympy as sp 
import pystencils as ps 

from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil

from lbmpy.creationfunctions import create_lb_update_rule

# Not sure if the following two are required but macroscopic_values_setter could offer a convenient way to initialize PDF fields (?)
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.boundaries import NoSlip, UBB, SimpleExtrapolationOutflow

from pystencils_walberla import CodeGeneration, generate_sweep, generate_pack_info_from_kernel

# No sure if this is required since for now creation of hardware-optimized boundary routines is beyond the scope of this project...
from lbmpy_walberla import generate_boundary

from pystencils.typing import TypedSymbol

# additional imports required for parabolic inflow UBB boundary condition:
from lbmpy_walberla.additional_data_handler import UBBAdditionalDataHandler

# Code generation:

with CodeGeneration() as ctx:
    # Parameters:
    data_type = 'float64' if ctx.double_accuracy else 'float32'
    stencil = LBStencil( Stencil.D2Q9 )
    omega = sp.Symbol( 'omega' )
    layout = 'fzyx'
    wall_velocity = sp.symbols( 'u_w, v_w')

    # PDF field: pdfs_tmp could remain unused because optimization is not intended to preserve portability of generated functions
    pdfs, pdfs_tmp = ps.fields( f'pdfs({stencil.Q}), pdfs_tmp({stencil.Q}): {data_type}[2D]', layout=layout )

    # Velocity and density output fields:
    velocity = ps.fields( f'velocity({stencil.D}): {data_type}[2D]', layout=layout )
    density = ps.fields( f'density(1): {data_type}[2D]', layout=layout )
    output = { 'velocity': velocity, 'density':density }

    # LBM optimisation: use common subexpression elimination (CSE) when possible to save computationa ressources (?)
    lbm_opt = LBMOptimisation( # cse_pdfs=True,
                               cse_global=True,
                               symbolic_field=pdfs,
                               symbolic_temporary_field=pdfs_tmp,
                               field_layout=layout )

    # Method setup
    lbm_config = LBMConfig( stencil=stencil,
                            method=Method.SRT,
                            relaxation_rate=omega,
                            compressible=True,
                            output=output )

    lbm_update_rule = create_lb_update_rule( lbm_config=lbm_config, lbm_optimisation=lbm_opt )
    lbm_method = lbm_update_rule.method

    # PDF initialization
    initial_rho = sp.Symbol( 'rho_0' )
    pdfs_setter = macroscopic_values_setter( lbm_method,
                                             initial_rho,
                                             velocity.center_vector,
                                             pdfs.center_vector )

    # Define hardware target:
    target = ps.Target.GPU if ctx.gpu else ps.Target.CPU 

    # non-constant inflow UBB:
    inflowFuncUBB = UBB( lambda *args: None, dim=stencil.D, data_type=data_type )
    inflowFuncUBBDataHandler = UBBAdditionalDataHandler( stencil, inflowFuncUBB )

    # Code generation:
    # LBM Sweep
    generate_sweep( ctx, 'D2Q9SRTSweep', lbm_update_rule, field_swaps=[( pdfs, pdfs_tmp )], target=target )

    # Pack info
    generate_pack_info_from_kernel( ctx, 'D2Q9SRTPackInfo', lbm_update_rule, target=target )

    # Macroscopic values setter (PDF field initializer)
    generate_sweep( ctx, 'D2Q9SRTInitialPdfsSetter', pdfs_setter, target=target )

    # Boundaries:
    generate_boundary( ctx, "D2Q9SRTNoSlip", NoSlip(), lbm_method, target=target )

    generate_boundary( ctx, "D2Q9SRTUBB", UBB( velocity=wall_velocity, dim=stencil.D, data_type=data_type ), lbm_method, target=target )

    generate_boundary( ctx, "D2Q9SRTSEO_E", SimpleExtrapolationOutflow( normal_direction=stencil[4], stencil=stencil ), lbm_method, target=target )

    # non-constant inflow boundary:

    generate_boundary( ctx, "D2Q9SRTInflowFunc", inflowFuncUBB, lbm_method, target=target, additional_data_handler=inflowFuncUBBDataHandler )
