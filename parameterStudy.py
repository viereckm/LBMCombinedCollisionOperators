import numpy as np
import matplotlib.pyplot as plt 

import PostProcessing as pp 

Res     = np.array( [ 10, 30, 60, 100, 133, 150 ] )
omegas  = np.array( [ 1.2, 1.4, 1.6 ] )
uMaxs   = np.array( [ 0.1, 0.15, 0.2 ] )

nRes    = len( Res )
nOmegas = len( omegas )
nuMaxs  = len( uMaxs )

C_D     = np.zeros( ( nRes, nOmegas, nuMaxs ) )
C_L     = np.zeros( ( nRes, nOmegas, nuMaxs ) )
C_D_osc = np.zeros( ( nRes, nOmegas, nuMaxs ) )
C_L_osc = np.zeros( ( nRes, nOmegas, nuMaxs ) )
Stouhal = np.zeros( ( nRes, nOmegas, nuMaxs ) )

steady  = Res <= 60
unsteady = Res >= 60

begin = 0
end = 0

for i in range( 0, nRes ):
    for j in range( 0, nOmegas ):
        for k in range( 0, nuMaxs ):

            Re = Res[ i ]
            omega = omegas[ j ]
            uMax = uMaxs[ k ]

            C_D_arr, C_L_arr, time_arr = pp.loadArrays( Re=Re, omega=omega, uMax=uMax, N=200 )
            
            if Re > 50:
                begin, end = pp.getOnePeriod( C_L_arr )
                # print( f"begin = {begin}\nend = {end}\n")
                C_D_mean, C_L_mean, C_D_Oscillation, C_L_Oscillation = pp.meanAndOscillation( C_D_arr[ begin:end ], C_L_arr[ begin:end ] )
            else:
                C_D_mean = C_D_arr[ 0 ]
                C_L_mean = C_L_arr[ 0 ]
                C_D_Oscillation = 0
                C_L_Oscillation = 0

            correction = 4. / 9.

            C_D[ i, j, k ]      = correction * C_D_mean
            C_L[ i, j, k ]      = correction * C_L_mean
            C_D_osc[ i, j, k ]  = correction * C_D_Oscillation
            C_L_osc[ i, j, k ]  = correction * C_L_Oscillation

            if Re > 50:
                nu = pp.nuFromOmega( omega )
                D = pp.ObstacleDimensionFromRe( Re, uMax, nu )
                uMax_corr = pp.velocityCorrection( Re, nu, D )
                Stouhal[ i, j, k ] = pp.computeStouhal( time_arr, begin, end, D, uMax_corr )



for j in range( 0, nOmegas ):

    # plot lift coefficient amplitude
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( '|max(C_L) - min(C_L)|')
    plt.ylim( ( -0.1, 1.1) )
    plt.title( f'omega = {omegas[j]:.2}')
    plt.grid()

    for k in range( 0, nuMaxs ):
        label = f'uMax = {uMaxs[k]:.2}'
        plt.plot( Res[ unsteady ], C_L_osc[ unsteady, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/CL_oscillation_omega{int(10*omegas[j])}.png', format='png')
    plt.close( fig )


    # plot drag coefficient amplitude
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( '|max(C_D) - min(C_D)|')
    plt.ylim( ( 0, 0.1) )
    plt.title( f'omega = {omegas[j]:.2}')
    plt.grid()

    for k in range( 0, nuMaxs ):
        label = f'uMax = {uMaxs[k]:.2}'
        plt.plot( Res[ unsteady ], C_D_osc[ unsteady, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/CD_oscillation_omega{int(10*omegas[j])}.png', format='png')
    plt.close( fig )

    
    # plot unsteady lift coefficient
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( 'C_L_mean')
    plt.ylim( ( -0.1, 0.1) )
    plt.title( f'omega = {omegas[j]:.2}')
    plt.grid()

    for k in range( 0, nuMaxs ):
        label = f'uMax = {uMaxs[k]:.2}'
        plt.plot( Res[ unsteady ], C_L[ unsteady, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/CL_unsteady_omega{int(10*omegas[j])}.png', format='png')
    plt.close( fig )
    
    # plot unsteady drag coefficient
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( 'C_D_mean')
    plt.ylim( ( 1.1, 1.45) )
    plt.title( f'omega = {omegas[j]:.2}')
    plt.grid()

    for k in range( 0, nuMaxs ):
        label = f'uMax = {uMaxs[k]:.2}'
        plt.plot( Res[ unsteady ], C_D[ unsteady, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/CD_unsteady_omega{int(10*omegas[j])}.png', format='png')
    plt.close( fig )

    # plot steady lift coefficient
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( 'C_L')
    plt.ylim( ( -0.1, 0.1) )
    plt.title( f'omega = {omegas[j]:.2}')
    plt.grid()

    for k in range( 0, nuMaxs ):
        label = f'uMax = {uMaxs[k]:.2}'
        plt.plot( Res[ steady ], C_L[ steady, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/CL_steady_omega{int(10*omegas[j])}.png', format='png')
    plt.close( fig )
    
    # plot steady drag coefficient
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( 'C_D')
    plt.ylim( ( 1.2, 2.8 ) )
    plt.title( f'omega = {omegas[j]:.2}')
    plt.grid()

    for k in range( 0, nuMaxs ):
        label = f'uMax = {uMaxs[k]:.2}'
        plt.plot( Res[ steady ], C_D[ steady, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/CD_steady_omega{int(10*omegas[j])}.png', format='png')
    plt.close( fig )

    # plot drag coefficient
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( 'C_D')
    plt.ylim( ( 1.0, 3.0 ) )
    plt.title( f'omega = {omegas[j]:.2}')
    plt.grid()

    for k in range( 0, nuMaxs ):
        label = f'uMax = {uMaxs[k]:.2}'
        plt.plot( Res, C_D[ :, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/CD_omega{int(10*omegas[j])}.png', format='png')
    plt.close( fig )

    # plot Stouhal number
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( 'St')
    plt.ylim( ( 0.1, 0.16 ) )
    plt.title( f'omega = {omegas[j]:.2}')
    plt.grid()

    for k in range( 0, nuMaxs ):
        label = f'uMax = {uMaxs[k]:.2}'
        plt.plot( Res[ unsteady ], Stouhal[ unsteady, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/St_omega{int(10*omegas[j])}.png', format='png')
    plt.close( fig )


for k in range( 0, nuMaxs ):

    # plot lift coefficient amplitude
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( '|max(C_L) - min(C_L)|')
    plt.ylim( ( -0.1, 1.1) )
    plt.title( f'uMax = {uMaxs[k]:.2}')
    plt.grid()

    for j in range( 0, nOmegas ):
        label = f'omega = {omegas[j]:.2}'
        plt.plot( Res[ unsteady ], C_L_osc[ unsteady, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/CL_oscillation_uMAx{int(100*uMaxs[k])}.png', format='png')
    plt.close( fig )


    # plot drag coefficient amplitude
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( '|max(C_D) - min(C_D)|')
    plt.ylim( ( 0, 0.1) )
    plt.title( f'uMax = {uMaxs[k]:.2}')
    plt.grid()

    for j in range( 0, nOmegas ):
        label = f'omega = {omegas[j]:.2}'
        plt.plot( Res[ unsteady ], C_D_osc[ unsteady, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/CD_oscillation_uMax{int(100*uMaxs[k])}.png', format='png')
    plt.close( fig )

    
    # plot unsteady lift coefficient
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( 'C_L_mean')
    plt.ylim( ( -0.1, 0.1) )
    plt.title( f'uMax = {uMaxs[k]:.2}')
    plt.grid()

    for j in range( 0, nOmegas ):
        label = f'omega = {omegas[j]:.2}'
        plt.plot( Res[ unsteady ], C_L[ unsteady, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/CL_unsteady_uMax{int(100*uMaxs[k])}.png', format='png')
    plt.close( fig )
    
    # plot unsteady drag coefficient
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( 'C_D_mean')
    plt.ylim( ( 1.1, 1.4) )
    plt.title( f'uMax = {uMaxs[k]:.2}')
    plt.grid()

    for j in range( 0, nOmegas ):
        label = f'omega = {omegas[j]:.2}'
        plt.plot( Res[ unsteady ], C_D[ unsteady, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/CD_unsteady_uMax{int(100*uMaxs[k])}.png', format='png')
    plt.close( fig )

    # plot steady lift coefficient
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( 'C_L')
    plt.ylim( ( -0.1, 0.1) )
    plt.title( f'uMax = {uMaxs[k]:.2}')
    plt.grid()

    for j in range( 0, nOmegas ):
        label = f'omega = {omegas[j]:.2}'
        plt.plot( Res[ steady ], C_L[ steady, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/CL_steady_uMax{int(100*uMaxs[k])}.png', format='png')
    plt.close( fig )
    
    # plot steady drag coefficient
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( 'C_D')
    plt.ylim( ( 1.2, 2.8 ) )
    plt.title( f'uMax = {uMaxs[k]:.2}')
    plt.grid()

    for j in range( 0, nOmegas ):
        label = f'omega = {omegas[j]:.2}'
        plt.plot( Res[ steady ], C_D[ steady, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/CD_steady_uMax{int(100*uMaxs[k])}.png', format='png')
    plt.close( fig )

    # plot drag coefficient
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( 'C_D')
    plt.ylim( ( 1.0, 3.0 ) )
    plt.title( f'uMax = {uMaxs[k]:.2}')
    plt.grid()

    for j in range( 0, nOmegas ):
        label = f'omega = {omegas[j]:.2}'
        plt.plot( Res, C_D[ :, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/CD_uMax{int(100*uMaxs[k])}.png', format='png')
    plt.close( fig )

    # plot Stouhal number
    fig = plt.figure()
    plt.xlabel( 'Re' )
    plt.ylabel( 'St')
    plt.ylim( ( 0.1, 0.16 ) )
    plt.title( f'uMax = {uMaxs[k]:.2}')
    plt.grid()

    for j in range( 0, nOmegas ):
        label = f'omega = {omegas[j]:.2}'
        plt.plot( Res[ unsteady ], Stouhal[ unsteady, j, k ], 'o-', label=label )

    plt.legend()
    plt.savefig( f'../figures/St_uMax{int(100*uMaxs[k])}.png', format='png')
    plt.close( fig )
