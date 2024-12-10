import numpy as np
import matplotlib.pyplot as plt

import PostProcessing as pp

Res     = np.array( [ 10, 30, 60, 100, 133, 150 ] )
Nxs     = np.array( [ 500, 1000 ] )
Nys     = np.array( [ 80, 160 ] )
Ds      = np.array( [ 10, 20 ] )

uMax    = 0.1

nNxs    = len( Nxs )
nRes = len( Res )

C_D     = np.zeros( ( nRes, nNxs ) )
C_L     = np.zeros( ( nRes, nNxs ) )
C_D_osc = np.zeros( ( nRes, nNxs ) )
C_L_osc = np.zeros( ( nRes, nNxs ) )
Stouhal = np.zeros( ( nRes, nNxs ) )

steady      = Res <= 60
unsteady    = Res >= 60

begin   = 0
end     = 0

for i in range( 0, nRes ):
    for j in range( 0, nNxs ):
        
        Re = Res[ i ]
        Nx = Nxs[ j ]
        Ny = Nys[ j ]
        D  = Ds[ j ]

        C_D_arr, C_L_arr, time_arr = pp.loadArraysFixedGrid( Re=Re, Nx=Nx, Ny=Ny, N=200 )

        if Re > 50:
            begin, end = pp.getOnePeriod( C_L_arr )
            C_D_mean, C_L_mean, C_D_Oscillation, C_L_Oscillation = pp.meanAndOscillation( C_D_arr[ begin:end ], C_L_arr[ begin:end ] )
        else:
            C_D_mean = C_D_arr[ 0 ]
            C_L_mean = C_L_arr[ 0 ]
            C_D_Oscillation = 0
            C_L_Oscillation = 0
    
        correction = 4. / 9.

        C_D[ i, j ]         = correction * C_D_mean
        C_L[ i, j ]         = correction * C_L_mean
        C_D_osc[ i, j ]     = correction * C_D_Oscillation
        C_L_osc[ i, j ]     = correction * C_L_Oscillation

        if Re > 50:
            Stouhal[ i, j ]     = pp.computeStouhal( time_arr, begin, end, D, uMax )

# plot lift coeff amplitude
fig = plt.figure()
plt.xlabel( 'Re' )
plt.ylabel( '|max(C_L) - min(C_L)|' )
plt.ylim( ( -0.1, 1.1 ) )
plt.title( f'uMax = {uMax:.2}' )
plt.grid()

for j in range( 0, nNxs ):
    label = f'Nx = {int(Nxs[j])}; Ny = {int(Nys[j])}'
    plt.plot( Res[ unsteady ], C_L_osc[ unsteady, j ], 'o-', label=label )

plt.legend()
plt.savefig( f'../figures/CL_oscillation_Grid.png', format='png' )
plt.close( fig ) 

# plot drag coeff amplitude
fig = plt.figure()
plt.xlabel( 'Re' )
plt.ylabel( '|max(C_D) - min(C_D)|' )
plt.ylim( ( 0, .1 ) )
plt.title( f'uMax = {uMax:.2}' )
plt.grid()

for j in range( 0, nNxs ):
    label = f'Nx = {int(Nxs[j])}; Ny = {int(Nys[j])}'
    plt.plot( Res[ unsteady ], C_D_osc[ unsteady, j ], 'o-',  label=label )

plt.legend()
plt.savefig( f'../figures/CD_oscillation_Grid.png', format='png' )
plt.close( fig ) 

# plot unsteady lift coeff
fig = plt.figure()
plt.xlabel( 'Re' )
plt.ylabel( 'C_L_mean' )
plt.ylim( ( -0.1, 0.1 ) )
plt.title( f'uMax = {uMax:.2}' )
plt.grid()

for j in range( 0, nNxs ):
    label = f'Nx = {int(Nxs[j])}; Ny = {int(Nys[j])}'
    plt.plot( Res[ unsteady ], C_L[ unsteady, j ], 'o-', label=label )

plt.legend()
plt.savefig( f'../figures/CL_unsteady_Grid.png', format='png' )
plt.close( fig ) 

# plot unsteady drag coeff
fig = plt.figure()
plt.xlabel( 'Re' )
plt.ylabel( 'C_D_mean' )
plt.ylim( ( 1.1, 1.45 ) )
plt.title( f'uMax = {uMax:.2}' )
plt.grid()

for j in range( 0, nNxs ):
    label = f'Nx = {int(Nxs[j])}; Ny = {int(Nys[j])}'
    plt.plot( Res[ unsteady ], C_D[ unsteady, j ], 'o-', label=label )

plt.legend()
plt.savefig( f'../figures/CD_unsteady_Grid.png', format='png' )
plt.close( fig ) 

# plot steady lift coeff
fig = plt.figure()
plt.xlabel( 'Re' )
plt.ylabel( 'C_L' )
plt.ylim( ( -0.1, 0.1 ) )
plt.title( f'uMax = {uMax:.2}' )
plt.grid()

for j in range( 0, nNxs ):
    label = f'Nx = {int(Nxs[j])}; Ny = {int(Nys[j])}'
    plt.plot( Res[ steady ], C_L[ steady, j ], 'o-', label=label )

plt.legend()
plt.savefig( f'../figures/CL_steady_Grid.png', format='png' )
plt.close( fig ) 

# plot steady drag coeff
fig = plt.figure()
plt.xlabel( 'Re' )
plt.ylabel( 'C_D' )
plt.ylim( ( 1.2, 2.8 ) )
plt.title( f'uMax = {uMax:.2}' )
plt.grid()

for j in range( 0, nNxs ):
    label = f'Nx = {int(Nxs[j])}; Ny = {int(Nys[j])}'
    plt.plot( Res[ steady ], C_D[ steady, j ], 'o-', label=label )

plt.legend()
plt.savefig( f'../figures/CD_steady_Grid.png', format='png' )
plt.close( fig ) 

# plot drag coeff
fig = plt.figure()
plt.xlabel( 'Re' )
plt.ylabel( 'C_D' )
plt.ylim( ( 1.0, 3.0 ) )
plt.title( f'uMax = {uMax:.2}' )
plt.grid()

for j in range( 0, nNxs ):
    label = f'Nx = {int(Nxs[j])}; Ny = {int(Nys[j])}'
    plt.plot( Res, C_D[ :, j ], 'o-', label=label )

plt.legend()
plt.savefig( f'../figures/CD_Grid.png', format='png' )
plt.close( fig ) 

# plot Stouhal number
fig = plt.figure()
plt.xlabel( 'Re' )
plt.ylabel( 'St' )
plt.ylim( ( 0.1, .16 ) )
plt.title( f'uMax = {uMax:.2}' )
plt.grid()

for j in range( 0, nNxs ):
    label = f'Nx = {int(Nxs[j])}; Ny = {int(Nys[j])}'
    plt.plot( Res[ unsteady ], Stouhal[ unsteady, j ], 'o-', label=label )

plt.legend()
plt.savefig( f'../figures/St_Grid.png', format='png' )
plt.close( fig ) 