import numpy as np 

def loadArrays( Re, omega, uMax, N=100 ):
    path = f'./results/BreuerGuoRe{int(Re)}omega{int(10*omega)}uMax{int(100*uMax)}/'
    C_D_file = path + 'C_D.csv'
    C_L_file = path + 'C_L.csv'
    time_file = path + 'time.csv'

    C_D_arr = np.genfromtxt( C_D_file, delimiter=',' )
    C_L_arr = np.genfromtxt( C_L_file, delimiter=',' )
    time_arr = np.genfromtxt( time_file, delimiter=',' )

    return C_D_arr[-N:], C_L_arr[-N:], time_arr[-N:]

def loadArraysFixedGrid( Re, Nx, Ny, N=100 ):
    path = f'./results/BreuerGuoNx{int(Nx)}Ny{int(Ny)}Re{int(Re)}uMax10/'
    C_D_file = path + 'C_D.csv'
    C_L_file = path + 'C_L.csv'
    time_file = path + 'time.csv'

    C_D_arr = np.genfromtxt( C_D_file, delimiter=',' )
    C_L_arr = np.genfromtxt( C_L_file, delimiter=',' )
    time_arr = np.genfromtxt( time_file, delimiter=',' )

    return C_D_arr[-N:], C_L_arr[-N:], time_arr[-N:]
    

def getOnePeriod( C_L_arr ):
    begin = 0
    end = 0
    monotonicity = np.sign( C_L_arr[1:] - C_L_arr[:-1])

    for i in range( 0, len( monotonicity ) ):
        if monotonicity[ i+1 ] < monotonicity[ i ]:
            begin = i 
            break

    for i in range( begin+1, len( monotonicity ) ):
        if monotonicity[ i+1 ] < monotonicity[ i ]:
            end = i
            break

    return begin, end

def meanAndOscillation( C_D_one_period, C_L_one_period ):
    C_D_mean = np.mean( C_D_one_period )
    C_L_mean = np.mean( C_L_one_period )

    C_D_Oscillation = np.amax( C_D_one_period ) - np.amin( C_D_one_period )
    C_L_Oscillation = np.amax( C_L_one_period ) - np.amin( C_L_one_period )

    return C_D_mean, C_L_mean, C_D_Oscillation, C_L_Oscillation

def nuFromOmega( omega ):
    return ( 1. / omega - 0.5 ) / 3.

def ObstacleDimensionFromRe( Re, uMax, nu ):
    D = int( Re * nu / uMax + 1 )
    
    if( D % 2 == 1 ):
        D = D + 1

    return D 

def velocityCorrection( Re, nu, D ):
    return Re*nu/D

def computeStouhal( time_arr, begin, end, D, uMax ):
    nSteps = time_arr[ end ] - time_arr[ begin ]
    frequency = 1 / nSteps
    return frequency * D / uMax





