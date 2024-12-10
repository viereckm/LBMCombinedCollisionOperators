target_sources( lbm
    PRIVATE
    CollisionKernels.hpp 
    CustomCollisionModels.cpp
    MixedCollisionSweep.hpp
    CodegenCollisionKernels.cpp
    D2Q9CodegenPackInfo.cpp
    CodegenPdfsSetterKernels.cpp
    D2Q9CodegenInitialPdfsSetter.cpp
    D2Q9CodegenNoSlip.cpp
    D2Q9CodegenUBB.cpp
    D2Q9CodegenSEO_E.cpp
    MixedCodegenCollisionSweep.cpp
    D2Q9InflowFunc.cpp
    InflowVelocitySetter.cpp
    InflowProfiles.hpp
    my_utility.hpp
    MomentumExchangeCalculation.cpp
    LiftAndDrag.cpp
    PressureDrop.cpp
    WriteFields.hpp
    vorticity.cpp
)

