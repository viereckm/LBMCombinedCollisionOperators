walberla_generate_target_from_python(
    NAME D2Q9CumulantSweep
    FILE D2Q9CumulantSweep.py
    OUT_FILES   D2Q9CumulantSweep.cpp               D2Q9CumulantSweep.h
                D2Q9CumulantPackInfo.cpp            D2Q9CumulantPackInfo.h
                D2Q9CumulantInitialPdfsSetter.cpp   D2Q9CumulantInitialPdfsSetter.h
                D2Q9CumulantNoSlip.cpp              D2Q9CumulantNoSlip.h
                D2Q9CumulantUBB.cpp                 D2Q9CumulantUBB.h
                D2Q9CumulantSEO_E.cpp               D2Q9CumulantSEO_E.h
                D2Q9CumulantInflowFunc.cpp          D2Q9CumulantInflowFunc.h
)
                
walberla_generate_target_from_python(
    NAME D2Q9SRTSweep
    FILE D2Q9SRTSweep.py
    OUT_FILES   D2Q9SRTSweep.cpp                D2Q9SRTSweep.h
                D2Q9SRTPackInfo.cpp             D2Q9SRTPackInfo.h 
                D2Q9SRTInitialPdfsSetter.cpp    D2Q9SRTInitialPdfsSetter.h
                D2Q9SRTNoSlip.cpp               D2Q9SRTNoSlip.h
                D2Q9SRTUBB.cpp                  D2Q9SRTUBB.h
                D2Q9SRTSEO_E.cpp                D2Q9SRTSEO_E.h
                D2Q9SRTInflowFunc.cpp            D2Q9SRTInflowFunc.h
)

walberla_generate_target_from_python(
    NAME D2Q9TRTSweep
    FILE D2Q9TRTSweep.py
    OUT_FILES   D2Q9TRTSweep.cpp                D2Q9TRTSweep.h
                D2Q9TRTPackInfo.cpp             D2Q9TRTPackInfo.h
                D2Q9TRTInitialPdfsSetter.cpp    D2Q9TRTInitialPdfsSetter.h
                D2Q9TRTNoSlip.cpp               D2Q9TRTNoSlip.h
                D2Q9TRTUBB.cpp                  D2Q9TRTUBB.h
                D2Q9TRTSEO_E.cpp                D2Q9TRTSEO_E.h
                D2Q9TRTInflowFunc.cpp           D2Q9TRTInflowFunc.h
)