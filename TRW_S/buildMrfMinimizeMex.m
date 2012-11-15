%mex command to build mrfMinimizeMex
mex src/mrfMinimizeMex.cpp src/ordering.cpp src/MRFEnergy.cpp src/treeProbabilities.cpp src/minimize.cpp -output mrfMinimizeMex -largeArrayDims
