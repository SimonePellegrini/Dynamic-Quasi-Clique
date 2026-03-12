This repository contains the codes for the paper "Detecting Large Quasi-cliques on Dynamic Networks". 

- CreditsAlgorithm: contains the implementation of the fully dynamic credits algorithm
- CreditsAlgorithmOnlyIns: contains the implementation of the incremental version of the credits algorithm
- DynamicMinHash: is an implementation of the dynamic minhash data structure cited in the paper
- DynamincNBSim: contains the implementation of the fully dynamic baseline
- TabulationHash: contains the implementation of the tabulation hashing used to compute the minhash signature
- GenerateSequences: contains the code used to generate the synthetic sequences cited in the paper
- main.c++ contains the code used to generate the experiments

To compile the code you need gcc 10 or above and to use the O3 optimizer.

g++ -o dyn -O3 main.c++

