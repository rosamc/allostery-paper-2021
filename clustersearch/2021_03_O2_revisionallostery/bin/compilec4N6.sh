module load gcc/6.2.0
module load python/3.6.0
source /home/rm335/myenvpy36/bin/activate
c++ -O2 -DNDEBUG -Wall -shared -std=c++11  -fPIC -I /n/groups/gunawardena/shared/boost-1.71.0/build/include/ -I /home/rm335/repos/polynomials/include/polynomial/ -I /home/rm335/libs/eigen-eigen-323c052e1731 -I /home/rm335/repos/sharedposstpNov19/GeneRegulatoryFunctions/utilsGRF  `python3 -m pybind11 --includes` ./CG_c4_N6_samesitesFalse.cpp -lmpfr -lmpc -o ./CG_c4_N6_samesitesFalse`python3-config --extension-suffix`
