#include <iostream>
#include <cassert>
#define ASSERT(expr) assert(expr)

#include "Common/binary_tree.h"


int main(int argc, char **argv) {
    COMM::Binary bin;
    bin.period =  0.00002737925;
    bin.ecc = 0.0;
    bin.m1 = 0.4;
    bin.m2 = 110.0;

    Float G = 0.00449830997959438;

    bin.calcSemiFromPeriod(G);
    std::cout<<"semi:"<<bin.semi<<" period:"<<bin.period<<std::endl;

    Float mtot=bin.m1+bin.m2;
    bin.semi = bin.periodToSemi(bin.period, mtot, G);
    std::cout<<"semi:"<<bin.semi<<" period:"<<bin.period<<std::endl;

    return 0;
}
