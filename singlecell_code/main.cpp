/*
 *  main_cell_dts.h
 *
 *  Copyright 2010 Colleen Clancy Lab. All rights reserved.
 *
 */

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <time.h>
#include <omp.h>


using namespace std;


typedef struct cell {
    double y0n[112];
    double highestV, lowestV, Apd90, DI, t1, t2;
    double highestdvdt, V90, t3, t4, s2APD, V80, Apd80;
    double dvdt, v, v_new, v_old, dvdt2;
    double currents[2];
    double I_Total, I_inj, I_start, I_end;
    double ui1[4], uj1[4], ui2[4], uj2[4], dun[4], tu[4];
    
} Cell;


const int tw = 1;
const int tl = 1;



typedef struct SimState1 {
    double t, tt;
    int tstep, counter, beat;
    Cell cellData[1][1];
    double drug;
    double cycleLength;
} SimState1;



#include "human_vclamp.h"
#include "integrate_rk1_halfdt.h"
#include "main_cell_dts.h"


int main() {
    
    double drug, cycleLength;
    int i, j;
    for( i = 1; i < 6; i+= 1 ){
        drug = 0.5 * (1.0e-6) * i; // uM for Na
        
        for( j = 1; j < 2; j+=20 ){
            
            cycleLength = 400. * j ;
            main_cell_ss_dts( drug, cycleLength  );
            
        }
    }
    
}
