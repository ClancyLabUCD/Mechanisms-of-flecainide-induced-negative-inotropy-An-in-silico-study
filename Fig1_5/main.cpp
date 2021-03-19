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

// for mkdir
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;


typedef struct cell {
    double y0n[112];
    double highestV, lowestV, Apd90, DI, t1, t2;
    double highestdvdt, V90, t3, t4, s2APD, V80, Apd80;
    double dvdt, v, v_new, v_old;
    double currents[2];
    double I_Total, I_inj, I_start, I_end;
    double ui1[4], uj1[4], ui2[4], uj2[4], dun[4], tu[4];
    
    
} Cell;


const int tw = 1;
const int tl = 165;


typedef struct SimState1 {
    double t, tt;
    int tstep, counter, beat;
    Cell *cellData[tw];
    double drug;
    double cycleLength;
    double t1c1, t1c2;
    
    double highestdvdt;
} SimState1;



#include "human_vclamp.h"
#include "integrate_rk1_halfdt.h"
#include "main_1D_line.h"




int main() {
    
    time_t startTime;
    time_t previousTime;
    
    startTime = time(NULL);
    previousTime = startTime;
    
    srand48(startTime);
    
    double drug, cycleLength;
    int i, j;
    
    
    
    omp_set_num_threads(12);
    
    
    
    for( i = 0; i < 6; i+= 1 ){
        drug = 0.5 * (1.0e-6) * i;
        
        
        
        for( j = 150; j < 160; j+=10 ){
            
            cycleLength = 60000. / j;
            main_1D_line( drug, cycleLength  );
            
            
        }
    }
}

