 /*
 *  main_cell_dts.h
 *
 *  Copyright 2020 Colleen Clancy Lab. All rights reserved.
 *
 */

#ifndef main_cell_ss_dts_H
#define main_cell_ss_dts_H


void main_cell_ss_dts( double drug, double cycleLength );

void main_cell_ss_dts( double drug, double cycleLength ){
    
    SimState1 theState_1;
    theState_1.drug = drug;
    theState_1.cycleLength = cycleLength;
    
    double p = 0;//par;  // Parameter array for passing nondefault conditions
    //// Initial conditions
    const double mo=1.405627e-3;
    const double ho= 9.867005e-1;
    const double jo=9.915620e-1;
    const double do1=7.175662e-6;
    const double fo=1.000681;
    const double fcaBjo=2.421991e-2;
    const double fcaBslo=1.452605e-2;
    const double xtoso=4.051574e-3;
    const double ytoso=9.945511e-1;
    const double xtofo=4.051574e-3;
    const double ytofo= 9.945511e-1;
    const double xkro=8.641386e-3;
    const double xkso= 5.412034e-3;
    const double RyRro=8.884332e-1;
    const double RyRoo=8.156628e-7;
    const double RyRio=1.024274e-7;
    const double NaBjo=3.539892;
    const double NaBslo=7.720854e-1;
    const double TnCLo=8.773191e-3;
    const double TnCHco=1.078283e-1;
    const double TnCHmo=1.524002e-2;
    const double CaMo=2.911916e-4;
    const double Myoco=1.298754e-3;
    const double Myomo=1.381982e-1;
    const double SRBo=2.143165e-3;
    const double SLLjo=9.566355e-3;
    const double SLLslo=1.110363e-1;
    const double SLHjo=7.347888e-3;
    const double SLHslo=7.297378e-2;
    const double Csqnbo= 1.242988;
    const double Ca_sro=0.1e-1;
    const double Najo=9.136;
    const double Naslo=9.136;
    const double Naio=9.136;
    const double Kio=120;
    const double Cajo=1.737475e-4;
    const double Caslo= 1.031812e-4;
    const double Caio=8.597401e-5;
    const double Vmo=-8.09763e+1;
    const double rtoso=0.9946;
    const double ICajuncinto=1;
    const double ICaslinto=0;
    const double C1o=0.0015;       //  []
    const double C2o=0.0244;       //  []
    const double C3o=0.1494;       //  []
    const double C4o=0.4071;       //  []
    const double C5o=0.4161;       //  []
    const double C7o=0.0001;       //  []
    const double C8o=0.0006;       //  []
    const double C9o=0.0008;       //  []
    const double C10o=0;           //  []
    const double C11o=0;           //  []
    const double C12o=0;           //  []
    const double C13o=0;           //  []
    const double C14o=0;           //  []
    const double C15o=0;           //  []
    const double O1o=0;            //  []
    const double O2o=0;            //  []
    const double C6o=1-(C1o+C2o+C3o+C4o+C5o+C7o+C8o+C9o+C10o+C11o+C12o+C13o+C14o+C15o+O1o+O2o);       //  []
    //WT Channel Markov Initial Conditions
    //Drug Free States of WT Channel
    const double IC3 = 0;
    const double IC2 = 0;
    const double IF = 0;
    const double IM1 = 0;
    const double IM2 = 0;
    const double C3 = 1;
    const double C2 = 0;
    const double O = 0;
    const double OS =0;
    const double BC3 = 0;
    const double BC2 = 0;
    const double BC1 = 0;
    const double BO = 0;
    
    //Charged Drug States of WT Channel
    const double DIC3 = 0;
    const double DIC2 =0;
    const double DIF = 0;
    const double DIM1 = 0;
    const double DIM2 = 0;
    const double DC3 = 0;
    const double DC2 =0;
    const double DC1 = 0;
    const double DO = 0;
    const double DOS = 0;
    const double DBC3 = 0;
    const double DBC2 = 0;
    const double DBC1 = 0;
    const double DBO = 0;
    //Neutral Drug States of WT Channel
    const double D_IC3 = 0;
    const double D_IC2 =0;
    const double D_IF = 0;
    const double D_IM1 = 0;
    const double D_IM2 = 0;
    const double D_C3 = 0;
    const double D_C2 =0;
    const double D_O = 0;
    const double D_OS = 0;
    const double D_C1 = 0;
    const double D_BC3 = 0;
    const double D_BC2 = 0;
    const double D_BC1 = 0;
    const double D_BO = 0;
    
    const double C1 = 1- ( O +  OS +  C3 +  C2 +  IC3 +  IC2 +  IF +  IM1 +  IM2 +  BC3 +  BC2 +  BC1 +  BO +
                          DO +  DOS +  DC1 +  DC2 +  DC3 +  DIC3 +  DIC2 +  DIF +  DIM1 +  DIM2 +  DBC3 +  DBC2 +  DBC1 +  DBO +
                          D_O +  D_OS +  D_C1 +  D_C2 +  D_C3 +  D_IC3 +  D_IC2 +  D_IF +  D_IM1 +  D_IM2 +  D_BC3 +  D_BC2 +  D_BC1 +  D_BO);
    // Markov Ikr
    const double Ir = 0;
    const double C3r = 0;
    const double C2r = 0;
    const double Okr = 0;
    
    const double I_drug = 0;
    const double C3_drug = 0;
    const double C2_drug = 0;
    const double O_drug = 0;
    const double C1_drug = 0;
    
    const double C1r = 1 - (Ir + C3r + C2r + Okr + I_drug + C3_drug + C2_drug + O_drug + C1_drug);
    
    const double RyR_drug = 0;
    
    // initial conditions for IKur
    const double rkuro = 0;
    const double skuro = 1.0;
   
    
    
    int sy0 = 112;
    
    double y0[] = { mo, ho, jo, do1, fo, fcaBjo, fcaBslo, xtoso, ytoso, xtofo,
        ytofo, xkro, xkso, RyRro, RyRoo, RyRio, NaBjo, NaBslo, TnCLo, TnCHco,
        TnCHmo, CaMo, Myoco, Myomo, SRBo, SLLjo, SLLslo, SLHjo, SLHslo, Csqnbo,
        Ca_sro, Najo, Naslo, Naio, Kio, Cajo, Caslo, Caio, Vmo, rtoso,
        1, C1o, C2o, C3o, C4o, C5o, C6o, C7o, C8o, C9o,
        C10o, C11o, C12o, C13o, C14o, C15o, O1o, rkuro, skuro, O, OS, C3, C2, C1, IC3, IC2, IF, IM1, IM2, BC3, BC2, BC1, BO, DO, DOS, DC1, DC2, DC3, DIC3, DIC2, DIF,
        DIM1, DIM2, DBC3, DBC2, DBC1, DBO, D_O, D_OS, D_C1, D_C2, D_C3, D_IC3, D_IC2, D_IF, D_IM1, D_IM2, D_BC3, D_BC2, D_BC1, D_BO, RyR_drug,
        Ir, C3r, C2r, Okr, C1r , I_drug, C3_drug, C2_drug, O_drug, C1_drug};
    
    SimState1 *S = &theState_1;
    Cell * theCell;
    int ww, ll;
    double * y0n;
    double * y0n2;

    S->beat = 0;
    theCell = &(S->cellData[0][0]);
    
    const char *RestingStates = "restingcells";

    FILE *fp3 = fopen(RestingStates, "r");
    
    if (fp3 == NULL) {
        
        cout << "From initial file" << endl;
        
    for ( int idy = 0; idy < sy0; idy++ ) {
        S->cellData[0][0].y0n[idy] = y0[idy];
        cout << idy << "\t" << y0[idy] << endl;
    }

        
    } else {
        
        cout << "From initial single cell file" << endl;
        fread( theCell, sizeof(Cell), 1, fp3 );
        fclose( fp3 );
    }

 
    
    const int fold = 128;
    double dt =  1. / fold;
    
    const int TotalBeats = 1500;
    double t_max = cycleLength * TotalBeats;
    

    double dv;
    
    //// Parameters for external modules
    // ECC and CaM modules
    
    time_t startTime;
    time_t currentTime, oldTime;
    
    startTime = time(NULL);
    currentTime = startTime;
    oldTime = startTime;
    
   
    char namefpy[200];
    
    sprintf( namefpy, "S1Result_%3.2f_%d.txt", drug*1.e6, (int)cycleLength );
    
    FILE *fpy;
    
    fpy = fopen(namefpy, "w");
    
    double tspan[2];
    
    double t_iter;
    
    
    
    double currents[10];
    
    for( t_iter = 0; t_iter < t_max; t_iter += dt ) {
        
        
        tspan[0] = t_iter ;
        tspan[1] = tspan[0] + dt ;
        
        
        y0n = theCell->y0n;
        
        
        theCell->v_old = y0n[38];
        theCell->v = y0n[38];
        
        integrate_rk1_halfdt( f, tspan, sy0, y0n, dt, p, currents , S );
        
        theCell->I_Total = (theCell->v_old - y0n[38]) / dt;
        
        
        dv=dt*(-theCell->I_Total);
        
        theCell->y0n[38] = theCell->v_old + dv;
        theCell->v_new = theCell->y0n[38];
        
        
        S->tt = tspan[1];
        S->t = fmod( S->tt, cycleLength );
        
        if( S->t == dt ) {
            S->beat += 1;
            currentTime = time(NULL);
            cout << "Now beat = " << S->beat << endl;
            cout << "Time Elapsed from previous beat = " << currentTime - oldTime << endl;
            cout << "Total Time Elapsed = " << currentTime - startTime << endl;

            oldTime = currentTime;
        }
        
        double dvdt = -theCell->I_Total; // dv / dt;
        
        // Find APD80 in the last beat.
        if ( S->beat == TotalBeats ){
            theCell = &S->cellData[0][0];
            if ( S->t == dt ) {
                theCell->lowestV = theCell->v;
                theCell->highestV = theCell->v;
                theCell->highestdvdt = 0;
                theCell->V80 = 0;
                theCell->t1 = 0;
                theCell->t2 = 0;
                
            }
            
            if ( theCell->highestdvdt < dvdt ) {
                
                theCell->highestdvdt = dvdt;
                theCell->t1 = S->t;
            }
            if ( (S->t < 10) & (theCell->highestV <= theCell->v_new) ) {
                theCell->highestV = theCell->v_new;
                theCell->V80 = theCell->highestV - 0.8 * ( theCell->highestV - theCell->lowestV );
               
            } else if ( theCell->v_new > theCell->V80 ) {
                theCell->t2 = S->t;
                theCell->Apd80 = S->t - theCell->t1;
                
            }
        }
        
        
        // Output Voltage to a file
        double N2, dui, duj;
        
        theCell = &(S->cellData[0][0]);
        y0n = theCell->y0n;
        
        double voltage = y0n[38];
        double cytoCa = y0n[37];
        double caSR = y0n[30];
        double caSRT = y0n[29]+y0n[30];
        double nai = y0n[33];
        theCell->dvdt = ( voltage - theCell->v_old ) / dt;
        
        // Computing direction change
        N2 = pow( ( 1. + theCell->dvdt * theCell->dvdt ), 0.5 );
        theCell->ui2[0] = 1. / N2;
        theCell->uj2[0] = theCell->dvdt / N2;
        dui = theCell->ui1[0] - theCell->ui2[0];
        duj = theCell->uj1[0] - theCell->uj2[0];
        theCell->dun[0] += pow( ( dui * dui + duj * duj ), 0.5 );
        theCell->ui1[0] = theCell->ui2[0];
        theCell->uj1[0] = theCell->uj2[0];
        
        
        //
        // Save if there is an significant change of direction or end of a cycle.
        if( S->beat > 0   && ( ( ( exp( 0.1 * ( S->tt - theCell->tu[0] ) ) * theCell->dun[0] ) > 0.001 ) || ( S->t < dt ) ) ) {
            
            fprintf(fpy, "%16.8f\t%12.9f\t%12.9f\t%12.9f\t%12.9f\t%12.9f\t", S->tt , voltage, cytoCa, caSR, caSRT, nai );
            for( int i = 0; i < 10; i++ ) {
                fprintf( fpy, "%12.9f\t", currents[i] );
            }
            fprintf( fpy, "\n" );
            
            theCell->dun[0] = 0.;
            theCell->tu[0] = S->tt;
            
        }
        
        
        
        
    }  // end t_iter
    
    theCell = &S->cellData[0][0];
    
    cout << "Apd80 = " << theCell->Apd80 << endl;
    
    char namefpy1[200];
    
    sprintf( namefpy1, "S1LastBeatcells_%3.2f_%d", drug*1.e6, (int)cycleLength );
    FILE *fp2 = fopen( namefpy1, "w");


    
    fwrite( theCell , sizeof(Cell), 1, fp2);
    fclose(fp2);
    
    char namefpy2[200];
    sprintf( namefpy2, "S1Apd80_%3.2f_%d.txt", drug*1.e6, (int)cycleLength );
    fp2 = fopen( namefpy2, "w");
    

    fprintf( fp2, "%12.8f\t%12.8f\t%12.8f\t", theCell->Apd80, theCell->lowestV, theCell->highestV );
    fclose(fp2);
  
    
}

#endif
