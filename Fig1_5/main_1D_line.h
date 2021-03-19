/*
 *  main_1D.h
 *
 *  Copyright 2010 Colleen Clancy Lab. All rights reserved.
 *
 */

#ifndef main_1D_line_H
#define main_1D_line_H

void main_1D_line( double drug, double cycleLength);

void main_1D_line( double drug, double cycleLength) {
    
    
    SimState1 theState_1;
    theState_1.drug = drug;
    theState_1.cycleLength = cycleLength;
    
    
    double p = 0;//par;  // Parameter array for passing nondefault conditions
    //// Initial conditions
    const double mo = 1.405627e-3;
    const double ho = 9.867005e-1;
    const double jo = 9.915620e-1;
    const double do1 = 7.175662e-6;
    const double fo = 1.000681;
    const double fcaBjo = 2.421991e-2;
    const double fcaBslo = 1.452605e-2;
    const double xtoso = 4.051574e-3;
    const double ytoso = 9.945511e-1;
    const double xtofo = 4.051574e-3;
    const double ytofo = 9.945511e-1;
    const double xkro = 8.641386e-3;
    const double xkso = 5.412034e-3;
    const double RyRro = 8.884332e-1;
    const double RyRoo = 8.156628e-7;
    const double RyRio = 1.024274e-7;
    const double NaBjo = 3.539892;
    const double NaBslo = 7.720854e-1;
    const double TnCLo = 8.773191e-3;
    const double TnCHco = 1.078283e-1;
    const double TnCHmo = 1.524002e-2;
    const double CaMo = 2.911916e-4;
    const double Myoco = 1.298754e-3;
    const double Myomo = 1.381982e-1;
    const double SRBo = 2.143165e-3;
    const double SLLjo = 9.566355e-3;
    const double SLLslo = 1.110363e-1;
    const double SLHjo = 7.347888e-3;
    const double SLHslo = 7.297378e-2;
    const double Csqnbo = 1.242988;
    const double Ca_sro = 0.1e-1; //5.545201e-1;
    const double Najo = 9.06;//8.80329;
    const double Naslo = 9.06;//8.80733;
    const double Naio = 9.06;//8.80853;
    const double Kio = 120;
    const double Cajo = 1.737475e-4;
    const double Caslo = 1.031812e-4;
    const double Caio = 8.597401e-5;
    const double Vmo = -8.09763e+1;
    const double rtoso = 0.9946;
    const double ICajuncinto = 1;
    const double ICaslinto = 0;
    const double C1o = 0.0015;       // []
    const double C2o = 0.0244;       // []
    const double C3o = 0.1494;       // []
    const double C4o = 0.4071;       // []
    const double C5o = 0.4161;       // []
    const double C7o = 0.0001;       // []
    const double C8o = 0.0006;       // []
    const double C9o = 0.0008;       // []
    const double C10o = 0;           // []
    const double C11o = 0;           // []
    const double C12o = 0;           // []
    const double C13o = 0;           // []
    const double C14o = 0;           // []
    const double C15o = 0;           // []
    const double O1o = 0;            // []
    const double O2o = 0;            // []
    const double C6o = 1 - ( C1o + C2o + C3o + C4o + C5o + C7o + C8o + C9o + C10o + C11o + C12o + C13o + C14o + C15o + O1o + O2o );       
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
        DIM1, DIM2, DBC3, DBC2, DBC1, DBO, D_O, D_OS, D_C1, D_C2, D_C3, D_IC3, D_IC2, D_IF, D_IM1, D_IM2, D_BC3, D_BC2, D_BC1, D_BO, 0,
        Ir, C3r, C2r, Okr, C1r, I_drug, C3_drug, C2_drug, O_drug, C1_drug };
    
    SimState1 *S = &theState_1;
    Cell * theCell;
    Cell *Allcells;
    Allcells = new Cell[ tw * tl ];
    
    int ww, ll;
    for( ww = 0; ww < tw; ww++ ) {
        theState_1.cellData[ww] = &Allcells[ ww*tw ];
    }
    
    double * y0n;
    double * y0n2;
    
    
    
    double dt = 1. / 128;
    const int fold = 1. / dt;
    
    const int TotalBeats = 1500;
    
    int iter_max = (cycleLength*TotalBeats) / dt;
    
    
    double dx = 0.01;
    double dy = 0.01;
    
    
    const double Df = 0.00154;
    double dv;
    
    //// Parameters for external modules
    // ECC and CaM modules
    
    time_t startTime;
    time_t currentTime, oldTime;
    
    startTime = time(NULL);
    currentTime = startTime;
    oldTime = startTime;
    
    double cv1 = 0;
    int cvcell1 = 15 ;
    int cvcell2 = (tl-15);
    
    
    char singlestart[200];
    
    sprintf( singlestart, "LastBeat_wBO/S1LastBeatcells_%3.2f_%d", drug*1.e6, (int)cycleLength );
    
    
    theCell = &(S->cellData[0][0]);
    
    FILE *fp3 = fopen(singlestart, "r");
    
    if (fp3 == NULL) {
        
        cout << "From initial file" << endl;
        
        
        for ( ll = 0; ll < tl; ll++ ) {
            for ( int idy = 0; idy < sy0; idy++ ) {
                S->cellData[0][ll].y0n[idy] = y0[idy];
            }
            S->cellData[0][ll].v = y0[38];
            
        }
        
        S->beat = 0;
        S->counter = 0;
        S->t1c1 = 1E10;
        S->t1c2 = 1E10;
        theCell->lowestV = theCell->v;
        theCell->highestV = theCell->v;
        S->highestdvdt = 0;
        //theCell->dvdt = 0;
        theCell->V90 = 0;
        theCell->t1 = 0;
        theCell->t2 = 0;
        
    } else {
        cout << "from single states " << endl;
        
        theCell = &(S->cellData[0][0]);
        
        fread( theCell, sizeof(Cell), 1, fp3);
        for ( ll = 0; ll < tl; ll++ ) {
            for ( int idy = 0; idy < sy0; idy++ ) {
                S->cellData[0][ll].y0n[idy] = y0[idy];
                
            }
        }
        
        
        fclose( fp3 );
        
        S->beat = 0;
        S->counter = 0;
        S->t1c1 = 1E99;
        S->t1c2 = 1E99;
        theCell->lowestV = theCell->v;
        theCell->highestV = theCell->v;
        S->highestdvdt = 0;
        //theCell->dvdt = 0;
        theCell->V90 = 0;
        theCell->t1 = 0;
        theCell->t2 = 0;
    }
    
    char OutputFolder[690] ;
    
    
    sprintf( OutputFolder,"caseW1_5targets_%3.2f_%d" , drug*1.e6, (int)cycleLength );
    
    
    
    mkdir( OutputFolder, 0777 );
    
    char namefpy[690];
    
    sprintf( namefpy, "%s/%s%3.2f%s", OutputFolder, "VmResult_", drug*1.e6, ".txt" );
    
    FILE *fpy = fopen(namefpy, "w");
    
    char namefpyca[690];
    
    sprintf( namefpyca, "%s/%s%3.2f%s", OutputFolder, "CaResult_", drug*1.e6, ".txt" );
    
    FILE *fpyca = fopen(namefpyca, "w");
    
    char nameecg[690];
    sprintf( nameecg, "%s/%s%3.2f%s", OutputFolder, "ecg_", drug*1.e6, ".txt" );
    
    FILE *ecg = fopen(nameecg, "w");
    
    char namefpcv[690];
    sprintf( namefpcv, "%s/%s%3.2f%s", OutputFolder, "CV_", drug*1.e6,  ".txt" );
    
    FILE *fpcv= fopen(namefpcv, "w");
    
    double tspan[2];
    
    double t_iter;
    
    double currents[10];
    

    
    double I_inj;
    
    int beatadded = 0;
    
    for( int iter = 1; iter <= iter_max; iter ++ ) {
        
        S->tt = iter * dt;
        
        S->t = fmod( S->tt, cycleLength );
        
        S->counter ++;
        
#pragma omp parallel default(none) private( ll, theCell, y0n , I_inj, currents ) shared( cout , S, dt, tspan, p, sy0 )
        {
#pragma omp for schedule(static, 1)
            
            
            for (ll = 0; ll < tl; ll+=1 ) {
                
                if( S->t < 0.5 && ( ll == 0 || ll == 1) ) {
                    I_inj = -300;
                } else {
                    I_inj=0.0;
                }
                
                
                theCell = &(S->cellData[0][ll]);
                y0n = theCell->y0n;
                
                theCell->v_old = y0n[38];
                theCell->v = y0n[38];
                
                integrate_rk1_halfdt( f, tspan, sy0, y0n, dt, p, currents , S );
                
                theCell->I_Total = (theCell->v_old - y0n[38]) / dt + I_inj;
                
                
            }
            
            
        } // end omp
        
        
#pragma omp parallel default(none) private( ll, theCell, dv) shared( cout , S, dt, sy0, dx ,dy, cvcell1, cvcell2  )
        {
#pragma omp for // schedule(static, 1)
            
            for ( ll = 0; ll < tl; ll++ ) {
                
                
                
                theCell = &(S->cellData[0][ll]);
                
                if(ll>0 && ll<(tl-1)) {
                    dv = dt * (-theCell->I_Total + Df * (S->cellData[0][ll-1].v - 2*theCell->v + S->cellData[0][ll+1].v) / (dx*dx));
                }
                else if (ll==0) {
                    dv = dt * (-theCell->I_Total + Df * (-theCell->v + S->cellData[0][1].v) / (dx*dx));
                }
                else if (ll==(tl-1)) {
                    dv=dt*(-theCell->I_Total + Df * ( -theCell->v + S->cellData[0][ll-1].v) / (dx*dx));
                }
                
                
                theCell->dvdt = fabs( dv / dt );
                theCell->y0n[38] = theCell->v_old + dv;
                theCell->v_new = theCell->y0n[38];
                
                
                if ( ll == cvcell1) {
                    theCell->v = theCell->y0n[38];
                    
                    if ( theCell->v > -40 ) {
                        
                        if (S->t1c1 > S->tt ) {
                            S->t1c1 = S->tt; // 1st beat arrived.
                            cout << cvcell1 << "th beat arrived (@ ms): "  << S->t1c1 << endl;
                        }
                        
                    }
                    
                  
                    if ( theCell->dvdt > S->highestdvdt ) {
                        S->highestdvdt = theCell->dvdt;
                    }
                }
                else if ( ll == cvcell2 ) {
                    
                    
                    theCell->v = theCell->y0n[38];
                    
                    if ( theCell->v > -40 ) {
                        
                        if (S->t1c2 > S->tt ) {
                            S->t1c2 = S->tt; // 2nd beat arrived.
                            cout << cvcell2 << "th beat arrived (@ ms): " << S->t1c2 << endl;
                        }
                        
                    }
                }
                
            }
            
            
            
        } // end openmp
        
        
        cv1 = dx * 1000* ( cvcell2 - cvcell1 ) / ( S->t1c2 - S->t1c1 ) ; //cm/s
        
        
        
        if( S->cellData[0][0].v > -40 && beatadded == 0 ) {
            S->beat += 1;
            cout << "Now beat = " << S->beat << endl;
            beatadded = 1;
            
            cout <<"highestdvdt in " << cvcell1 << "th cell = " << S->highestdvdt  << endl;

            cout << "CV = " << cv1 << " cm/s " << endl;

            fprintf (fpcv, "%d\t%8.5f\t%8.5f\n", S->beat, cv1, S->highestdvdt);

            S->t1c1 = 1E99;
            S->t1c2 = 1E99;
            S->highestdvdt = 0;

        }
        
        if( S->cellData[0][0].v <= -40 ) {
            beatadded = 0;
        }
        
        
       
        
        
        if ( (S->counter % 10000 ) == 0 ) {
            
            cout << "tt: " << S->tt  << ", runtime: " << (time(NULL) - startTime)/60 << " min " << endl;
            
        }
        
        
        if ( (S->counter % 100 ) == 0  && S->beat >= (TotalBeats - 100)) {
            
            
            fprintf (fpy, "%8.6e\t", S->tt);
            fprintf (fpyca, "%8.6e\t", S->tt);
            
            for (ll=0; ll< tl; ll+=1){
                fprintf(fpy, "%8.6f\t", S->cellData[0][ll].y0n[38]);
                fprintf(fpyca, "%8.6f\t", S->cellData[0][ll].y0n[37]);
            }
            
            fprintf(fpy, "\n");
            fprintf(fpyca, "\n");
            
            double Ev=0;
            double r0 = tl * 0.01 / 2 / 3.14159;
            double theta1, theta2, xll1, yll1, xll2, yll2;
            
            for (ll=15; ll<=(tl-15); ll++){
                                Ev = ( Ev
                                      + ( S->cellData[0][ll].y0n[38] - S->cellData[0][ll+1].y0n[38] )
                                      * ( 1. / ( ( tl - ll - 1 ) * dx + 2 )
                                         - 1. / ( ( tl - ll ) * dx + 2 ) ) );
                

            }

            fprintf(ecg, "%10.3f\t%10.6f\n", S->tt, Ev);
        }
        
        //
        
        
    }  // end t_iter
    
    
    fclose (fpy);
    fclose (ecg);
    fclose (fpyca);
    fclose (fpcv);
    
    theCell = &S->cellData[0][0];
    

    
    char namefpy2[690];
    
    sprintf( namefpy2, "%s/%s%3.2f%s", OutputFolder, "S1LastBeat_", drug*1.e6, ".txt" );
    
    FILE *fpy2 = fopen(namefpy2, "w");
    
  
    fwrite( Allcells, sizeof( Cell ), tw*tl, fpy2 );
    fclose(fpy2);
    
   
    
    delete[] Allcells;
    
    
}


#endif
