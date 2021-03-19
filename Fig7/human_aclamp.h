#ifndef f_H
#define f_H

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
// #include <omp.h>

using namespace std;

void f( double dt0, double t, double * y, int p, double * ydot, double * currents , SimState1 *S);
double mod( double a, double b);

double mod( double a, double b) {
	double c = a / b;
	return ( ( c - floor( c ) ) * b );
}

void f( double dt0, double t, double * y, int p, double * ydot, double * currents , SimState1 *S){
   
    double cycleLength = S->cycleLength;

    //// Model Parameters
    //// EPI or ENDO?
    const int epi=1;
    //// AF
    const double AF=0;
    //// ISO
    const double ISO=0;
    //// Right ATRIUM
    const double RA=1;
    
    const double Ach = 0.0;
    
    double drug_onNa = S->drug;
    
    double drug_onKr = S->drug;

    double drug_onRyR = S->drug;

    double drug_onIto = S->drug;
    
    double drug_onCaL = S->drug;
    
    
    // Constants
    const double R = 8314.;       // [J/kmol*K]
    const double Frdy = 96485.;   // [C/mol]
    const double Temp = 310.;     // [K]
    const double FoRT = Frdy / R / Temp;
    const double Cmem = 1.1e-10;   // [F] membrane capacitance 1.3810e-10;//
    const double Qpow = ( Temp - 310. ) / 10. ;
    
    const double pi = 3.141592653589793;
    
    // Cell geometry
    const double cellLength = 100.;     // cell length [um]113;//100
    const double cellRadius = 10.25;   // cell radius [um]12;//10.25
    const double junctionLength = 160e-3;  //  junc length [um]
    const double junctionRadius = 15e-3;   //  junc radius [um]
    const double distSLcyto = 0.45;    //  dist. SL to cytosol [um]
    const double distJuncSL = 0.5;  //  dist. junc to SL [um]
    const double DcaJuncSL = 1.64e-6;  //  Dca junc to SL [cm^2/sec]
    const double DcaSLcyto = 1.22e-6; //  Dca SL to cyto [cm^2/sec]
    const double DnaJuncSL = 1.09e-5;  //  Dna junc to SL [cm^2/sec]
    const double DnaSLcyto = 1.79e-5;  //  Dna SL to cyto [cm^2/sec]
    const double Vcell = pi * cellRadius * cellRadius * cellLength * 1e-15;    //  [L]
    const double Vmyo = 0.65 * Vcell;
    const double Vsr = 0.035 * Vcell;
    const double Vsl = 0.02 * Vcell;
    const double Vjunc = 1. * 0.0539 * .01 * Vcell;
    const double SAjunc = 20150. * pi * 2. * junctionLength * junctionRadius;  //  [um^2]
    const double SAsl = pi * 2. * cellRadius * cellLength;          //  [um^2]
    
    const double J_ca_juncsl = 1. / 1.2134e12; //  [L/msec] = 8.2413e-13
    const double J_ca_slmyo = 1. / 2.68510e11; //  [L/msec] = 3.2743e-12
    const double J_na_juncsl = 1. / ( 1.6382e12 / 3. * 100. ); //  [L/msec] = 6.1043e-13
    const double J_na_slmyo = 1. / ( 1.8308e10 / 3. * 100. );  //  [L/msec] = 5.4621e-11
    
    //  Fractional currents in compartments
    const double Fjunc = 0.11;
    const double Fsl = 1. - Fjunc;
    const double Fjunc_CaL = 0.9;
    const double Fsl_CaL = 1. - Fjunc_CaL;
    
    //  Fixed ion concentrations
    const double Cli = 15.;   //  Intracellular Cl  [mM]
    const double Clo = 150.;  //  Extracellular Cl  [mM]
    const double Ko = 5.4;   //  Extracellular K   [mM]
    const double Nao = 140.;  //  Extracellular Na  [mM]
    const double Cao = 1.8;  //  Extracellular Ca  [mM]
    const double Mgi = 1.;    //  Intracellular Mg  [mM]
    
    //  Nernst Potentials
    const double ena_junc = ( 1. / FoRT ) * log( Nao / y[31] );     //  [mV]
    const double ena_sl = ( 1. / FoRT ) * log( Nao / y[32] );       //  [mV]
    const double ek = ( 1. / FoRT ) * log( Ko / y[34]);	        //  [mV]
    const double eca_junc = ( 1. / FoRT / 2. ) * log( Cao / y[35] );   //  [mV]
    const double eca_sl = ( 1. / FoRT / 2. ) * log( Cao / y[36] );     //  [mV]
    const double ecl = ( 1. / FoRT ) * log( Cli / Clo );            //  [mV]
    
    //  Na transport parameters
    
    // //
    const double GNa = 16. * ( 1. - 0.1 * AF ) * 0.8;  //  [mS/uF]
    const double GNaB = 0.597e-3;    //  [mS/uF]
    const double IbarNaK = 1.26 *1.2;     //  [uA/uF]
    const double KmNaip = 11. * ( 1. - 0.25 * ISO );         //  [mM]11
    const double KmKo = 1.5;         //  [mM]1.5
    const double Q10NaK = 1.63;
    const double Q10KmNai = 1.39;
    
    // //  K current parameters
    const double pNaK = 0.01833;
    const double gkp = 0.002;
    
    //  Cl current parameters
    const double GClCa = 0.0548;   //  [mS/uF]
    const double GClB = 9e-3;        //  [mS/uF]
    const double KdClCa = 100e-3;    //  [mM]
    
    //  I_Ca parameters
    double IC50_flecCaL = 27.1e-6;//10e-6;//20e-6;
    double factor_FlecCaL = 1/(1+(drug_onCaL/IC50_flecCaL));
    
    const double pNa = ( 1. + 0.5 * ISO ) * ( 1. - 0.5 * AF ) * 0.75e-8 * 0.8 * factor_FlecCaL;       //  [cm/sec]
    const double pCa = ( 1. + 0.5 * ISO ) * ( 1. - 0.5 * AF ) * 2.7e-4 * 0.8 * factor_FlecCaL;       //  [cm/sec]
    const double pK = ( 1. + 0.5 * ISO ) * ( 1. - 0.5 * AF ) * 1.35e-7 * 0.8 * factor_FlecCaL;        //  [cm/sec]
    const double Q10CaL = 1.8;
    
    // //  Ca transport parameters
    const double IbarNCX = ( 1. + 0.4 * AF ) * 3.15;      //  [uA/uF]5.5 before - 9 in rabbit
    const double KmCai = 3.59e-3;    //  [mM]
    const double KmCao = 1.3;        //  [mM]
    const double KmNai = 12.29;      //  [mM]
    const double KmNao = 87.5;       //  [mM]
    const double ksat = 0.27;        //  [none]
    const double nu = 0.35;          //  [none]
    const double Kdact = 0.384e-3;   //  [mM] 0.256 rabbit
    const double Q10NCX = 1.57;      //  [none]
    const double IbarSLCaP =  0.0471; //  IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
    const double KmPCa = 0.5e-3;     //  [mM]
    const double GCaB = 6.0643e-4;    //  [uA/uF] 3
    const double Q10SLCaP = 2.35;    //  [none]
    
    //  SR flux parameters
    const double Q10SRCaP = 2.6;          //  [none]
    const double Vmax_SRCaP = 5.3114e-3;  //  [mM/msec] (286 umol/L cytosol/sec)
    const double Kmf = ( 2.5 - 1.25 * ISO ) * 0.246e-3;          //  [mM] default
    const double Kmr = 1.7;               //  [mM]L cytosol
    const double hillSRCaP = 1.787;       //  [mM]
    const double ks = 25.;                 //  [1/ms]
    const double koCa = 10. + 20. * AF + 10. * ISO * ( 1. - AF );               //  [mM^-2 1/ms]   // default 10   modified 20
    const double kom = 0.06;              //  [1/ms]
    const double kiCa = 0.5;              //  [1/mM/ms]
    const double kim = 0.005;             //  [1/ms]
    const double ec50SR = 0.45;           //  [mM]
    
    //  Buffering parameters
    //  koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
    const double Bmax_Naj = 7.561;       //  [mM] //  Na buffering
    const double Bmax_Nasl = 1.65;       //  [mM]
    const double koff_na = 1e-3;         //  [1/ms]
    const double kon_na = 0.1e-3;        //  [1/mM/ms]
    const double Bmax_TnClow = 70e-3;    //  [mM]                      //  TnC low affinity
    const double koff_tncl = ( 1. + 0.5 * ISO ) * 19.6e-3;    //  [1/ms]
    const double kon_tncl = 32.7;        //  [1/mM/ms]
    const double Bmax_TnChigh = 140e-3;  //  [mM]                      //  TnC high affinity
    const double koff_tnchca = 0.032e-3; //  [1/ms]
    const double kon_tnchca = 2.37;      //  [1/mM/ms]
    const double koff_tnchmg = 3.33e-3;  //  [1/ms]
    const double kon_tnchmg = 3e-3;      //  [1/mM/ms]
    const double Bmax_CaM = 24e-3;       //  [mM] **? about setting to 0 in c-code**   //  CaM buffering
    const double koff_cam = 238e-3;      //  [1/ms]
    const double kon_cam = 34.;           //  [1/mM/ms]
    const double Bmax_myosin = 140e-3;   //  [mM]                      //  Myosin buffering
    const double koff_myoca = 0.46e-3;   //  [1/ms]
    const double kon_myoca = 13.8;       //  [1/mM/ms]
    const double koff_myomg = 0.057e-3;  //  [1/ms]
    const double kon_myomg = 0.0157;     //  [1/mM/ms]
    const double Bmax_SR = 19. * .9e-3;     //  [mM] (Bers text says 47e-3) 19e-3
    const double koff_sr = 60e-3;        //  [1/ms]
    const double kon_sr = 100.;           //  [1/mM/ms]
    const double Bmax_SLlowsl = 37.4e-3 * Vmyo / Vsl;        //  [mM]    //  SL buffering
    const double Bmax_SLlowj = 4.6e-3 * Vmyo / Vjunc * 0.1;    //  [mM]    // Fei *0.1!!! junction reduction factor
    const double koff_sll = 1300e-3;     //  [1/ms]
    const double kon_sll = 100.;          //  [1/mM/ms]
    const double Bmax_SLhighsl = 13.4e-3 * Vmyo / Vsl;       //  [mM]
    const double Bmax_SLhighj = 1.65e-3 * Vmyo / Vjunc * 0.1;  //  [mM] // Fei *0.1!!! junction reduction factor
    const double koff_slh = 30e-3;       //  [1/ms]
    const double kon_slh = 100.;          //  [1/mM/ms]
    const double Bmax_Csqn = 140e-3 * Vmyo / Vsr;            //  [mM] //  Bmax_Csqn = 2.6;      //  Csqn buffering
    const double koff_csqn = 65.;         //  [1/ms]
    const double kon_csqn = 100.;         //  [1/mM/ms]
    
    
    // //  Membrane Currents

    ydot[0] = 0;
	ydot[1] = 0;
	ydot[2] = 0;
  


    double Q10 = 3. ;
	
	double Tfactor = 1.0 / pow( Q10 , ( ( 37.0 - ( Temp - 273 ) ) / 10.0 ) );
	double pH = 7.4;
	double pKa = 9.3;
	
	double portion = 1.0 / ( 1 + pow( 10 , ( pH - pKa ) ) );
	double diffusion = 5500.;
	
    
	double drug_charged = drug_onNa * portion;
	double drug_neutral = drug_onNa * ( 1 - portion );
	double dd = -0.7;
	
	
	////Transition Rates for WT Channel
	double a11= Tfactor*8.5539/(7.4392e-2*exp(-y[38]/17.0)+ 2.0373e-1*exp(-y[38]/150));
	double a12= Tfactor*8.5539/(7.4392e-2*exp(-y[38]/15.0)+ 2.0373e-1*exp(-y[38]/150));
	double a13= Tfactor*8.5539/(7.4392e-2*exp(-y[38]/12.0)+ 2.0373e-1*exp(-y[38]/150));
	double b11= Tfactor*7.5215e-2*exp(-y[38]/20.3);
	double b12= Tfactor*2.7574*exp(-(y[38]-5)/20.3);
	double b13= Tfactor*4.7755e-1*exp(-(y[38]-10)/20.3);
	
	double a3 = Tfactor*5.1458e-6*exp(-y[38]/8.2471);
	double b3=Tfactor*6.1205*exp((y[38])/13.542);
	
	
	double a2= Tfactor*(13.370*exp(y[38]/43.749));
	double b2= ((a13*a2*a3)/(b13*b3));
	
	double a4 = 0*a2;
	double b4 = 0*a3;
	double a5= 0*a2;
	double b5 = 0*a3;
	
    double mu1;
	double mu2;
	
	
    mu1 = 2.0462e-7; mu2 = 8.9731e-4;

	
	
	double ax = 3.4229e-2*a2;
	double bx = 1.7898e-2*a3;
	
	////Flec Transition Rates
	double ax1 = 5.7839e-05 * ax;
	double bx1 =  1.6689e-08* bx;
	double a13c = 3.6324e-03 *a13;
	double a22 = 1.4847e+03 *a2;
	double b33 =  1.7352e-06* b3;
	double a33 = 6.7505e-05 * a3;
	double a44 =  2.4135e+00* a2;
	double b44 =  4.9001e-02* a3;
	double a55 = 0;
    double b55 = 0;
	
	
	
	double ax2 = 2.6126e-01 * ax;
	double a13n = 2.6452e+00 * a13;
	double a_22 =  4.2385e+01 * a2;
	double b_33 = 2.1181e+00 * b3;
	double a_44 =  1.0326e-03 * a2;
	double b_44 = 2.1378e-02 * a3;
	
	double a_55 = 0;
	double b_55 = 0;
	const double kd0 = 11.2*(1e-6);
	double kd_open=kd0*exp( (dd* y[38] * Frdy) /(R * Temp));
    
	
	// charged drug
	double kon=drug_charged*diffusion;
	double koff=kd_open*diffusion;
	double kcoff = koff;
	double kcon = kon;
	
    //bursting drug
    double kbon = kon;
    double kboff = 95.2165 *(1e-6)  *   exp( (dd*y[38]*Frdy) /(R*Temp))   *diffusion;	//This expression gives at 30uM, 70% block at -100mV Nagamoto, 2000
    double kcbon = kbon;
    double kcboff = kboff;
    
  
    double b13c, b22;
    
	if (drug_onNa ==0 || drug_charged ==0 ){b13c = 0;}
	else{b13c = (b13*kcon*koff*a13c)/(kon*kcoff*a13);}
	
	if (b13c ==0){b22 = 0;}
	else {b22=(a13c*a22*a33)/(b13c*b33);}
	
    
    
	// neutral drug
	double k_on = drug_neutral*diffusion;
	double k_off=400*(1e-6)*diffusion;
	double ki_on=k_on/2;
	double ki_off=5.4*(1e-6)*diffusion;
	double kc_on=k_on/2;
	double kc_off=800*(1e-6)*diffusion;
	
    double a_33, b13n, b_22, bx2;
    
	if (drug_onNa ==0 || drug_neutral ==0 ){a_33 = 0;}
	else {a_33 = (ki_off*a3*kc_on*b_33)/(ki_on*kc_off*b3);}
	
	if (drug_onNa ==0 || drug_neutral ==0){b13n = 0;}
	else {b13n = (b13*kc_on*a13n*k_off)/(kc_off*a13*k_on);}
	
	if (b13n==0){b_22 =0;}
	else {b_22 = (a_33*a13n*a_22)/(b_33*b13n);}
	
	if (drug_onNa ==0 || drug_neutral ==0){bx2 = 0;}
	else{bx2 = (bx*k_on*ax2*ki_off)/(ax*ki_on*k_off);}
	

	
	const double coef_O = (b13 + a2  + mu1 + kon + k_on + ax);
	const double coef_C1 = (b12 + b3 + a13  + mu1 + kcon + kc_on);
	const double coef_C2 = (b11 + b3 + a12  + mu1 + kcon + kc_on);
	const double coef_C3 = (b3 + a11  + mu1 + kcon + kc_on);
	const double coef_IC3 = (a11 + a3 + ki_on);
	const double coef_IC2 = (b11 + a3 + a12 + ki_on);
	const double coef_IF = (b12 + b2 + a3 + a4 + ki_on);
	const double coef_IM1 = (b4 + a5);
	const double coef_IM2 = b5;
	const double coef_OS = (bx + ki_on);
	
	const double coef_BO = (mu2 + b13 + kbon + k_on);
	const double coef_BC3 = (mu2 + a11 + kcbon + kc_on);
	const double coef_BC2 = (mu2 + b11 + a12 + kcbon + kc_on);
	const double coef_BC1 = (mu2 + b12 + a13 + kcbon + kc_on);
	
	const double coef_DO = (koff + b13c + a22 + ax1 + mu1);
	const double coef_DC1 = (kcoff + b12 + b33 + a13c + mu1);
	const double coef_DC2 = (kcoff + b11 + b33 + a12 + mu1);
	const double coef_DC3 = (kcoff+ b33 + a11 + mu1);
	const double coef_DOS = bx1;
	const double coef_DIC3 = (a11 + a33);
	const double coef_DIC2 = (a33 + b11 + a12);
	const double coef_DIF = (a33 + b12 + a44 + b22);
	const double coef_DIM1 = ( b44 + a55 );
	const double coef_DIM2 = b55 ;
	
	const double coef_DBO = (kboff + b13 + mu2);
	const double coef_DBC3 = (kcboff + a11 + mu2);
	const double coef_DBC2 = (kcboff + b11 + a12 + mu2);
	const double coef_DBC1 = (kcboff + b12 + a13 + mu2);
	
	const double coef_D_O = (k_off + b13n + a_22 + ax2 + mu1);
	const double coef_D_C1 = (kc_off + b12 + b_33 + a13n + mu1);
	const double coef_D_C2 = (kc_off + b11 + b_33 + a12 + mu1);
	const double coef_D_C3 = (kc_off + b_33 + a11 + mu1);
	const double coef_D_OS = (bx2 + ki_off);
	const double coef_D_IC3 = (a_33 + a11 + ki_off);
	const double coef_D_IC2 = (a_33 + b11 + a12 + ki_off);
	const double coef_D_IF = (a_33 + a_44 + b_22 + b12 + ki_off);
	const double coef_D_IM1 = (b_44 + a_55);
	const double coef_D_IM2 = b_55;
	
	const double coef_D_BO = (k_off + b13 + mu2);
	const double coef_D_BC3 = (kc_off + a11 + mu2);
	const double coef_D_BC2 = (kc_off + b11 + a12+ mu2);
	const double coef_D_BC1 = (kc_off + b12 + a13 + mu2);
	
	double yo59, yo60, yo61, yo62, yo63, yo64, yo65, yo66, yo67, yo68, yo69, yo70, yo71, yo72, yo73, yo74, yo75, yo76, yo77, yo78, yo79, yo80, yo81, yo82, yo83, yo84, yo85, yo86, yo87, yo88;
	double yo89, yo90, yo91, yo92, yo93, yo94, yo95, yo96, yo97, yo98, yo99, yo100;
	
	double O = y[59];
	double C1 = y[60];
	double C2 = y[61];
	double C3 = y[62];
	double IC3 = y[63];
	double IC2 = y[64];
	double IF = y[65];
	double IM1 = y[66];
	double IM2 = y[67];
	double OS = y[68];
	
    double BO = y[69];
	double BC3 = y[70];
	double BC2 = y[71];
	double BC1 = y[72];
	
    double DO = y[73];
	double DC1 = y[74];
	double DC2 = y[75];
	double DC3 = y[76];
	double DOS = y[77];
	double DIC3 = y[78];
	double DIC2 = y[79];
	double DIF = y[80];
	double DIM1 = y[81];
	double DIM2 = y[82];
	
    double DBO = y[83];
	double DBC3 = y[84];
	double DBC2 = y[85];
	double DBC1 = y[86];
	
    double D_O = y[87];
	double D_C1 = y[88];
	double D_C2 = y[89];
	double D_C3 = y[90];
	double D_OS = y[91];
	double D_IC3 = y[92];
	double D_IC2 = y[93];
	double D_IF = y[94];
	double D_IM1 = y[95];
	double D_IM2 = y[96];
	
    double D_BO = y[97];
	double D_BC3 = y[98];
	double D_BC2 = y[99];
	double D_BC1 = y[100];
	
    double dt = 1.0 * dt0;
    
	double dtinv = 1. / dt;
	double myeps = 1.e-100;
	double err1 = 1.;
	int iter1 = 0;
	
	while ( err1 > myeps && iter1 < 1000 ) {
		yo59 = O;
		yo60 = C1;
		yo61 = C2;
		yo62 = C3;
		yo63 = IC3;
		yo64 = IC2;
		yo65 = IF;
		yo66 = IM1;
		yo67 = IM2;
		yo68 = OS;
		yo69 = BO;
		yo70 = BC3;
		yo71 = BC2;
		yo72 = BC1;
		yo73 = DO;
		yo74 = DC1;
		yo75 = DC2;
		yo76 = DC3;
		yo77 = DOS;
		yo78 = DIC3;
		yo79 = DIC2;
		yo80 = DIF;
		yo81 = DIM1;
		yo82 = DIM2;
		yo83 = DBO;
		yo84 = DBC3;
		yo85 = DBC2;
		yo86 = DBC1;
		yo87 = D_O;
		yo88 = D_C1;
		yo89 = D_C2;
		yo90 = D_C3;
		yo91 = D_OS;
		yo92 = D_IC3;
		yo93 = D_IC2;
		yo94 = D_IF;
		yo95 = D_IM1;
		yo96 = D_IM2;
		yo97 = D_BO;
		yo98 = D_BC3;
		yo99 = D_BC2;
		yo100 = D_BC1;
		
		//Drug Free States
		O = ( y[59] + dt * ( a13 * C1 + b2 * IF + koff * DO + k_off * D_O + bx * OS + mu2 * BO) ) / ( 1. + dt * coef_O );		// O
		C1 = ( y[60] + dt * (a12 * C2 + a3 * IF + b13 * O  + kcoff * DC1 + kc_off * D_C1 + mu2 * BC1 ) ) / (1. + dt * coef_C1 );	// C1
		C2 = ( y[61] + dt * (a11 * C3 + a3 * IC2 + b12 * C1  + kcoff * DC2 + kc_off * D_C2 + mu2 * BC2) ) / ( 1. + dt * coef_C2 ); //C2
		C3 = ( y[62] + dt * (a3 * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 + mu2 * BC3 ) )/( 1.+ dt * coef_C3 ); //C3
		IC3 = ( y[63] + dt * (b3 * C3 + b11 * IC2 + ki_off * D_IC3 ) ) / ( 1. + dt * coef_IC3 );
		IC2 = ( y[64] + dt * (a11 * IC3 + b3 * C2 + b12 * IF + ki_off * D_IC2 ) ) / ( 1. + dt * coef_IC2 );
		IF = ( y[65] + dt * (a12 * IC2 + b3 * C1 + b4 * IM1 + a2 * O + ki_off * D_IF) ) / (1. + dt * coef_IF );
		IM1 = ( y[66] + dt * (a4 * IF + b5 * IM2) ) / (1. + dt * coef_IM1 );
		IM2 = ( y[67] + dt * (a5 * IM1 ) ) / ( 1. + dt * coef_IM2 );
		OS = ( y[68] + dt * (ax * O + ki_off * D_OS) ) / (1. + dt * coef_OS );
		
        BO = ( y[69] + dt * (mu1 * O + a13 * BC1 + kboff* DBO + k_off * D_BO) ) / (1. + dt * coef_BO );
		BC3 = ( y[70] + dt * (mu1 * C3 + b11 * BC2 + kcboff * DBC3 + kc_off * D_BC3) ) / (1. + dt * coef_BC3 );
		BC2 = ( y[71] + dt * (mu1 * C2 + a11 * BC3 + b12 * BC1 + kcboff * DBC2 + kc_off * D_BC2) ) / (1. + dt * coef_BC2 );
		BC1 = ( y[72] + dt * (mu1 * C1 + a12 * BC2 + b13 * BO + kcboff * DBC1 + kc_off * D_BC1) ) / (1. + dt * coef_BC1 );
		
		
		//Charged Drug Bound States
		DO = ( y[73] + dt * (kon * O + a13c * DC1 + bx1 * DOS + b22 * DIF + mu2 * DBO ) ) / ( 1. + dt * coef_DO );
		DC1 = ( y[74] + dt * (kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO  + mu2 * DBC1 ) ) /( 1.+ dt * coef_DC1 );
		DC2 = ( y[75] + dt * (kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1  + mu2 * DBC2 )) / ( 1. + dt * coef_DC2 );
		DC3 = ( y[76] + dt * (kcon * C3 + b11 * DC2 + a33 * DIC3  + mu2 * DBC3 ) ) / ( 1. + dt * coef_DC3 );
		DOS = ( y[77] + dt * (ax1 * DO )) / (1. + dt * coef_DOS );
		DIC3 = ( y[78] + dt * (b33 * DC3 + b11 * DIC2 ) ) / ( 1. + dt * coef_DIC3 );
		DIC2 = ( y[79] + dt * (b33 * DC2 + a11 * DIC3 + b12 * DIF )) /(1. +dt * coef_DIC2 );
		DIF = ( y[80] + dt * (b33 * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO ))/( 1. + dt * coef_DIF );
		DIM1 = ( y[81] + dt * (a44 * DIF + b55 * DIM2 )) / ( 1. + dt * coef_DIM1 );
		DIM2 = ( y[82] + dt * (a55 * DIM1 )) /( 1. + dt * coef_DIM2 );
		
        DBO = ( y[83] + dt * (kbon * BO + a13 * DBC1 + mu1 * DO )) / ( 1. + dt * coef_DBO );
		DBC3 = ( y[84] + dt * (kcbon * BC3 + b11 * DBC2 + mu1 * DC3 )) / (1. + dt * coef_DBC3 );
		DBC2 = ( y[85] + dt * (kcbon * BC2 + a11 * DBC3 + b12 * DBC1 + mu1 * DC2 )) / ( 1. + dt* coef_DBC2 );
		DBC1 = ( y[86] + dt * (kcbon * BC1 + a12 * DBC2 + b13 * DBO + mu1 * DC1 )) / (1. + dt* coef_DBC1 );
		
		//Need to check these
		//Neutral Drug Bound States
		D_O = ( y[87] + dt * (k_on * O + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS + mu2 * D_BO )) / ( 1. + dt * coef_D_O );
		D_C1 = ( y[88] + dt * (kc_on * C1 + a12 * D_C2 + a_33 * D_IF + b13n * D_O + mu2 * D_BC1  )) / ( 1. + dt * coef_D_C1 );
		D_C2 = ( y[89] + dt * (kc_on * C2 + a11 * D_C3 + a_33 * D_IC2 + b12 * D_C1 + mu2 * D_BC2 )) / ( 1. + dt * coef_D_C2 );
		D_C3 = ( y[90] + dt * (kc_on * C3 + a_33 * D_IC3 + b11 * D_C2 + mu2 * D_BC3 )) / ( 1. + dt * coef_D_C3 );
		D_OS = ( y[91] + dt * (ax2 * D_O + ki_on * OS )) / ( 1. + dt * coef_D_OS );
		D_IC3 = ( y[92] + dt * (b_33 * D_C3 + b11 * D_IC2 + ki_on * IC3 )) / ( 1. + dt * coef_D_IC3 );
		D_IC2 = ( y[93] + dt * (b_33 * D_C2 + a11 * D_IC3 + b12 * D_IF + ki_on * IC2 )) / ( 1. + dt * coef_D_IC2 );
		D_IF = ( y[94] + dt * (b_33 * D_C1 + a12 * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF )) / ( 1. + dt * coef_D_IF );
		D_IM1 = ( y[95] + dt * (a_44 * D_IF + b_55 * D_IM2 )) / ( 1. + dt * coef_D_IM1 );
		D_IM2 = ( y[96] + dt * (a_55 * D_IM1 )) / ( 1. + dt * coef_D_IM2 );
		
        D_BO = ( y[97] + dt * (k_on * BO + a13 * D_BC1 + mu1 * D_O )) / ( 1. + dt * coef_D_BO  );
		D_BC3 = ( y[98] + dt * (kc_on * BC3 + b11 * D_BC2 + mu1 * D_C3 )) / ( 1. + dt * coef_D_BC3  );
		D_BC2 = ( y[99] + dt * (kc_on * BC2 + a11 * D_BC3 + b12 * D_BC1 + mu1 * D_C2 )) / ( 1. + dt * coef_D_BC2 );
		D_BC1 = ( y[100] + dt * (kc_on * BC1 + a12 * D_BC2 + b13 * D_BO + mu1 * D_C1 )) / ( 1. + dt* coef_D_BC1 );
        
		err1 = ( fabs( yo59 - O )
				+ fabs( yo60 - C1 )
				+ fabs( yo61 - C2 )
				+ fabs( yo62 - C3 )
				+ fabs( yo63 - IC3 )
				+ fabs( yo64 - IC2 )
				+ fabs( yo65 - IF )
				+ fabs( yo66 - IM1 )
				+ fabs( yo67 - IM2 )
				+ fabs( yo68 - OS )
				+ fabs( yo69 - BO )
				+ fabs( yo70 - BC3 )
				+ fabs( yo71 - BC2 )
				+ fabs( yo72 - BC1 )
				+ fabs( yo73 - DO )
				+ fabs( yo74 - DC1 )
				+ fabs( yo75 - DC2 )
				+ fabs( yo76 - DC3 )
				+ fabs( yo77 - DOS )
				+ fabs( yo78 - DIC3 )
				+ fabs( yo79 - DIC2 )
				+ fabs( yo80 - DIF )
				+ fabs( yo81 - DIM1 )
				+ fabs( yo82 - DIM2 )
				+ fabs( yo83 - DBO )
				+ fabs( yo84 - DBC3 )
				+ fabs( yo85 - DBC2 )
				+ fabs( yo86 - DBC1 )
				+ fabs( yo87 - D_O )
				+ fabs( yo88 - D_C1 )
				+ fabs( yo89 - D_C2 )
				+ fabs( yo90 - D_C3 )
				+ fabs( yo91 - D_OS )
				+ fabs( yo92 - D_IC3 )
				+ fabs( yo93 - D_IC2 )
				+ fabs( yo94 - D_IF )
				+ fabs( yo95 - D_IM1 )
				+ fabs( yo96 - D_IM2 )
				+ fabs( yo97 - D_BO )
				+ fabs( yo98 - D_BC3 )
				+ fabs( yo99 - D_BC2 )
				+ fabs( yo100 - D_BC1 )
				
				);
		iter1++;
	}
	
	ydot[59] = ( O - y[59] ) * dtinv;//O
	ydot[60] = (C1 - y[60] ) * dtinv;//C1
	ydot[61] = ( C2 - y[61] ) * dtinv;//C2
	ydot[62] = ( C3 - y[62] ) * dtinv;//C3
	ydot[63] = ( IC3 - y[63] ) * dtinv;//IC3
	ydot[64] = ( IC2 - y[64] ) * dtinv;//IC2
	ydot[65] = ( IF - y[65] ) * dtinv;//IF
	ydot[66] = ( IM1 - y[66] ) * dtinv;//IM1
	ydot[67] = ( IM2 - y[67] ) * dtinv;//IM2
	ydot[68] = ( OS - y[68] ) * dtinv;//OS
    
	ydot[69] = ( BO - y[69] ) * dtinv;//BO
	ydot[70] = ( BC3 - y[70] ) * dtinv;//BC3
	ydot[71] = ( BC2 - y[71] ) * dtinv;//BC2
	ydot[72] = ( BC1 - y[72] ) * dtinv;//BC1
    
	ydot[73] = ( DO - y[73] ) * dtinv;//DO
	ydot[74] = ( DC1 - y[74] ) * dtinv;//DC1
	ydot[75] = ( DC2 - y[75] ) * dtinv;//DC2
	ydot[76] = ( DC3 - y[76] ) * dtinv;//DC3
	ydot[77] = ( DOS - y[77] ) * dtinv;//DOS
	ydot[78] = ( DIC3 - y[78] ) * dtinv;//DIC3
	ydot[79] = ( DIC2 - y[79] ) * dtinv;//DIC2
	ydot[80] = ( DIF - y[80] ) * dtinv;//DIF
	ydot[81] = ( DIM1 - y[81] ) * dtinv;//DIM1
	ydot[82] = ( DIM2 - y[82] ) * dtinv;//DIM2
    
	ydot[83] = ( DBO - y[83] ) * dtinv;//DBO
	ydot[84] = ( DBC3 - y[84] ) * dtinv;//DBC3
	ydot[85] = ( DBC2 - y[85] ) * dtinv;//DBC2
	ydot[86] = ( DBC1- y[86] ) * dtinv;//DBC1
    
	ydot[87] = ( D_O - y[87] ) * dtinv;//D_O
	ydot[88] = ( D_C1 - y[88] ) * dtinv;//D_C1
	ydot[89] = ( D_C2 - y[89] ) * dtinv;//D_C2
	ydot[90] = ( D_C3 - y[90] ) * dtinv;//D_C3
	ydot[91] = ( D_OS - y[91] ) * dtinv;//D_OS
	ydot[92] = ( D_IC3 - y[92] ) * dtinv;//D_IC3
	ydot[93] = ( D_IC2 - y[93] ) * dtinv;//D_IC2
	ydot[94] = ( D_IF - y[94] ) * dtinv;//D_IF
	ydot[95] = ( D_IM1 - y[95] ) * dtinv;//D_IM1
	ydot[96] = ( D_IM2 - y[96] ) * dtinv;//D_IM2
    
	ydot[97] = ( D_BO - y[97] ) * dtinv;//D_BO
	ydot[98] = ( D_BC3 - y[98] ) * dtinv;//D_BC3
	ydot[99] = ( D_BC2 - y[99] ) * dtinv;//D_BC2
	ydot[100] = ( D_BC1 - y[100] ) * dtinv;//D_BC1

	
	double I_Na_junc = Fjunc*GNa*(y[59]+y[69])*(y[38]-ena_junc);
	double I_Na_sl = Fsl*GNa*(y[59]+y[69])*(y[38]-ena_sl);
	double I_Na = I_Na_junc + I_Na_sl;

    // end
    
    //IK-Ach
    
    const double gkach2 = 1. /(1 + pow((0.03/Ach),2.1) );
    double I_kach = gkach2 * (0.08 + ( 0.04/(1 + exp( (y[38]+91)/12 ) ) ) ) * (y[38]-ek);
    
    if (Ach == 0.) { I_kach = 0.0; }
    else { I_kach = I_kach; }
    
    //  I_nabk: Na Background Current
    const double I_nabk_junc = Fjunc * GNaB * ( y[38] - ena_junc );
    const double I_nabk_sl = Fsl * GNaB * ( y[38] - ena_sl );
    const double I_nabk = I_nabk_junc + I_nabk_sl;
    
    //  I_nak: Na/K Pump Current
    const double sigma = ( exp( Nao / 67.3 ) - 1. ) / 7.;
    const double fnak = 1. / ( 1. + 0.1245 * exp( -0.1 * y[38] * FoRT ) + 0.0365 * sigma * exp(- y[38] * FoRT ) );
    const double I_nak_junc = 1. * Fjunc * IbarNaK * fnak * Ko / ( 1. + pow( ( KmNaip / y[31] ) , 4 ) ) / ( Ko + KmKo );
    const double I_nak_sl = 1. * Fsl * IbarNaK * fnak * Ko / ( 1. + pow( ( KmNaip / y[32] ) , 4 ) ) / ( Ko + KmKo );
    const double I_nak = I_nak_junc + I_nak_sl;
    
    // //  I_kr: Rapidly Activating K Current
    double I_kr;
    const double markov_ikr = 1;
   
    double IC50 = 1.5 * (1e-6);
    double factor_flec = 1/(1+(drug_onKr/IC50));

   
    // Markov IKr
    
    const double GKr = 0.024 * 1.3 * factor_flec ;
    
    double ae=exp(24.335+0.0112*y[38]-25.914);
    double be=exp(13.688-0.0603*y[38]-15.707);
    double ai=exp(30.061-0.0312*y[38]-33.243);
    double ain=exp(22.746-25.914);
    double bin=exp(13.193-15.707);
    double bi=exp(30.016+0.0223*y[38]-30.888)*pow((5.4/Ko),0.4);   //inactivation
    double aa=exp(22.098+0.0365*y[38]-25.914);  //activation
    double bb=exp(7.313-0.0399*y[38]-15.707);  //deactivation
    
 
        
//        double dO = (ai*I + aa*C1  - O*(bi + bb));
        ydot[102] = (ai * y[106] + aa * y[103]  - y[102] * (bi + bb ));
        
//        double dC1 = (ain*C2 + bb*O - C1*(aa + bin));
        ydot[103] = (ain * y[104] + bb * y[102] - y[103] * (aa + bin));
        
//        double dC2 = (ae*C3 + bin*C1 - C2*(be + ain));
        ydot[104] = (ae * y[105]+ bin * y[103] - y[104] * (be + ain));
        
//        double dC3 = (be*C2 - ae*C3);
        ydot[105] = (be * y[104] - ae* y[105]);
        
//        double dI = (bi*O  - I*ai);
        ydot[106] = (bi * y[102]  - y[106] * ai );
    
        ydot[107] = 0;

        ydot[108] = 0;

        ydot[109] = 0;

        ydot[110] = 0;

       ydot[111] = 0;
        
        
        
        I_kr = GKr * sqrt(Ko/5.4) * y[102] * (y[38] - ek);
    
    
    
    // //  I_ks: Slowly Activating K Current
    const double markov_iks = 0;
   
    
    const double eks = (1/FoRT)*log((Ko+pNaK*Nao)/( y[34] +pNaK* y[33] ));
    
    double  gks_junc, gks_sl, xsss,  tauxs ,I_ks_junc , I_ks_sl , I_ks , alpha, beta, gamma, delta, teta, eta, psi, omega, O2;

    if( markov_iks == 0 ) {
         gks_junc = 1. * ( 1. + 1. * AF + 2. * ISO ) * 0.0035  ;
         gks_sl = 1. * ( 1. + 1. * AF + 2. * ISO ) * 0.0035  ; // FRA
         xsss = 1. / ( 1. + exp( -( y[38] + 40. * ISO + 3.8 ) / 14.25 ) ); //  fitting Fra
         tauxs = 990.1 / ( 1. + exp( -( y[38] + 40. * ISO + 2.436 ) / 14.12 ) );
        ydot[12]  = ( xsss - y[12] ) / tauxs;
         I_ks_junc = Fjunc * gks_junc * y[12] * y[12] * ( y[38] - eks );
         I_ks_sl = Fsl * gks_sl * y[12] * y[12] * ( y[38] - eks );
         I_ks = I_ks_junc + I_ks_sl;
    } else {
         gks_junc = 1. * 0.0065;
         gks_sl = 1. * 0.0065; // FRA
         alpha = 3.98e-4 * exp( 3.61e-1 * y[38] * FoRT );
         beta = 5.74e-5 * exp( -9.23e-2 * y[38] * FoRT );
         gamma = 3.41e-3 * exp( 8.68e-1 * y[38] * FoRT );
         delta = 1.2e-3 * exp( -3.3e-1 * y[38] * FoRT );
         teta = 6.47e-3;
         eta = 1.25e-2 * exp( -4.81e-1 * y[38] * FoRT );
         psi = 6.33e-3 * exp( 1.27 * y[38] * FoRT );
         omega = 4.91e-3 * exp( -6.79e-1 * y[38] * FoRT );
        
        ydot[41] = -4. * alpha * y[41] + beta * y[42] ;
        ydot[42] = 4. * alpha * y[41] - ( beta + gamma + 3. * alpha ) * y[42] + 2. * beta * y[43] ;
        ydot[43] = 3. * alpha * y[42] - ( 2. * beta + 2. * gamma + 2. * alpha ) * y[43] + 3. * beta * y[44] ;
        ydot[44] = 2. * alpha * y[43] - ( 3. * beta + 3. * gamma + alpha ) * y[44] + 4. * beta * y[45] ;
        ydot[45] = 1. * alpha * y[43] - ( 4. * beta + 4. * gamma ) * y[45] + delta * y[49] ;
        ydot[46] = gamma * y[42] - ( delta + 3. * alpha ) * y[46] + beta * y[47] ;
        ydot[47] = 2. * gamma * y[43] + 3. * alpha * y[46] - ( delta + beta + 2. * alpha + gamma ) * y[47] + 2. * beta * y[48] + 2. * delta * y[50] ;
        ydot[48] = 3. * gamma * y[44] + 2. * alpha * y[47] - ( delta + 2. * beta + 1. * alpha + 2. * gamma ) * y[48] + 3. * beta * y[49] + 2. * delta * y[51] ;
        ydot[49] = 4. * gamma * y[45] + 1. * alpha * y[48] - ( delta + 3. * beta + 0 * alpha + 3. * gamma ) * y[49] + 2. * delta * y[52] ;
        ydot[50] = 1. * gamma * y[47] - ( 2. * delta + 2. * alpha ) * y[50] + beta * y[51] ;
        ydot[51] = 2. * gamma * y[48] + 2. * alpha * y[50] - ( 2. * delta + beta + 1. * alpha + gamma ) * y[51] + 2. * beta * y[52] + 3. * delta * y[53] ;
        ydot[52] = 3. * gamma * y[49] + 1. * alpha * y[51] - ( 2. * delta + 2. * beta + 2. * gamma ) * y[52] + 3. * delta * y[54] ;
        ydot[53] = 1. * gamma * y[51] - ( 3. * delta + 1. * alpha ) * y[53] + beta * y[54] ;
        ydot[54] = 2. * gamma * y[52] + 1. * alpha * y[53] - ( 3. * delta + 1. * beta + 1. * gamma ) * y[54] + 4. * delta * y[55] ;
        ydot[55] = 1. * gamma * y[54] - ( 4. * delta + teta ) * y[55] + eta * y[56] ;
         O2 = 1. - ( y[41] + y[42] + y[43] + y[44] + y[45] + y[46] + y[48] + y[47] + y[49] + y[50] + y[51] + y[52] + y[53] + y[54] + y[55] + y[56] );
        ydot[56] = 1. * teta * y[55] - ( eta + psi ) * y[56] + omega * O2;
         I_ks_junc = Fjunc * gks_junc * ( y[56] +O2 ) * ( y[38] - eks );
         I_ks_sl = Fsl * gks_sl * ( y[56] + O2 ) * ( y[38] - eks );
         I_ks = I_ks_junc + I_ks_sl;
    }
    // end
    
    // I_kp: Plateau K current
    const double kp_kp = 1. / ( 1. + exp( 7.488 - y[38] / 5.98 ) );
    const double I_kp_junc = Fjunc * gkp * kp_kp * ( y[38] - ek );
    const double I_kp_sl = Fsl * gkp * kp_kp * ( y[38] - ek );
    const double I_kp = I_kp_junc + I_kp_sl;
    
    // //  I_to: Transient Outward K Current (slow and fast components)
    //  modified for human myocytes
    double IC50_FlecIto = 15.2 * (1e-6);
    double factor_FlecIto = 1/(1+(drug_onIto/IC50_FlecIto));
    
    const double GtoFast = ( 1.0 - 0.7 * AF ) * 0.165 * factor_FlecIto; // nS/pF maleckar; // human atrium
    
    // 11/12/09; changed Itof to that from maleckar/giles/2009; removed I_tos
    // atrium
    // equations for activation;
    const double xtoss = ( (1.) / ( 1. + exp( -( y[38] + 1.0 ) / 11.0 ) ) );
    const double tauxtof = 3.5 * exp( -( pow( ( y[38] / 30.0 ) , 2.0 ) ) ) + 1.5;
    
    // equations for inactivation;
    const double ytoss = ( (1.0) / ( 1. + exp( ( y[38] + 40.5 ) / 11.5 ) ) ) ;
    const double tauytof = 25.635 * exp( -( pow( ( ( y[38] + 52.45 ) / 15.8827 ) , 2.0 ) ) ) + 24.14;// 14.14
    
    ydot[9]  = ( xtoss- y[9] ) / tauxtof;
    ydot[10]  = ( ytoss - y[10] ) / tauytof;
    const double I_tof = 1.0 * GtoFast * y[9] * y[10] * ( y[38] - ek );
    
    const double I_to = 1. * I_tof;
    
    // //  I_kur: Ultra rapid delayed rectifier Outward K Current
    // Equation for IKur; from Maleckar et al. 2009 - EG
    // atrium
    // equations for activation;
    const double Gkur = 1. * ( 1.0 - 0.5 * AF ) * ( 1. + 2. * ISO ) * 0.045 * ( 1. + 0.2 * RA ); // nS/pF maleckar 0.045
    const double xkurss = ( (1.) / ( 1 + exp( ( y[38] + 6. ) / (-8.6) ) ) );
    const double tauxkur = 9. / ( 1. + exp( ( y[38] + 5. ) / 12.0 ) ) + 0.5;
    
    // equations for inactivation;
    const double ykurss = ( (1.) / ( 1. + exp( ( y[38] + 7.5 ) / 10.0 ) ) );
    const double tauykur = 590. / ( 1. + exp( ( y[38] + 60. ) / 10.0 ) ) + 3050.;
    
    ydot[57]  = ( xkurss - y[57] ) / tauxkur;
    ydot[58]  = ( ykurss - y[58] ) / tauykur;
    const double I_kur = 1. * Gkur * y[57] * y[58] * ( y[38] - ek );
    
    // //  I_ki: Time-Independent K Current
    const double aki = 1.02 / ( 1. + exp( 0.2385 * ( y[38] - ek - 59.215 ) ) );
    const double bki = ( 0.49124 * exp( 0.08032 * ( y[38] + 5.476 - ek ) ) + exp( 0.06175 * ( y[38] - ek - 594.31 ) ) ) / ( 1. + exp( -0.5143 * ( y[38] - ek + 4.753 ) ) );
    const double kiss = aki / ( aki + bki );
    
    // I_ki =1* 0.35*sqrt(Ko/5.4)*kiss*( y[38] -ek);
    // SVP 11/11/09
    // multiplieD IK1 by 0.15 to scale it to single cell isolated atrial cell
    // resting potential
    const double I_ki = ( 1. + 1. * AF ) * 0.0525 * sqrt( Ko / 5.4 ) * kiss * ( y[38] - ek );
    
    //  I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
    const double I_ClCa_junc = Fjunc * GClCa / ( 1. + KdClCa / y[35] ) * ( y[38] - ecl );
    const double I_ClCa_sl = Fsl * GClCa / ( 1. + KdClCa / y[36] ) * ( y[38] - ecl );
    const double I_ClCa = I_ClCa_junc + I_ClCa_sl;
    const double I_Clbk = GClB * ( y[38] - ecl );
    
    const double GClCFTR = 0;// 4.9e-3*ISO;     //  [mS/uF]
    const double I_ClCFTR = GClCFTR * ( y[38] - ecl );
    
    // //  I_Ca: L-type Calcium Current
    const double dss = 1. / ( 1. + exp( -( y[38] + 3 * ISO + 9 ) / 6 ) ); // in Maleckar v1/2=-9 S=6 (mV); Courtemanche v1/2=-9 S=5.8 (mV)
    const double taud = 1. * dss * ( 1. - exp( -( y[38] + 3 * ISO + 9 ) / 6 ) ) / ( 0.035 * ( y[38] + 3 * ISO + 9 ) );
    const double fss = 1. / ( 1. + exp( ( y[38] + 3. * ISO + 30. ) / 7. ) ) + 0.2 / ( 1. + exp( ( 50. - y[38] - 3. * ISO ) / 20. ) ); //  in Maleckar v1/2=-27.4 S=7.1 (mV); Courtemanche v1/2=-28 S=6.9 (mV)
    const double tauf = 1. / ( 0.0197 * exp( - pow( ( 0.0337 * ( y[38] + 3. * ISO + 25. ) ) , 2 ) ) + 0.02 );
    ydot[3]  = (dss- y[3] )/taud;
    ydot[4]  = (fss- y[4] )/tauf;
    ydot[5]  = 1.7 * y[35] * ( 1. - y[5] ) - 1. * 11.9e-3 * y[5] ; //  fCa_junc   koff!!!!!!!!
    ydot[6]  = 1.7 * y[36] * ( 1. - y[6] ) - 1. * 11.9e-3 * y[6] ; //  fCa_sl
    // const double fcaCaMSL = 0.1 / ( 1. + ( 0.01 / y[36] ) );
    // const double fcaCaj = 0.1 / ( 1. + ( 0.01 / y[35] ) );
    const double fcaCaMSL = 0;
    const double fcaCaj = 0;
    const double ibarca_j = pCa * 4. * ( y[38] * Frdy * FoRT ) * ( 0.341 * y[35] * exp( 2 * y[38] * FoRT ) - 0.341 * Cao ) / ( exp( 2 * y[38] * FoRT ) - 1. );
    const double ibarca_sl = pCa * 4. * ( y[38] * Frdy * FoRT ) * ( 0.341 * y[36] * exp( 2 * y[38] * FoRT ) - 0.341 * Cao ) / ( exp( 2 * y[38] * FoRT ) - 1. );
    const double ibark = pK * ( y[38] * Frdy * FoRT ) * ( 0.75 * y[34] * exp( y[38] * FoRT ) - 0.75 * Ko ) / ( exp( y[38] * FoRT ) - 1. );
    const double ibarna_j = pNa * ( y[38] * Frdy * FoRT ) * ( 0.75 * y[31] * exp( y[38] * FoRT ) - 0.75 * Nao ) / ( exp( y[38] * FoRT ) - 1. );
    const double ibarna_sl = pNa * ( y[38] * Frdy * FoRT ) * ( 0.75 * y[32] * exp( y[38] * FoRT ) - 0.75 * Nao ) / ( exp( y[38] * FoRT ) - 1. );
    const double I_Ca_junc = ( Fjunc_CaL * ibarca_j * y[3] * y[4] *( ( 1. - y[5] ) + fcaCaj ) * pow( Q10CaL , Qpow ) ) * 0.45;
    const double I_Ca_sl = ( Fsl_CaL * ibarca_sl * y[3] * y[4] * ( ( 1. - y[6] ) + fcaCaMSL ) * pow( Q10CaL , Qpow ) ) * 0.45;
    const double I_Ca = I_Ca_junc + I_Ca_sl;
    const double I_CaK = ( ibark * y[3] * y[4] * ( Fjunc_CaL * ( fcaCaj + ( 1. - y[5] ) ) + Fsl_CaL * ( fcaCaMSL + ( 1. - y[6] ) ) ) * pow( Q10CaL , Qpow ) ) * 0.45;
    const double I_CaNa_junc = ( Fjunc_CaL * ibarna_j * y[3] * y[4] * ( ( 1. - y[5] ) + fcaCaj ) * pow( Q10CaL , Qpow ) ) * 0.45;
    const double I_CaNa_sl = ( Fsl_CaL * ibarna_sl * y[3] * y[4] * ( ( 1. - y[6] ) + fcaCaMSL ) * pow( Q10CaL , Qpow ) ) * 0.45;
    const double I_CaNa = I_CaNa_junc + I_CaNa_sl;
    const double I_Catot = I_Ca + I_CaK + I_CaNa;
    
    //  I_ncx: Na/Ca Exchanger flux
    const double Ka_junc = 1. / ( 1. + ( Kdact * Kdact / ( y[35] * y[35] ) ) );
    const double Ka_sl = 1. / ( 1. + ( Kdact * Kdact / ( y[36] * y[36] ) ) );
    const double s1_junc = exp( nu * y[38] * FoRT ) * y[31] * y[31] * y[31] * Cao;
    const double s1_sl = exp( nu * y[38] * FoRT ) * y[32] * y[32] * y[32] * Cao;
    const double s2_junc = exp( ( nu - 1. ) * y[38] * FoRT ) * Nao * Nao * Nao * y[35] ;
    const double s3_junc = ( KmCai * Nao * Nao * Nao * ( 1. + ( y[31] * y[31] * y[31] / ( KmNai * KmNai * KmNai ) ) )
                            + KmNao * KmNao * KmNao * y[35] * ( 1. + y[35] / KmCai )
                            + KmCao * y[31] * y[31] * y[31] + y[31] * y[31] * y[31] * Cao
                            + Nao * Nao * Nao * y[35] );
    const double s2_sl = exp( ( nu - 1. ) * y[38] * FoRT ) * Nao * Nao * Nao * y[36] ;
    // const double s3_sl = KmCai * Nao * Nao * Nao * ( 1. + ( y[32] / KmNai ) ^3) + KmNao ^3 * y[36] * ( 1. + y[36] / KmCai ) + KmCao * y[32] ^3+ y[32] ^3 * Cao + Nao ^3 * y[36] ;
    const double s3_sl = ( KmCai * Nao * Nao * Nao * ( 1. + ( y[32] * y[32] * y[32] / ( KmNai * KmNai * KmNai ) ) )
                          + KmNao * KmNao * KmNao * y[36] * ( 1. + y[36] / KmCai )
                          + KmCao * y[32] * y[32] * y[32]
                          + y[32] * y[32] * y[32] * Cao
                          + Nao * Nao * Nao * y[36] );
    
    
    const double I_ncx_junc = Fjunc * IbarNCX * pow( Q10NCX , Qpow ) * Ka_junc * ( s1_junc - s2_junc ) / s3_junc / ( 1. + ksat * exp( ( nu - 1. ) * y[38] * FoRT ) );
    const double I_ncx_sl = Fsl * IbarNCX * pow( Q10NCX , Qpow ) * Ka_sl * ( s1_sl - s2_sl ) / s3_sl / ( 1. + ksat * exp( ( nu - 1. ) * y[38] * FoRT ) );
    const double I_ncx = I_ncx_junc + I_ncx_sl;
    
    //  I_pca: Sarcolemmal Ca Pump Current
    const double I_pca_junc = Fjunc * pow( Q10SLCaP , Qpow ) * IbarSLCaP * pow( y[35] , 1.6 ) / ( pow( KmPCa , 1.6 ) + pow( y[35] , 1.6 ) );
    const double I_pca_sl = Fsl * pow( Q10SLCaP , Qpow ) * IbarSLCaP * pow( y[36] , 1.6 ) / ( pow( KmPCa , 1.6 ) + pow( y[36] , 1.6 ) );
    const double I_pca = I_pca_junc + I_pca_sl;
    
    //  I_cabk: Ca Background Current
    const double I_cabk_junc = Fjunc * GCaB * ( y[38] - eca_junc );
    const double I_cabk_sl = Fsl * GCaB * ( y[38] - eca_sl);
    const double I_cabk = I_cabk_junc + I_cabk_sl;
    
    // //  SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
    
    double kon_flec = drug_onRyR * 5500;
    
    double IC50_flecRyR = 20e-6; //test 5e-6, 11e-6, 17e-6 27.1e-6??
    double koff_flecRyR = IC50_flecRyR * 5500;
    
    const double MaxSR = 15.;
    const double MinSR = 1.;
    const double kCaSR = MaxSR - ( MaxSR - MinSR ) / ( 1. + pow( ( ec50SR / y[30] ) , 2.5 ) );
    const double koSRCa = (1.) * koCa / kCaSR;//
    const double kiSRCa = kiCa * kCaSR;
    const double RI = 1. - y[13] - y[14] - y[15] - y[101];
    ydot[13]  = ( kim * RI - kiSRCa * y[35] * y[13] ) - ( koSRCa * y[35] * y[35] * y[13] - kom * y[14] );   //  R
    
    ydot[14]  = ( (koSRCa * y[35] * y[35] * y[13] - kom * y[14] ) - ( kiSRCa * y[35] * y[14] - kim * y[15] )+koff_flecRyR*y[101]-kon_flec*y[14] );//  O
    
    ydot[15]  = ( kiSRCa * y[35] * y[14] - kim * y[15] ) - ( kom * y[15] - koSRCa * y[35] * y[35] * RI);   //  I
    
    ydot[101] = ((kon_flec * y[14] - koff_flecRyR * y[101]));
    
    const double J_SRCarel = ks * (y[14]+0.2*y[101])  * ( y[30] - y[35] );          //  [mM/ms]
    
    const double J_serca = ( 1.0 * pow( Q10SRCaP, Qpow ) * Vmax_SRCaP * ( pow( ( y[37] / Kmf ) , hillSRCaP ) - pow( ( y[30] / Kmr ) , hillSRCaP ) )
                            / ( 1. + pow( ( y[37] / Kmf ) , hillSRCaP ) + pow( ( y[30] / Kmr ) , hillSRCaP ) ) );
    const double J_SRleak = (1.) * ( 1.0 + 0.25 * AF ) * 5.348e-6 * ( y[30] - y[35] );           //    [mM/ms]
    
    
    // //  Sodium and Calcium Buffering
    ydot[16]  = kon_na* y[31] *(Bmax_Naj- y[16] )-koff_na* y[16] ;        //  NaBj      [mM/ms]
    ydot[17]  = kon_na* y[32] *(Bmax_Nasl- y[17] )-koff_na* y[17] ;       //  NaBsl     [mM/ms]
    
    //  Cytosolic Ca Buffers
    ydot[18]  = kon_tncl* y[37] *(Bmax_TnClow- y[18] )-koff_tncl* y[18] ;            //  TnCL      [mM/ms]
    ydot[19]  = kon_tnchca* y[37] *(Bmax_TnChigh- y[19] - y[20] )-koff_tnchca* y[19] ; //  TnCHc     [mM/ms]
    ydot[20]  = kon_tnchmg*Mgi*(Bmax_TnChigh- y[19] - y[20] )-koff_tnchmg* y[20] ;   //  TnCHm     [mM/ms]
    ydot[21]  = kon_cam* y[37] *(Bmax_CaM- y[21] )-koff_cam* y[21] ;                 //  CaM       [mM/ms]
    ydot[22]  = kon_myoca* y[37] *(Bmax_myosin- y[22] - y[23] )-koff_myoca* y[22] ;    //  Myosin_ca [mM/ms]
    ydot[23]  = kon_myomg*Mgi*(Bmax_myosin- y[22] - y[23] )-koff_myomg* y[23] ;      //  Myosin_mg [mM/ms]
    ydot[24]  = kon_sr* y[37] *(Bmax_SR- y[24] )-koff_sr* y[24] ;                    //  SRB       [mM/ms]
    const double J_CaB_cytosol = ydot[18] + ydot[19] + ydot[20] + ydot[21] + ydot[22] + ydot[23] + ydot[24]; // sum(ydot(19:25));
    
    //  Junctional and SL Ca Buffers
    ydot[25]  = kon_sll* y[35] *(Bmax_SLlowj- y[25] )-koff_sll* y[25] ;       //  SLLj      [mM/ms]
    ydot[26]  = kon_sll* y[36] *(Bmax_SLlowsl- y[26] )-koff_sll* y[26] ;      //  SLLsl     [mM/ms]
    ydot[27]  = kon_slh* y[35] *(Bmax_SLhighj- y[27] )-koff_slh* y[27] ;      //  SLHj      [mM/ms]
    ydot[28]  = kon_slh* y[36] *(Bmax_SLhighsl- y[28] )-koff_slh* y[28] ;     //  SLHsl     [mM/ms]
    const double J_CaB_junction =  ydot[25] + ydot[27] ;
    const double J_CaB_sl =  ydot[26] + ydot[28] ;
    
    // //  Ion concentrations
    //  SR Ca Concentrations
    ydot[29]  = kon_csqn * y[30] * ( Bmax_Csqn - y[29] ) - koff_csqn * y[29] ;       //  Csqn      [mM/ms]
    ydot[30]  = J_serca - ( J_SRleak * Vmyo / Vsr + J_SRCarel ) - ydot[29] ;         //  Ca_sr     [mM/ms] // Ratio 3 leak current
    //   ydot[30] =0;
    
    //  Sodium Concentrations
    const double I_Na_tot_junc = I_Na_junc + I_nabk_junc + 3 * I_ncx_junc + 3 * I_nak_junc + I_CaNa_junc ;   //  [uA/uF]
    const double I_Na_tot_sl = I_Na_sl + I_nabk_sl + 3 * I_ncx_sl + 3 * I_nak_sl + I_CaNa_sl ;   //  [uA/uF]
    const double I_Na_tot_sl2 = 3 * I_ncx_sl + 3 * I_nak_sl + I_CaNa_sl;   //  [uA/uF]
    const double I_Na_tot_junc2 = 3 * I_ncx_junc + 3 * I_nak_junc + I_CaNa_junc;   //  [uA/uF]
    
    ydot[31]  = -I_Na_tot_junc * Cmem / ( Vjunc * Frdy ) + J_na_juncsl / Vjunc * ( y[32] - y[31] ) - ydot[16] ;
    ydot[32]  = ( -I_Na_tot_sl * Cmem / ( Vsl * Frdy ) + J_na_juncsl / Vsl * ( y[31] - y[32] )
                 + J_na_slmyo / Vsl * ( y[33] - y[32] ) - ydot[17] );
    // FluxNaSL= ydot[32] ;
    //   ydot[31]  = 0;
    //   ydot[32]  = 0;
    ydot[33]  = J_na_slmyo/Vmyo*( y[32] - y[33] );             //  [mM/msec]
    //   ydot[33] =0;
    
    //  Potassium Concentration
    const double I_K_tot = I_to + I_kr + I_ks + I_ki - 2 * I_nak + I_CaK + I_kp + I_kur + I_kach;     //  [uA/uF] // SVP: added IKur
    //   ydot[34]  = 0; // -I_K_tot*Cmem/(Vmyo*Frdy);           //  [mM/msec]
    ydot[34]  =  -I_K_tot*Cmem/(Vmyo*Frdy);
    
    //  Calcium Concentrations
    const double I_Ca_tot_junc = I_Ca_junc + I_cabk_junc + I_pca_junc - 2 * I_ncx_junc;                   //  [uA/uF]
    const double I_Ca_tot_sl = I_Ca_sl + I_cabk_sl + I_pca_sl - 2 * I_ncx_sl;            //  [uA/uF]
    ydot[35]  = ( -I_Ca_tot_junc * Cmem / ( Vjunc * 2 * Frdy ) + J_ca_juncsl / Vjunc * ( y[36] - y[35] )
                 - J_CaB_junction + ( J_SRCarel ) * Vsr / Vjunc + J_SRleak * Vmyo / Vjunc );  //  Ca_j
    ydot[36]  = ( -I_Ca_tot_sl * Cmem / ( Vsl * 2 * Frdy ) + J_ca_juncsl / Vsl * ( y[35] - y[36] )
                 + J_ca_slmyo / Vsl * ( y[37] - y[36] ) - J_CaB_sl );   //  Ca_sl
    //   ydot[35] =0;
    //   ydot[36] =0;
    //   ydot[37]  = -J_serca*Vsr/Vmyo-J_CaB_cytosol;// +J_ca_slmyo/Vmyo*( y[36] - y[37] );    //  [mM/msec]
    ydot[37]  = -J_serca * Vsr / Vmyo - J_CaB_cytosol + J_ca_slmyo / Vmyo * ( y[36] - y[37] );
    //   ydot[37] =0;
    

    // //  Membrane Potential
    // //
    const double I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          //  [uA/uF]
    const double I_Cl_tot = I_ClCa + I_Clbk + I_ClCFTR;                        //  [uA/uF]
    const double I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
    const double I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
    //  ydot[38]  = -(I_Ca_tot+I_K_tot+I_Na_tot-I_app);
    ydot[38]  = -(I_tot);
  
    const double vmax =  ydot[38] ;
    //  ----- END EC COUPLING MODEL ---------------
    //  adjust output depending on the function call
    
    //    if (nargin == 3)
    //        output = ydot;
    //    elseif (nargin == 4) & strcmp(runType,'ydot')
    //        output = ydot;
    //    elseif (nargin == 4) & strcmp(runType,'rates')
    //        output = r;
    //    elseif (nargin == 4) & strcmp(runType,'currents')
    //        // currents = [I_Na I_nabk I_nak I_kr I_ks I_kp I_tos I_tof I_ki I_ClCa I_Clbk I_Catot I_ncx I_pca I_cabk J_serca*Vmyo/Vsr];
    //        // currents = [I_Na I_tof I_tos I_kr I_ks I_ClCa I_Catot J_SRCarel*Vsr/Vmyo J_SRleak RI I_ncx];
    //        //      currents= [I_Catot J_SRCarel*2*Vsr*Frdy/Cmem];
    //        //  currents=[I_Na I_Catot I_ncx I_nak I_kr I_ks I_tof I_tos I_ki];
    //        currents=[];
    //        output = currents;
    //    end
    currents[0] = I_Na ;
    currents[1] = I_Catot;
    currents[2] = I_kr;
    currents[3] = I_ks;
    currents[4] = I_to;
    currents[5] = I_kur;
    currents[6] = I_ki;
    currents[7] = I_nak;
    currents[8] = I_ncx;
    currents[9] = I_kach;

}

#endif
