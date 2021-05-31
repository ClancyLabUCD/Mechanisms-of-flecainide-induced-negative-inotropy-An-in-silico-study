/*
 *  integrate_rk1_halfdt.h
 *  
 *
 *  Created by Mao-Tsuen Jeng on 8/22/20.
 *
 *
 *  Use RK1 to compute ODE integration
 *  Integrate with order 1 accuracy
 *
 */

#ifndef integrate_rk1_H
#define integrate_rk1_H

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
// #include <omp.h>

void integrate_rk1_halfdt( void (*f)( double, double, double *, int,  double *,  double * , SimState1 * ), double * tspan, int sy0, double * y0, double dt, int p, double * mycurrents, SimState1 *S );

//void integrate_rk2( void (*f)( double , double , double *, double * , double * , pars_rec * ), double * tspan, int sy0, double * y0, double dt, double * pin , pars_rec * pars1 );

// f : function that evaluate ydots ( 1st derivatives ).
// tspan : time interval for integration
// dt : time step size
// y0 : initial value
// t1 : array of saved time.
// sy0 : size of y0

using namespace std;

void integrate_rk1_halfdt( void (*f)( double, double, double *, int,  double * ,  double *, SimState1 * ), double * tspan, int sy0, double * y0, double dt, int p, double * mycurrents , SimState1 *S ){

// void integrate_rk2( void (*f)( double , double , double *, double *, double * , pars_rec * ), double * tspan, int sy0, double * y0, double dt, double * pin , pars_rec * pars1 ){
	
	cout.precision(16);
	
	double tol = 0.25;
	int counter, counter_sy;
	double sy;
	double t, tn, err_y;
	double y1[sy0], y2[sy0], z[sy0], s, sm;
    double yh[sy0], dth;
	double dy1[sy0], dy2[sy0];
	double ydot[sy0];
	// char runType[] = "ydot";
	int i, j, k;
	//	int st = floor( 0.5 + ( tspan[1] - tspan[0] ) / ( dt * sy ) );
	// cout << st << "\t" << dt << "\t" << sy << endl;
	
	for( i = 0; i < sy0; i++ ){
		y1[i] = y0[i];
        ydot[i] = 0;
	}
	t = tspan[0];
    dth = 0.5 * dt;
	
	f( dth, t, y1, p, ydot, mycurrents , S );
	
	for( i = 0; i < sy0; i++ ) {
		dy1[i] = ydot[i] * dt; // k1
		y2[i] = y1[i] + dy1[i];
        yh[i] = y1[i] + 0.5 * dy1[i];
	}
    // cout << y0[54] << "\t" << y1[54] << "\t" << ydot[54] << "\t" << y2[54] << endl;
	
	
	f( dth, (t+dth), yh, p, ydot, mycurrents , S );
	
	for( i = 0; i < sy0; i++ ) {
		z[i] = yh[i] + dth * ydot[i];
	}
	
	
	t += dt;
	// Result using RK2
	sm = 1;
	
	for( i = 0; i < sy0; i++ ) {
		
		if ( z[i] != y2[i] && z[i] == z[i] ) {
			s = pow( tol * dt / ( 2 * fabs( z[i] - y2[i] ) ) , 0.5 );
			// cout << i << "\t" << z[i] << "\t" << y2[i] << "\t" << ydot[i] << endl;
			
			if( sm > s && s == s ) {
				
				sm = s;
				if( 0.1 > sm ) {
					sm = 0.1;
					// cout << z[i] << "\t" << y2[i] << endl;
				}
				// cout << s[i] << "\t";
				
			} else if ( s != s ) { // if s is NAN
				// cout << s << endl;
				sm = 0.1;
			}
			//cout << endl;
		} else if ( z[i] != z[i] ){
			cout << i << "\t" << z[i] << endl;  // z[i] is NAN
            cout << y1[i] << "\t" << y2[i] << endl;
            exit(0);
            sm = 0.1;
		}
	}
	if( sm >= 1 ){
		for( i = 0; i < sy0; i++ ) {
			// cout << "dt = " << dt << endl;
			y1[i] = z[i];
			// y1[i] = y2[i];		
		}
	} else {
		
		sy = 2 * ceil( 1./ sm );
		//	sy = 1;
		// cout << "Divided by " << sy << " parts.  sm = " << sm  << " at t = " << t << endl;
		
		dt = dt / sy;
		t = tspan[0];
		for( i = 0; i < sy0; i++ ){
			y1[i] = y0[i];
		}
		
		for( counter_sy = 0; counter_sy < sy; counter_sy ++ ) {
			
			// f( y1,  ydot , Istim , currents, celltype, gendertype, gterm , dt );
			f( dt, t, y1, p, ydot, mycurrents , S );
			
			for( i = 0; i < sy0; i++ ) {
				y1[i] = y1[i] + dt * ydot[i] ;
			}
			
			t += dt;
			
		}
		// cout << endl;
	}
	
	for( i = 0; i < sy0; i++ ) {
		y0[i] = y1[i];
		// y1[i] = y2[i];		
	}
}

#endif
