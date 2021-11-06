#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>
#include <sys/time.h>
// #include <mpi.h>
#define const_fields 1
#define Np           2
#define dt           1e-3
#define pusher_name  3     // ['Boris relativistic': 0, 'Vay relativ':1, 'HC relativ':2, 'Lapenta_Markidis':3]

#define nsteps 1e5
#define PI 4*atan(1)


// Basic definitions  // 

// Constants (cgs units)
#define  mp    1.673e-24         // ion (proton) mass
#define  mp2e  1836.15267389     // proton to electron mass ratio
#define  qe    4.803e-10         // elementary charge
#define  c     2.998e10          // speed of light
#define  kB    1.3807e-16        // Boltzmann constant 

// Plasma properties //
#define  n0    1e10              // plasma number density
#define  T     1e6              // ion temperature

// Box
#define  lz  1e3           // z dimension
#define  lx  1e1           // x dimension
#define  ly  1e1           //  y dimension

// Magnetic field, alfven speed
#define B0  1e4                        // ambient magnetic field, in +z direction

// frequencies
#define omegac  qe*B0/(mp*c)          // ion gyrofrequency

//-----------------------------------------------------------//

// Adopted code units //
// r0 (length): skin depth, vA/omegac = c/omegap
//#define r0 = vA/omegac;
// m0 (mass): ion mass, mi
//#define m0 = mp;
// t0 (time): inverse gyrofrequency, 1/Omegac
//#define t0 = 1./omegac;
// v0 (velocity): alfven speed, vA
//#define v0 = vA;
// B0 (magnetic field): B0
// E0 (electric field): vA*B0/c
//#define E0 = vA*B0/c;


static inline long double WTime(void){
    //timing the simulation for efficiency
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
}

long double dotP(long double vect_A[3], long double vect_B[3]){
//Estimate dot product of two vector arrays.
    int product = 0;
    // Loop for calculate cot product
    for (int i = 0; i < 3; i++)
        product = product + vect_A[i] * vect_B[i];
    return product;
}


void crossP(long double vect_A[3], long double vect_B[3], long double cross_P[3]){
 // Estimate cross product of two vector arrays.
 
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
}




/*  Relativistic requations */
long double gamma_u( long double p[3], long double c2){
    //    Lorentz factor
    // p: dimensionless momentum, i.e. gamma*v
    //c2: (vA/c)^2, vA the Alfven speed
    long double p2, gL;

     p2  = dotP(p,p);
     gL  = sqrtl(1.0L + p2/c2);
    return gL;
}


 long double gamma_v( long double v[3], long double c2){
    //    Lorentz factor
    // v: dimensionless velocity, i.e. v/vA
    //c2: (vA/c)^2, vA the Alfven speed
    long double v2, gL;

    v2  = dotP(v,v);
    gL  = 1.0L /sqrtl(1.0L - v2/c2);
    return gL;
}

 long double Ekin(long double p[3], long double c2){
    // Kinetic energy of particles 
    // (in units of mp*vA^2)
    long double gL;
    long double Ekinetic;

     gL       = gamma_u(p, c2);
     Ekinetic = (gL-1.0L) * c2;

    return Ekinetic;
}


void pm(long double v[3], long double gL, long double p[3]){
// Momentum 
// (normalized to m*vA) from velocity

     for(int i=0;i<3;i++){
        p[i] = v[i]*gL;
     }
}



/*
void interpolate_fields(){

    // Get E, B at new position
      call get_vec(bp, particle(ipart)%igrid, 
           particle(ipart)%self%x,particle(ipart)%self%time,b)
      call get_vec(vp, particle(ipart)%igrid, 
           particle(ipart)%self%x,particle(ipart)%self%time,vfluid)
      call get_vec(jp, particle(ipart)%igrid, 
           particle(ipart)%self%x,particle(ipart)%self%time,current)
      e(1) = -vfluid(2)*b(3)+vfluid(3)*b(2) + particles_eta*current(1)
      e(2) = vfluid(1)*b(3)-vfluid(3)*b(1) + particles_eta*current(2)
      e(3) = -vfluid(1)*b(2)+vfluid(2)*b(1) + particles_eta*current(3)

}
*/



void get_fields(long double E[3], long double B[3], long double x[3]){
    E[0] = 0.;
    E[1] = 0.;
    E[2] = 0.;
    B[0] = 0.;
    B[1] = 0.;
    B[2] = 1.0L;
  
    
}

void Vay(long double v[3], long double E[3], long double B[3], long double c2,  long double upart_new[3] ){

   long double  tau[3], emom[3], tmp[3], tmp_new[3] ;
   long double uprime[3], vprime[3], upart[3];
   long double sigma, gL;

   gL = gamma_v(v, c2);  //lorentz at start of time step
   pm(v, gL, upart);     // momentum at start of time step


   for(int i=0;i<3;i++){ 
       emom[i] = E[i] * dt /(2.0L );
    }

    for(int i=0;i<3;i++){
         tau[i]  =  B[i] * dt / (2.0L ); 
    }
 
    // Estimate crossP product of u, tau. 
    // Return tmp vector.
    crossP(upart, tau, tmp);
    for(int i=0;i<3;i++){
        uprime[i] = upart[i] + 2.0*emom[i] + tmp[i]/gL;
    }

    // New lorentz factor
    long double gL_n = gamma_u(uprime, c2);

    sigma = powl(gL_n,2) - dotP(tau,tau);   
    long double gL_n2 = sqrtl((sigma + sqrtl(powl(sigma,2) 
            + 4.0L * (dotP(tau,tau) + powl((dotP(uprime,tau)),2)))) / 2.0L);
    
    crossP(uprime,tau,tmp_new);

    // New Momentum
    for(int i=0;i<3;i++){
        upart_new[i] = (uprime[i] 
                   + dotP(uprime,tau)*tau[i]/powl(gL_n2,2) 
                   + tmp_new[i]/gL_n2) / (1.0L + dotP(tau,tau)/powl(gL_n2,2) );
    }
}


void Higuera_Cary(int part, int ntime, long double v[3], long double E[3], long double B[3], long double c2, long double upart_new[3] ){

   long double  tau[3], emom[3], tmp[3], tmp_new[3] ;
   long double uprime[3], vprime[3], upart[3];
   long double gL, gL_n, gL_n2, sigma;

   
   gL = gamma_v(v, c2);  // Lorentz at start of time step
   
   pm(v, gL, upart);     // momentum at start of time step

   if (part==0 && ntime<20){
       printf("p =%.100Lg\n",upart[0]);
   }

   for(int i=0;i<3;i++){ 
       emom[i] = E[i] * dt /(2.0L);
    }

    for(int i=0;i<3;i++){
         tau[i]  = B[i] * dt / (2.0L ); 
    }
 
    // Estimate crossP product of u, tau. 
    // Return tmp vector.
    //crossP(upart, tau, tmp);
    for(int i=0;i<3;i++){
        uprime[i] = upart[i] + emom[i];
    }

// New lorentz factor
    gL_n = gamma_u(uprime, c2);
    

    sigma = powl(gL_n,2) - dotP(tau,tau);   
    gL_n2 = sqrtl((sigma + sqrtl(powl(sigma,2) 
            + 4.0L * (dotP(tau,tau) + powl((dotP(uprime,tau)),2)))) / 2.0L);
    
    crossP(uprime,tau,tmp_new);

    // New Momentum
    for(int i=0;i<3;i++){
        upart_new[i] = (uprime[i] 
                   + dotP(uprime,tau)*tau[i]/powl(gL_n2,2) 
                   + tmp_new[i]/gL_n2) / (1.0L + dotP(tau,tau)/powl(gL_n2,2)) + emom[i] + tmp_new[i]/gL_n2;
    }
}


void Boris_relativ( long double v[3], long double E[3], long double B[3], long double c2, long double upart_new[3] ){

   long double  tau[3], emom[3], tmp[3], tmp_new[3] ;
   long double uprime[3], vprime[3], upart[3];
   long double sigma[3], uplus_1[3], uplus[3];

   long double gL = gamma_v(v, c2);
   

    for(int i=0;i<3;i++){       
        upart[i] = v[i]*gL;       // momentum at start of time step
    }

   for(int i=0;i<3;i++){ 
       emom[i] =  E[i] * dt /(2.0L );
    }

  
 
    // Estimate crossP product of u, tau. 
    // Return tmp vector.
    //crossP(upart, tau, tmp);
    for(int i=0;i<3;i++){
        uprime[i] = upart[i] + emom[i];
    }

// New lorentz factor
    long double gL_n = gamma_u(uprime, c2);

    for(int i=0;i<3;i++){
        tau[i]  =  B[i] * dt / (2.0L * gL_n ); 
    }
    for(int i=0;i<3;i++){
        sigma[i] = 2*tau[i]/(1.0L + dotP(tau,tau)); 
    } 

    crossP(uprime,tau,tmp);
    for(int i=0;i<3;i++){
        uplus_1[i] = uprime[i] + tmp[i]; 
    } 
    crossP(uplus_1,sigma,tmp_new);

    for(int i=0;i<3;i++){
        uplus[i] = uprime[i] + tmp_new[i]; 
    } 

    

    // New Momentum
    for(int i=0;i<3;i++){
        upart_new[i] = uplus[i] + emom[i];
    }
}


void Lapenta_Markidis( long double v[3], long double E[3], long double B[3], long double c2, long double upart[3] ){
   
   long double   upartk[3], vbar[3];
   long double  tmp[3], Fk[3], C1[3], C2[3] ;
   long double  dupartk[3];

   long double gL = gamma_v(v, c2);
   

    for(int i=0;i<3;i++){       
        upart[i] = v[i]*gL;       // momentum at start of time step
        upartk[i] = upart[i];
    }


    /* Start of the nonlinear cycle */

    long double  abserr = 1.0;
    long double  tol    = 1e-14;
    int     nkmax  = 10;
    int     nk     = 0; //


    do{
        long double  J11, J12, J13,J21, J22,  J23, J31, J32, J33, Det;
        long double  iJ11, iJ12, iJ13,iJ21, iJ22,  iJ23, iJ31, iJ32, iJ33;
        long double  gL_new;

        nk     =  nk+1;
        gL_new = gamma_u(upartk, c2);

        for(int i=0;i<3;i++){  
            vbar[i] = (upart[i] + upartk[i])/(gL_new + gL);
        }
        
        crossP(vbar,B,tmp);

        // Compute residual vector
        for(int i=0;i<3;i++){ 
            Fk[i] = upartk[i] - upart[i] - dt * (E[i] + tmp[i]);
        }

         // Compute auxiliary coefficients
        for(int i=0;i<3;i++){
            C1[i] = (gL_new + gL - (upartk[i]*(upartk[i] + upart[i])) / (gL_new*c2) )/ powl((gL + gL_new),2) ;
            C2[i] = -( upartk[i] / (gL_new*c2) )/ ((gL_new + gL),2) ;
        }

    
        // Compute Jacobian
          J11 = 1. - dt * (C2[0] * (upartk[1] + upart[1]) * B[0] - C2[1] * (upartk[2] + upart[2]) * B[1]) ;
          J12 = - dt*(C1[1] * B[2] - C2[1] * (upartk[2] + upart[2]) * B[1]) ;
          J13 = - dt * (C2[2] * (upartk[1] + upart[1]) * B[2] - C1[2] * B[1]) ;
          J21 = - dt * (- C1[0] * B[2] + C2[0] * (upartk[2] + upart[2]) * B[0]) ;
          J22 = 1. - dt * (- C2[1] * (upartk[0] + upart[0]) * B[2] + C2[1] * (upartk[2] + upart[2]) * B[0]) ;
          J23 = - dt * (- C2[2] * (upartk[0] + upart[0]) * B[2] + C1[2] * B[0]) ;
          J31 = - dt * (C1[0] * B[1] - C2[0] * (upartk[1] + upart[1]) * B[0]) ;
          J32 = - dt * (C2[1] * (upartk[0] + upart[0]) * B[1] - C1[1] * B[0]) ;
          J33 = 1. - dt * (C2[2] * (upartk[0] + upart[0]) * B[1] - C2[2] * (upartk[1] + upart[1]) * B[0]) ;

          // Compute inverse Jacobian
           Det = J11*J22*J33 + J21*J32*J13 + J31*J12*J23 - J11*J32*J23 - J31*J22*J13 - J21*J12*J33;
           iJ11 = (J22*J33 - J23*J32) / Det ;
           iJ12 = (J13*J32 - J12*J33) / Det ;
           iJ13 = (J12*J23 - J13*J22) / Det ;
           iJ21 = (J23*J31 - J21*J33) / Det ;
           iJ22 = (J11*J33 - J13*J31) / Det ;
           iJ23 = (J13*J21 - J11*J23) / Det ;
           iJ31 = (J21*J32 - J22*J31) / Det ;
           iJ32 = (J12*J31 - J11*J32) / Det ;
           iJ33 = (J11*J22 - J12*J21) / Det ;

        // Compute new upartk = upartk - J^(-1) * F(upartk)
          dupartk[0] = - (iJ11 * Fk[0] + iJ12 * Fk[1] + iJ13 * Fk[2]);
          dupartk[1] = - (iJ21 * Fk[0] + iJ22 * Fk[1] + iJ23 * Fk[2]);
          dupartk[2] = - (iJ31 * Fk[0] + iJ32 * Fk[1] + iJ33 * Fk[2]);

        // Check convergence
        for(int i=0;i<3;i++){
          upartk[i] +=  dupartk[i] ;
        }
          abserr = sqrtl(dotP(dupartk, dupartk));
    } while(abserr > tol && nk < nkmax); // End of while -> end of cycle

    // Update velocity 
    for(int i=0;i<3;i++){
        upart[i] = upartk[i];
    }
   
    
}





void pusher(int part, int ntime,long double x[3], long double v[3], long double E[3], long double B[3], long double Vth, long double c2, long double x_new[3], long double v_new[3]){
    /*
    1) First, the position of a particle is advanced by half time step (1st-order Euler scheme)
    2) Then, the velociy is updated by full time step using the positions at half time step
    3) Finally, the new position at full step size is computed using the new velocity
    */
   long double x_mid[3];
   long double upart_new[3];

    // half time-step for Euler step
     long double dth = 0.5*dt;
    
    // position x_{i+1/2} at time t_{i+1/2} + dt/2 (1st-order Euler scheme)
    for(int i=0;i<3;i++){
        x_mid[i] = x[i] + v[i]*dth;
    }
 
   // TO DO: Get E, B at new position (x^{n+1/2})
   // interpolate_fields();
    get_fields(E, B, x);

   // 'Kick' particle (update velocity) based on the chosen integrator
    switch(pusher_name) {
    case 1 :
        Vay(v, E, B, c2, upart_new);
        break;
        
    case 0  :
        Boris_relativ(v, E, B, c2, upart_new);
        break; 
        
    case 2  :
        Higuera_Cary(part, ntime,v, E, B, c2, upart_new);
        break;        

    case 3 :
        Lapenta_Markidis(v, E, B, c2,  upart_new);
        break;          
    }

    long double gL_new = gamma_u(upart_new, c2);
    // full position update at t_{i+1} = t_{i+1/2} + dt/2 = t_i + dt
   for(int i=0;i<3;i++){
        v_new[i]  = upart_new[i]/gL_new; 
        x_new[i]  = x_mid[i] + v_new[i]*dth;
        
   }

}
/* End of Relativistic requations*/





void fill_vector(long double X[3][Np], long double x[3], int particle_number){
    //used for shared-memory OpenMP caculations
    x[0] = X[0][particle_number];
    x[1] = X[1][particle_number];
    x[2] = X[2][particle_number];
}

void update_matrix(long double X[3][Np], long double x[3], int particle_number){
    //used for shared memory OpenMP calculations
	X[0][particle_number] = x[0];
    X[1][particle_number] = x[1];
    X[2][particle_number] = x[2];
}



 long double gaussian_number(long double Vth, long double vA){
    //obtain a gaussian number from two random numbers on a range of [0,1] inclusive.
    //this procedure CAN produce two gaussian numbers, but we only need one...
    //also, as we multiply by a standard deviation to get Maxwellian
    long double number1 = ( long double)rand() / ( long double)RAND_MAX ;
    long double number2 = ( long double)rand() / ( long double)RAND_MAX ;
    
    return Vth/(vA*sqrtl(3.)) * sqrtl(-2. * log(number1)) * sin(2. * PI * number2) ;
}

void fill_matrices(long double X[3][Np], long double V[3][Np], long double vA, long double Vth){
    //initialise the matrices for position and velocity 

    for (int j = 0; j< 3; j++){ 
        for (int i = 0; i < Np; i++) {
            V[j][i] = gaussian_number(Vth, vA);
            X[j][i] = lz* ((long double)rand() / (long double)RAND_MAX) - (lz / 2.0L);
           
        }
    }
}

int main(){
    //initialise position and velocity matrices
     long double X[3][Np];
     long double V[3][Np];


    long double omegap   = sqrtl((4*PI*n0*qe*qe)/mp);     // ion plasma frequency
    long double vA       = B0/sqrtl(4*PI*n0*mp);          // alfven speed
    long double Vth      = sqrtl(kB*T/mp);                // ion thermal speed
    long double c2       = powl((c/vA),2);
    
    //printf("%.100Lg",Vth);


    //FILL MATRICES:
    fill_matrices(X,V, vA, Vth);
      
    //Openfile
	FILE *fout1 = NULL;
	fout1 = fopen("particle_tracer.csv", "w");

	//TIMEING
    long double time1 = WTime();

    //TIME LOOP:
    for (int nt = 0; nt < nsteps; nt++){
    
#pragma omp parallel for

        //PARTICLE LOOP:
        for (int particle = 0; particle < Np; particle++){
            long double x[3]; 
            long double v[3];
            long double x_new[3];
            long double v_new[3];
            long double p_new[3];
            long double E[3];
            long double B[3];

            fill_vector(X,x,particle);
            fill_vector(V,v,particle);
            if (nt==0){
                fprintf(fout1,"%d, %g, %.100Lg, %.100Lg, %.100Lg, %.100Lg, %.100Lg, %.100Lg\n",particle,nt*dt,x[0],x[1],x[2],v[0],v[1],v[2]); 
                }
            
            get_fields(E, B, x);
            pusher( particle, nt,x, v, E, B, Vth, c2,  x_new,v_new);
            update_matrix(X,x_new,particle); 
            update_matrix(V,v_new,particle); 
            fprintf(fout1,"%d, %g, %.100Lg, %.100Lg, %.100Lg, %.100Lg, %.100Lg, %.100Lg\n",particle,(nt+1)*dt,x_new[0],x_new[1],x_new[2],v_new[0],v_new[1],v_new[2]); 
        }
    }

    fclose(fout1); 

    long double time2 = WTime();
    printf("TIME = %.100Lg\n",time2-time1);

    return 0;
} 