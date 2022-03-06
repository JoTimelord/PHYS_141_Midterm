/*
 * LEAPINT.C: program to integrate hamiltonian system using leapfrog.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAXPNT 40000            /* maximum number of points */
#define PI 3.141592654

double M,G,R,m;

void leapstep();                /* routine to take one step */

void accel();                   /* accel. for harmonic osc. */

double rand_0_1();              /* generate random numbers */

void initial();                 /* generate initial conditions */

double g();                     /* generate q distribution */

void printstate();              /* print out system state   */

void printenergy();             /* verify that the energy of the system does not change */

void main(argc, argv)
int argc;
char *argv[];
{
    int n, mstep, nout, nstep;
    double eta, tmax, episqr;
    double x[MAXPNT], y[MAXPNT], z[MAXPNT], r[MAXPNT], V[MAXPNT], v[MAXPNT], w[MAXPNT], u[MAXPNT], E[MAXPNT], tnow, dt;
    double Eanalytic;

    /* define values of global variables */
    M=pow(10.0,11)*4.30091*pow(10,-3)*pow(1.02201*pow(10,-6),2); /* in new mass */
    G=1; /* in parsec^3/(newmass*yr^2) */
    R=1.5*1000; /* in parsec */
    m=M/20000.0; /* mass per star */

    /* set Aaserth Parameters */
    eta=0.02;
    n=20000; 
    episqr=0.25;
    Eanalytic=-3*PI/64*G*M*M/R; /* in newmass*(parsec/yr)^2 */

    /* next, set integration parameters */
    mstep = 5;                /* number of steps to take  */
    nout = 1;                   /* steps between outputs    */
    dt = 2000;            /* timestep for integration, in year */
    tmax=dt*100000;


    /* store init.dat */
    FILE *fp;
    fp=fopen("init.dat","w+");

    /* first, set up initial conditions */
    initial(r,x,y,z,V,v,w,u,E,n); /* set initial vel and posi of all points */
    fprintf(fp,"%-7s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s\n","index","r component","x component","y component","z component","speed","vx component","vy component","vz component","Energy");
    for (int i=0;i<n;i++) 
    {
        fprintf(fp,"%-7d%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f%-14.4E\n",i,r[i],x[i],y[i],z[i],V[i],u[i],v[i],w[i],E[i]); /* this prints out the index of the star, r, x, y, z, V, v, w, u, v component of the star */
    }
    fclose(fp); 


//    /* write inital condition files for aaserth process */
//    FILE *fp2;
//    fp2=fopen("init_aaserth.data","w+");
//    fprintf(fp2,"%-7d%-14.4f%-14.4f%-14.4f%-14.4f\n",n,eta,dt,tmax,episqr);
//    for (int i=0;i<n;i++)
//    {
//        fprintf(fp2,"%-14.4E%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f\n",m,x[i],y[i],z[i],u[i],v[i],w[i]);
//    } /* mass, x, y, z, vx, vy, vz */
//    fclose(fp2);
//
    tnow = 0.0;                 /* set initial time         */

//    FILE *fp3;
//    fp3=fopen("initplummer.data","w+");
//    printstate(x, y, z, u, v, w, n, tnow, fp3);
//    fclose(fp3);
//
    /* This is the file to check that the energy produced matches that of equation 6 (for problem a and b) */
    FILE *fp4;
    fp4=fopen("energy.dat","w+");
    fprintf(fp4,"This file contains the energy and totol momentum calculated from leapfrog code with time step as 5 years and total steps as 5.\n");
    fprintf(fp4,"This is the analytic energy of the plummer sphere over time: %14.4e\n",Eanalytic);
    fprintf(fp4,"The energy is in newmass*(parsec/yr)^2, and the momentum is in newmass*parsec/yr.\n");
    fprintf(fp4,"This is the check for momentum and energy conservation, as well as the fact that leapfrog produces energy as analytic form predicts.\n");
    fprintf(fp4,"%-14.4s%-14.4s%-14.4s\n","timenow","Energy","tot_moment");

    /* now, loop performing integration */
    FILE *fp5;
    fp5=fopen("evolution.dat","w+");
    fprintf(fp5,"%-7d%-14.4f\n",n,tnow);

    for (nstep = 0; nstep < mstep; nstep++) {   /* loop mstep times in all  */
        if (nstep % nout == 0)          /* if time to output state  */
        {
            printenergy(E, fp4, u, v, w, tnow, n);   
            printstate(x, y, z, u, v, w, n, fp5); /* then call output routine */
        }
        leapstep(r, x, y, z, V, w, v, u, E, n, dt); /* take integration step    */
        tnow = tnow + dt;           /* and update value of time */
    }
    if (mstep % nout == 0) /* if last output wanted    */
    {
        printenergy(E, fp4, u, v, w, tnow, n); 
        printstate(x, y, z, u, v, w, n, tnow, fp5); /* output last step */
    }              
    fclose(fp4);
    fclose(fp5);
}

/* set up initial conditions 
 * The positions are in kiloparsecs
 * The energies are in (km/s)^2
 * The velocities are in (km/s)
 * */
void initial(r,x,y,z,V,v,w,u,E,n)
double r[];
double x[];
double y[];
double z[];
double V[];
double v[];
double u[];
double w[];
double E[];
int n;
{
    for (int i=0;i<n;i++){
        double q;
        double X1, X2, X3, X4, X5, X6, X7;
        double Ve;
        X1=rand_0_1();
        r[i]=sqrt(R*R*pow(X1,2.0f/3.0f)/(1-pow(X1,2.0f/3.0f))); /* radius in parsec */
        X2=rand_0_1();
        z[i]=(1-2*X2)*r[i]; /*  z component */
        X3=rand_0_1();
        x[i]=pow(r[i]*r[i]-z[i]*z[i],0.5)*cos(2*PI*X3); /* x component */
        y[i]=pow(r[i]*r[i]-z[i]*z[i],0.5)*sin(2*PI*X3); /* y component */
        Ve=pow(2,0.5)*pow(1+r[i]*r[i]/(R*R),-0.25)*sqrt(G*M/R); /* escape velocity in parsec/yr */
        X4=rand_0_1();
        X5=rand_0_1();
        while (0.1*X5>=g(X4)) {
            X4=rand_0_1();
            X5=rand_0_1();
        }
        q=X4;
        X6=rand_0_1();
        X7=rand_0_1();
        V[i]=q*Ve; /* velocity component */
        w[i]=(1-2*X6)*V[i]; /* w component */
        u[i]=pow(V[i]*V[i]-w[i]*w[i],0.5)*cos(2*PI*X7); /* u component */
        v[i]=pow(V[i]*V[i]-w[i]*w[i],0.5)*sin(2*PI*X7); /* v component */
    }
    
    /* calculate energy */
    for (int j=0;j<n;j++) 
    {
        double U=0;
        for (int k=0;k<n;k++)
        {
            if (k>j){
                double r_jk=sqrt(pow(x[k]-x[j],2)+pow(y[k]-y[j],2)+pow(z[k]-z[j],2));
                U=U-G*m/r_jk;
            }
        }
        E[j]=U+V[j]*V[j]/2;
    }
}

/*
 * Probability distribution functions
 */
double g(q)
double q;
{
   return q*q*pow(1-q*q, 3.5); 
}

double rand_0_1(void)
{
    return (double) rand()/((double) RAND_MAX);
}

/*
 * LEAPSTEP: take one step using the leap-from integrator, formulated
 * as a mapping from t to t + dt.  WARNING: this integrator is not
 * accurate unless the timestep dt is fixed from one call to another.
 */

void leapstep(r, x, y, z, V, w, v, u, E, n, dt)
double r[];                 /* r-positions of all points  */
double x[];                 /* x-positions of all points  */
double y[];                 /* y-positions of all points  */
double z[];                 /* z-positions of all points  */
double V[];                 /* velocities of all points */
double w[];                 /* z velocities of all points */
double v[];                 /* y velocities of all points */
double u[];                 /* x velocities of all points */
double E[];                 /* energy of all points */
int n;                      /* number of points         */
double dt;                  /* timestep for integration */
{
    int i;
    double a[n],ax[n],ay[n],az[n];
    accel(a, ax, ay, az, r, x, y, z, n); /* call acceleration code   */
    for (i = 0; i < n; i++)         /* loop over all points...  */
    {
        V[i] = V[i] + 0.5 * dt * a[i];      /* advance vel by half-step */
        w[i] = w[i] + 0.5 * dt * az[i];      
        u[i] = u[i] + 0.5 * dt * ax[i];      
        v[i] = v[i] + 0.5 * dt * ay[i];      
    }
    for (i = 0; i < n; i++)         /* loop over points again...*/
    {
        r[i] = r[i] + dt * V[i];        /* advance pos by full-step */
        x[i] = x[i] + dt * u[i];        
        y[i] = y[i] + dt * v[i];        
        z[i] = z[i] + dt * w[i];        
    }
    accel(a, ax, ay, az, r, x, y, z, n);             /* call acceleration code   */
    for (i = 0; i < n; i++)             /* loop over all points...  */
    {
        V[i] = V[i] + 0.5 * dt * a[i];      /* advance vel by half-step */
        w[i] = w[i] + 0.5 * dt * az[i];      
        u[i] = u[i] + 0.5 * dt * ax[i];      
        v[i] = v[i] + 0.5 * dt * ay[i];      
    }
    /* calculate energy */
    for (int j=0;j<n;j++) 
    {
        double U=0;
        for (int k=0;k<n;k++)
        {
            if (k>j){
                double r_jk=sqrt(pow(x[k]-x[j],2)+pow(y[k]-y[j],2)+pow(z[k]-z[j],2));
                U=U-G*m/r_jk;
            }
        }
        E[j]=U+V[j]*V[j]/2;
    }
}

/*
 * ACCEL: compute accelerations for harmonic oscillator(s).
 */

void accel(a, ax, ay, az, r, x, y, z, n)
double a[];                 /* accelerations of points  */
double ax[];                 /* x accelerations of points  */
double ay[];                 /* y accelerations of points  */
double az[];                 /* z accelerations of points  */
double r[];                 /* r acceleration of points */
double x[];                 /* x acceleration of points */
double y[];                 /* y acceleration of points */
double z[];                 /* z acceleration of points */
int n;                      /* index of points         */
{
    for (int i=0;i<n;i++){ /* calculate acceleration for every point i */
        a[i]=0.0;
        ax[i]=0.0;
        ay[i]=0.0;
        az[i]=0.0;
        for (int j=0;j<n;j++){ /* calculate acceleration exerted on i by every other point */
            if (i!=j){
                double rij=sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2)+pow(z[i]-z[j],2)); /* in parsec */
                double aij=-G*m/(rij*rij);
                a[i]=a[i]+aij; /* in parsec per year */
                az[i]=az[i]+aij*(z[j]-z[i])/rij;
                ax[i]=ax[i]+aij*(x[j]-x[i])/rij;
                ay[i]=ay[i]+aij*(y[j]-y[i])/rij;
            }
        }
    }
}

/*
 * PRINTSTATE: output system state variables.
 */

void printstate(x, y, z, u, v, w, n, fp)
double x[];                 /* positions of all points  */
double y[];                 
double z[];                 
double u[];                 /* velocities of all points */
double v[];                 
double w[];                 
int n;
FILE *fp;                   /* the file name to store everything in */ 
{
    for (int i=0;i<n;i++){
        fprintf(fp,"%-14.4f%-14.4E%-14.4E%-14.4E%-14.4E%-14.4E%-14.4E\n",m,x[i],y[i],z[i],u[i],v[i],w[i]);
    }
}



/* 
 * printenergy: output system energy to make sure it stays the same.
 */

void printenergy(E, fp, u, v, w, tnow, n)
double E[];
FILE *fp;
double u[];
double v[];
double w[];
double tnow;
int n;
{
    double Energy=0;
    double uMomentum=0;
    double vMomentum=0;
    double wMomentum=0;
    double MomentumTotal=0;
    for (int i=0;i<n;i++)
    {
        Energy=Energy+E[i]; /* (parsec/year)^2 */
        uMomentum=uMomentum+m*u[i];
        vMomentum=vMomentum+m*v[i];
        wMomentum=wMomentum+m*w[i];
    }
    Energy=Energy*m; /* in newmass*(parsec/yr)^2 */
    MomentumTotal=sqrt(pow(uMomentum,2.0)+pow(wMomentum,2.0)+pow(vMomentum,2.0));
    fprintf(fp, "%-14.4f%-17.7E%-17.7E\n",tnow,Energy,MomentumTotal);
}

