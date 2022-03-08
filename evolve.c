 #include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

FILE *input1;
FILE *input2;
FILE *output;

#define MAXPNT 3000
#define kND 3

/* define global arrays to store positions */
double positions[MAXPNT][kND],velocities[MAXPNT][kND];
double r1[kND],r2[kND],v1[kND],v2[kND];

/* define global arrays to store accel */
double accel[MAXPNT][kND];
double a1[kND],a2[kND];


double M,Rinit;
double G=1.0; 
double Rmin=25;   
double e=0.6; //eccentricity 
double a,b; //axeses of ellipse 

void readInit();
void readCore();
void coreaccel();
void diskparticleaccel();
void diskparticlesInteract();
void printstate();
// void printEnergy();


int main(argc, argv)
int argc;
char *argv[];
{
    M=pow(10,11); /* in solar mass */
    /* to keep G as 1, convert M into a new mass unit */
    M=M*4.30091*pow(10,-6);

    int n, mstep, nout, nstep;
    double eta, tmax, episqr;
    double dt;
    double tnow;
    int particleN;

    mstep = 5;            
    nout = 1;                 
    dt = pow(10,8);            
    tmax = dt*40;
    particleN = 297*2;
    tnow=-16.4;

    input1 = fopen("initdisk.dat", "r");
    input2 = fopen("initcore.dat", "r");
    output = fopen("evolution.dat", "w+");

    readInit(particleN);
    readCore();
    diskparticlesInteract();

    for (nstep = 0; nstep < mstep; nstep++) {   /* loop mstep times in all  */
        if (nstep % nout == 0)          /* if time to output state  */
        {
            printstate(particleN); /* then call output routine */
        }
        diskparticlesInteract(); /* take integration step    */
        tnow = tnow + 1;        /* and update value of time */
    }
    if (mstep % nout == 0) /* if last output wanted    */
    {
        printstate(particleN); /* output last step */
    }              

    fclose(input1);
    fclose(input2);
    fclose(output);

    return 0;
}

void readCore()
{
    int j;
/* the following reads in the core masses */
    for (j=0;j<kND;j++)
    {
        fscanf(input2,"%lf",&r1[kND]);
    }
    for (j=0;j<kND;j++)
    {
        fscanf(input2,"%lf",&v1[kND]);
    }
    for (j=0;j<kND;j++)
    {
        fscanf(input2,"%lf",&r2[kND]);
    }
    for (j=0;j<kND;j++)
    {
        fscanf(input2,"%lf",&v2[kND]);
    }

    printf("The core positions are: %-14.4E\n",r1[0]);

}


void readInit(n)
int n;
{
    int i,j;
    /* this reads the initial positions and velocities of the disk particles */
    for (i=0;i<n;i++)
    {
        for (j=0;j<kND;j++)
        {
            fscanf(input1,"%lf",&positions[i][j]);
        }
        for (j=0;j<kND;j++)
        {
            fscanf(input1,"%lf",&velocities[i][j]);
        }
    }
}


void getrelCofM(r1,r2,r3,n)
int n;
double r1[],r2[],r3[];
{

    int i;
    for (i = 0; i < n; i++)	{		
	r3[i] = r1[i]-r2[i];
    }

}

/* Compute the acceleration of the core masses */
void coreaccel()
{
    int j;
    double arelativ[kND];
    double dist;
    
    dist=pow(r1[0]-r2[0],2)+pow(r1[1]-r2[1],2)+pow(r1[2]-r2[2],2);
    dist=sqrt(dist);
    for (j=0;j<kND;j++)
    {
       arelativ[j]=-G*M*2*r1[j]/pow(dist,3);
       a1[j]=arelativ[j]/2;
       a2[j]=-a1[j];
    }
    
}


/* Compute the acceleration of the disk particles */
void diskparticleaccel(n)                
int n;
{
    int i,j;
    /* Initialize acceleration arrays */
    for (i=0;i<n;i++){ 
        for (j=0;j<kND;j++){
            accel[i][j]=0;
        } 
    }

    for (i=0;i<n;i++){ 
        double ri1=0;
        double ri2=0;

        for (j=0;j<kND;j++)
        {
            ri1=ri1+pow(r1[j]-positions[i][j],2);
        }
        ri1=sqrt(ri1);
        for (j=0;j<kND;j++)
        {
            ri2=ri2+pow(r2[j]-positions[i][j],2);
        }
        ri2=sqrt(ri2);

        for (j=0;j<kND;j++)
        {
           accel[i][j]=accel[i][j]-G*M*(positions[i][j]-r1[j])/pow(ri1,3)-G*M*(positions[i][j]-r2[j])/pow(ri2,3);
        }
    }
}


/* This is the leapstep function */
void diskparticlesInteract(n,dt)           
int n;                      
double dt;                  
{
    int i,j;
    coreaccel();
    diskparticleaccel(n);
    /* advance disk particles */
    for (i = 0; i < n; i++)        
    {   
        for (j=0;j<kND;j++)
        {
            velocities[i][j]=velocities[i][j]+0.5*dt*accel[i][j];
        }
    }
    for (i = 0; i < n; i++)        
    {
        for (j=0;j<kND;j++)
        {
            positions[i][j]=positions[i][j]+dt*velocities[i][j];
        }
    }

    /* advance core masses */
    for (j=0;j<kND;j++)
    {
        v1[j]=v1[j]+0.5*dt*a1[j];
        v2[j]=v2[j]+0.5*dt*a2[j];
    }
    for (j=0;j<kND;j++)
    {
        r1[j]=r1[j]+dt*v1[j];
        r2[j]=r2[j]+dt*v2[j];
    }

    coreaccel();
    diskparticleaccel(n);      
    
    /* finish the steps */
    for (i = 0; i < n; i++)        
    {   
        for (j=0;j<kND;j++)
        {
            velocities[i][j]=velocities[i][j]+0.5*dt*accel[i][j];
        }
    }
    for (j=0;j<kND;j++)
    {
        v1[j]=v1[j]+0.5*dt*a1[j];
        v2[j]=v2[j]+0.5*dt*a2[j];
    }
}

/* Generate a dat file for glnemo */
void printstate(n)
int n;
{
    int i;
    for (i=0;i<n;i++)
    {
        fprintf(output,"%-14.4E%-14.4E%-14.4E\n",positions[i][0],positions[i][1],positions[i][2]);
    }
    fprintf(output,"%-14.4E%-14.4E%-14.4E\n",r1[0],r1[1],r1[2]);
    fprintf(output,"%-14.4E%-14.4E%-14.4E\n",r2[0],r2[1],r2[2]);
}


//void printEnergy(E, fp, u, v, w, tnow, n)
//double E[];
//FILE *fp;
//double u[];
//double v[];
//double w[];
//double tnow;
//int n;
//{
//    double Energy=0;
//    double uMomentum=0;
//    double vMomentum=0;
//    double wMomentum=0;
//    double MomentumTotal=0;
//    for (int i=0;i<n;i++)
//    {
//        Energy=Energy+E[i]; /* (parsec/year)^2 */
//        uMomentum=uMomentum+m*u[i];
//        vMomentum=vMomentum+m*v[i];
//        wMomentum=wMomentum+m*w[i];
//    }
//    Energy=Energy*m1; /* in newmass*(parsec/yr)^2 */
//    MomentumTotal=sqrt(pow(uMomentum,2.0)+pow(wMomentum,2.0)+pow(vMomentum,2.0));
//    fprintf(fp, "%-14.4f%-17.7E%-17.7E\n",tnow,Energy,MomentumTotal);
//}
//
