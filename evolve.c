#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

FILE *input1;
FILE *input2;
FILE *outputFile;

#define PI 3.141592654
#define MAXPNT 3000
#define kND 3

/* define global arrays to store positions */
double positions[MAXPNT][kND],velocities[MAXPNT][kND];
double r1[kND],r2[kND],v1[kND],v2[kND];


double M,Rinit;
double G=1.0; 
double Rmin=25;   
double e=0.6; //eccentricity 
double a,b; //axeses of ellipse 

void readInit();
void getrelCofM();
void accel();
//void getInitialVelocities();
//void moveParticles();
void evolveGalaxies();
//void Interact();



int main(argc, argv)
int argc;
char *argv[];
{
    M=pow(10,11); /* in solar mass */
    /* to keep G as 1, convert M into a new mass unit */
    M=M*4.30091*pow(10,-6);

    int particleN;
    double r[MAXPNT],theta[MAXPNT]; 
    double x[MAXPNT],y[MAXPNT],z[MAXPNT]; 
    double Rinit,e,rp,w1,w2,i1,i2; 
    int n, mstep, nout, nstep;
    double eta, tmax, episqr;

    mstep = 5;            
    nout = 1;                 
    dt = 2000;            
    tmax = dt*100000;
    particleN = 297*2;

    input1 = fopen("initdisk.dat", "r");
    input2 = fopen("initcore.dat", "r");
    
    readInit(particleN);

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
    
    /* the following reads in the core masses */
    for (j=0;j<kND;j++)
    {
        fscanf(input2,"%lf",&r1[kND]);
        fscanf(input2,"%lf",&v1[kND]);
    }
    for (j=0;j<kND;j++)
    {
        fscanf(input2,"%lf",&r2[kND]);
        fscanf(input2,"%lf",&v2[kND]);
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

/* not sure if I should add "ar3" which is the acc of relative 
center of mass
*/
void accel(a, ar1, ar2, r1, , r2, n)
double a[];                
double ar1[];                
double ar2[];                         
double r1[];                 
double r2[];                
int n;                   
{
    for (int i=0;i<n;i++){ 
        a[i]=0.0;
        ar1[i]=0.0;
        ar2[i]=0.0;
        for (int j=0;j<n;j++){ 
            if (i!=j){
                double rij=sqrt(pow(r1[i]-r1[j],2)+pow(r2[i]-r2[j],2); 
                double aij=-G*m1*m2/(rij*rij);
                a[i]=a[i]+aij; 
                ar1[i]=ar1[i]-aij*(ar1[j]-ar1[i])/rij;
                ar2[i]=ar2[i]-aij*(ar2[j]-ar2[i])/rij;
            }
        }
    }
}




/*void moveParticles()
{
//not sure if this is needed
}*/



/* not sure if I should add "r3" which is the relative 
center of mass
*/

void evolveGalaxy(r1,r2,v1,v2,n,dt)
double r1[],r2[],v1[],v2[];
int n;
int dt;
{
int i;
    double a[n],ar1[n],ar2[n];
    accel(a, ar1, ar2, r1, r2, n); 
    for (i = 0; i < n; i++)         
    {
        v1[i] = v1[i] + 0.5 * dt * ar1[i];          
        v2[i] = v2[i] + 0.5 * dt * ar2[i];      
    }
    for (i = 0; i < n; i++)        
    {
        r1[i] = r1[i] + dt * v1[i];       
        r2[i] = r2[i] + dt * v2[i];             
    }
    accel(a, ar1,ar2,r1,r2,n);             
    for (i = 0; i < n; i++)             
    {
        v1[i] = v1[i] + 0.5 * dt * ar1[i];          
        v2[i] = v2[i] + 0.5 * dt * ar2[i];    
    }




//I will update the code so that it takes intoi acount 
//test particles instead of just r1 r2

/* void diskparticleaccel(a, ax, ay, az, r, x, y, z, n)
double a[];                \
double ax[];                 
double ay[];                 
double az[];             
double r[];                
double x[];                
double y[];                 
double z[];              
int n;                            
{
    for (int i=0;i<n;i++){ 
        a[i]=0.0;
        ax[i]=0.0;
        ay[i]=0.0;
        az[i]=0.0;
        for (int j=0;j<n;j++){ 
            if (i!=j){
                double rij=sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2)+pow(z[i]-z[j],2)); 
                double aij=-G*m/(rij*rij);
                a[i]=a[i]+aij; 
                az[i]=az[i]+aij*(z[j]-z[i])/rij;
                ax[i]=ax[i]+aij*(x[j]-x[i])/rij;
                ay[i]=ay[i]+aij*(y[j]-y[i])/rij;
            }
        }
    }
}


void diskparticlesInteract(r, x, y, z, V, w, v, u, E, n, dt)
double r[];                
double x[];                 
double y[];                
double z[];                
double V[];               
double w[];                
double v[];               
double u[];               
double E[];               
int n;                      
double dt;                  
{
    int i;
    double a[n],ax[n],ay[n],az[n];
    accel(a, ax, ay, az, r, x, y, z, n); 
    for (i = 0; i < n; i++)        
    {
        V[i] = V[i] + 0.5 * dt * a[i];     
        w[i] = w[i] + 0.5 * dt * az[i];      
        u[i] = u[i] + 0.5 * dt * ax[i];      
        v[i] = v[i] + 0.5 * dt * ay[i];      
    }
    for (i = 0; i < n; i++)        
    {
        r[i] = r[i] + dt * V[i];        
        x[i] = x[i] + dt * u[i];        
        y[i] = y[i] + dt * v[i];        
        z[i] = z[i] + dt * w[i];        
    }
    accel(a, ax, ay, az, r, x, y, z, n);             
    for (i = 0; i < n; i++)            
    {
        V[i] = V[i] + 0.5 * dt * a[i];     
        w[i] = w[i] + 0.5 * dt * az[i];      
        u[i] = u[i] + 0.5 * dt * ax[i];      
        v[i] = v[i] + 0.5 * dt * ay[i];      
    }
*/


}
