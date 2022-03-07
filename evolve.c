#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>


double e = 0.6; //eccentricity 
double a, b ; //axeses of ellipse 
FILE * inputFile;
G = 1; 
m1  = 1;
m2 = 1;

#define MAXPNT 3000            /* maximum number of points */
#define PI 3.141592654

double M,Rinit;
double G=1.0; 
double Rmin=25;   

void readInit();
void getCofM();
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
 M=pow(10,11); 

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
    tmax=dt*100000;

inputFile = fopen("initc.data", "r");


}


void readInit(n,Rmin,r1,r2,x,y,z)
int n; 
double Rmin;
double r1[], r2[], x[], y[], z[];
{
    fscanf(inputFile, "%i", &n);
    fscanf(inputFile, "%lf", &Rmin);

    //will probably have to do a loop like in Aeresth 
    fscanf(inputFile, "%lf", &r1);
    fscanf(inputFile, "%lf", &r2);
    fscanf(inputFile, "%lf", &x);
    fscanf(inputFile, "%lf", &y);
    fscanf(inputFile, "%lf", &z);
    printf("\n");
}

//void getCofM(x,y,z, r1,r2)
//double x[],y[],z[],r1[],r2[];
//{
/* will get Center of Mass of each Galaxy*/


//}




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
                ar1[i]=ar1[i]+aij*(ar1[j]-ar1[i])/rij;
                ar2[i]=ar2[i]+aij*(ar2[j]-ar2[i])/rij;
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
