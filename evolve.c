 #include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>


double e = 0.6; //eccentricity 
double a, b ; //axeses of ellipse 
FILE * inputFile1;
FILE * inputFile2;
//these values are only as of right now 
m1  = 1;
m2 = 1;
m3 = 1;

#define MAXPNT 3000            
#define PI 3.141592654

double M,Rinit;
double G=1.0; 
double Rmin=25;   

void readInit();
void readInitdisk();
void readInitcore();
void getCofM();
void getrelCofM();
void diskparticleaccel();
void diskparticlesInteract();
void printEnergy();

int main(argc, argv)
int argc;
char *argv[];
{
 M=pow(10,11); 

    M=M*4.30091*pow(10,-6);
    

   /* int particleN;
    double r[MAXPNT],theta[MAXPNT]; 
    double x[MAXPNT],y[MAXPNT],z[MAXPNT]; 
    double Rinit,e,rp,w1,w2,i1,i2; 
     int n, mstep, nout, nstep;
    double eta, tmax, episqr;*/

    mstep = 5;            
    nout = 1;                 
    dt = 2000;            
    tmax=dt*100000;

inputFile1 = fopen("initdisk.data", "r");
inputFile2 = fopen("initcore.data", "r");


}

/* I assume this part is useless now??
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
*/

void readInitdisk(d1x,d1y,d1z,d1vx,d1vy,d1vz,d2x,d2y,d2z,d2vx,d2vy,d2vz)
double d1x[], d1y[], d1z[], d1vx[], d1vy[], d1vz[];
double d2x[], d2y[], d2z[], d2vx[], d1vy[], d2vz[];
{
    int i = 0;
    for (i=0; i >= 297; i++){
    i = i+1;
    fscanf(inputFile1, "%lf", &d1x);
    fscanf(inputFile1, "%lf", &d1y);
    fscanf(inputFile1, "%lf", &d1z);
    fscanf(inputFile1, "%lf", &d1vx);
    fscanf(inputFile1, "%lf", &d1vy);
    fscanf(inputFile1, "%lf", &d1vz);
    printf("\n");
}
    for (i=298; i >= 594; i++){
    i=i+1;
    fscanf(inputFile1, "%lf", &d2x);
    fscanf(inputFile1, "%lf", &d2y);
    fscanf(inputFile1, "%lf", &d2z);
    fscanf(inputFile1, "%lf", &d2vx);
    fscanf(inputFile1, "%lf", &d2vy);
    fscanf(inputFile1, "%lf", &d2vz);
    printf("\n");
}


}

void readInitcore(g1x,g1y,g1z,g1vx,g1vy,g1vz,g2x,g2y,g2z,g2vx,g2vy,g2vz)
double g1x[], g1y[], g1z[], g1vx[], g1vy[], g1vz[];
double g2x[], g2y[], g2z[], g2vx[], g1vy[], g2vz[];
{
    int i = 0;
    for (i=0; i >=1 ; i++){
    i = i+1;
    fscanf(inputFile1, "%lf", &g1x);
    fscanf(inputFile1, "%lf", &g1y);
    fscanf(inputFile1, "%lf", &g1z);
    fscanf(inputFile1, "%lf", &g1vx);
    fscanf(inputFile1, "%lf", &g1vy);
    fscanf(inputFile1, "%lf", &g1vz);
    printf("\n");
}
    for (i=1; i >=2; i++){
    i=i+1;
    fscanf(inputFile1, "%lf", &g2x);
    fscanf(inputFile1, "%lf", &g2y);
    fscanf(inputFile1, "%lf", &g2z);
    fscanf(inputFile1, "%lf", &g2vx);
    fscanf(inputFile1, "%lf", &g2vy);
    fscanf(inputFile1, "%lf", &g2vz);
    printf("\n");
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






void diskparticleaccel(ax, ay, az, g1x, g1y, g1z, g2x,g2y,g2z,x, y, z, n)                
double ax[];                 
double ay[];                 
double az[];             
double g1x[];
double g1y[];
double g1z[];
double g2x[];
double g2y[];
double g2z[];             
double x[];                
double y[];                 
double z[];              
int n;                            
{
    for (int i=0;i<n;i++){ 
        ax[i]=0.0;
        ay[i]=0.0;
        az[i]=0.0;
        for (int j=0;j<n;j++){ 
            if (i!=j){
                double rij1=sqrt(pow(x[i]-g1x[j],2)+pow(y[i]-g1y[j],2)+pow(z[i]-g1z[j],2)); 
                double rij2=sqrt(pow(x[i]-g2x[j],2)+pow(y[i]-g2y[j],2)+pow(z[i]-g2z[j],2)); 

                ax[i]=ax[i]- G*((-m1)*(x[i]-g1x[j])/rij1*rij1*rij1)*((-m2)*(x[i]-g2x[j])/rij2*rij2*rij2);
                ay[i]=ay[i]- G*((-m1)*(y[i]-g1y[j])/rij1*rij1*rij1)*((-m2)*(y[i]-g2y[j])/rij2*rij2*rij2);
                az[i]=az[i]- G*((-m1)*(z[i]-g1z[j])/rij1*rij1*rij1)*((-m2)*(z[i]-g2z[j])/rij2*rij2*rij2);
            }
        }
    }
}


void diskparticlesInteract(x, y, z, w, v, u, E, n, dt, g1x,g1y,g1z,g2x,g2y,g2z)           
double x[];                 
double y[];                
double z[];                            
double w[];                
double v[];               
double u[];               
double E[];               
double g1x[];
double g1y[];
double g1z[];
double g2x[];
double g2y[];
double g2z[];
int n;                      
double dt;                  
{
    int i;
    double ax[n],ay[n],az[n];
    diskparticleaccel(ax, ay, az, g1x, g1y, g1z, g2x,g2y,g2z,x, y, z, n);
    for (i = 0; i < n; i++)        
    {   
        w[i] = w[i] + 0.5 * dt * az[i];      
        u[i] = u[i] + 0.5 * dt * ax[i];      
        v[i] = v[i] + 0.5 * dt * ay[i];      
    }
    for (i = 0; i < n; i++)        
    {

        x[i] = x[i] + dt * u[i];        
        y[i] = y[i] + dt * v[i];        
        z[i] = z[i] + dt * w[i];        
    }
    diskparticleaccel(ax, ay, az, g1x, g1y, g1z, g2x,g2y,g2z,x, y, z, n);      
    for (i = 0; i < n; i++)            
    {  
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

void printEnergy(E, fp, u, v, w, tnow, n)
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
    Energy=Energy*m1; /* in newmass*(parsec/yr)^2 */
    MomentumTotal=sqrt(pow(uMomentum,2.0)+pow(wMomentum,2.0)+pow(vMomentum,2.0));
    fprintf(fp, "%-14.4f%-17.7E%-17.7E\n",tnow,Energy,MomentumTotal);
}
