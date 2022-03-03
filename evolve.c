#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* Variables will go here
...
...
...
...
*/ 
double e = 0.6; //eccentricity 
double a, b ; //axeses of ellipse 
FILE * inputFile;



void readInit();
void getCofM();
void getrelCofM();
void getAcc();
void getInitialVelocities();
//void moveParticles();
void evolveGalaxies();
void Interact();



int main(argc, argv)
int argc;
char *argv[];
{
inputFile = fopen("initc.data", "r");


}


void readInit()
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

void getCofM(x,y,z, r1,r2)
double x[],y[],z[],r1[],r2[];
{
/* will get Center of Mass of each Galaxy*/


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

void getAcc()
{

}


void getInitialVelocities(x,y,z)

{
// v squared = r * a?
}


/*void moveParticles()
{
//not sure if this is needed
}*/

void evolveGalaxies(x,y,z,vx,vy,vz,t)
double x[],y[],z[],vx[],vy[],vz[],t[];
{

/* Im thinking of evolving them usinf the leapfrog code with
the center of mass, like advancing r3 as the positions with 
respect to R1 and R2. But IDK. 
Maybe Rugnna Kutta works better...*/ 

}