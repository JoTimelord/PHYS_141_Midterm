 /* This the initial setup for the two galaxies in mice project */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAXPNT 3000            /* maximum number of points */
#define PI 3.141592654

double M,Rinit;
double G=1.0; 
double Rmin=25; /* in kiloparsec */
double tapo;

void angular();
void position();
void corepossep();
void corepos(); 
void printdat();
void checkinit();

int main(argc, argv)
int argc;
char *argv[];
{
    M=pow(10,11); /* in solar mass */
    G=4.301*pow(10,-6)*pow(1.022*pow(10,-9),2); /* now G is in solar mass,kp,calendar year */

    
    /* set up initial positions */
    int particleN;
    double r[MAXPNT],theta[MAXPNT];
    double r1[MAXPNT],r2[MAXPNT]; /* find out the initial positions of the core masses */
    double v1[MAXPNT],v2[MAXPNT]; /* find out the initial velocity of the core masses */
    double x[MAXPNT],y[MAXPNT],z[MAXPNT]; /* the position vector w.r.t. the center of mass of the disk particles */
    double vx[MAXPNT],vy[MAXPNT],vz[MAXPNT];
    double Rinit,e,rp,w1,w2,i1,i2; /* the Kepler orbits parameters */

    e=0.6;
    rp=Rmin;
    w1=-90/PI;
    w2=-90/PI;
    i1=15/PI;
    i2=60/PI;
    particleN=297*2;  


    angular(r,theta);
    corepossep(r1,r2,rp,e,v1,v2);
    position(x,y,z,vx,vy,vz,r,theta,particleN,r1,r2,i1,i2,w1,w2,v1,v2);
    
    FILE *fp1;
    fp1=fopen("initdisk.dat","w+");
    FILE *fp2;
    fp2=fopen("initcore.dat","w+");
    printdat(r1,r2,v1,v2,x,y,z,vx,vy,vz,fp1,fp2,particleN);
    fclose(fp1);
    fclose(fp2);
    /* the following is to check the initial conditions */
    FILE *fp3;
    fp3=fopen("initplummer.dat","w+");
    checkinit(r1,r2,x,y,z,fp3,particleN);
    fclose(fp3);
    


    return 0;
}

void angular(r,theta)
double r[];
double theta[];
{
    int ringn=11;
    int innern=12;
    int outern=42;
    for (int i=0;i<ringn;i++)
    {
        double radius=(0.2+0.05*i)*Rmin;
        double the=2*PI/(innern+3*i);
        r[i]=radius;
        theta[i]=the;
    }
}


/* when applying the rotation matrix, we assume that the pericenter axis is the x-axis, and that the x-y plane is the orbit plane. */

void position(x,y,z,vx,vy,vz,r,theta,n,r1,r2,i1,i2,w1,w2,v1,v2)
double x[];
double y[];
double z[];
double vx[];
double vy[];
double vz[];
double r[];
double theta[];
int n;
double r1[];
double r2[];
double i1;
double i2;
double w1;
double w2;
double v1[];
double v2[];
{
    int ringn=11;
    int innern=12;
    int outern=42;
    double esoft=0.2;
    
    int count=0;
   
    /* first galaxy */
    for (int i=0;i<ringn;i++)
    {
        for (int j=0;j<12+3*i;j++)
        {
            double rb[3]; /* before rotation, the disk particles are on x-y plane */
            double vb[3];
            double vel;

            /* the speed of the disk particle */
            vel=sqrt(G*M*r[i]/(pow(r[i],2)+pow(esoft,2)));

            /* set up the disk particle originally in the x-y plane */
            rb[0]=r[i]*cos(j*theta[i]);
            rb[1]=r[i]*sin(j*theta[i]);
            rb[2]=0;
            vb[0]=-vel*sin(j*theta[i]);
            vb[1]=vel*cos(j*theta[i]);
            vb[2]=0;

            /* suppose the axis of intersection is along x-axis */
            /* first rotate around the x-axis by -i1 */

            double rc[3]; /* after the first rotation */
            double vc[3];

            rc[0]=rb[0];
            rc[1]=rb[1]*cos(-i1)-rb[2]*sin(-i1);
            rc[2]=rb[1]*sin(-i1)+rb[2]*cos(-i1);
            vc[0]=vb[0];
            vc[1]=vb[1]*cos(-i1)-vb[2]*sin(-i1);
            vc[2]=vb[1]*sin(-i1)+vb[2]*cos(-i1);

            /* now rotate the disk particles around z-axis by -w1 so that the angle between the pericenter axis and the intersection axis is w1 */
            double rd[3]; /* after the second rotation */
            double vd[3]; 
            
            rd[0]=cos(-w1)*rc[0]-sin(-w1)*rc[1];
            rd[1]=sin(-w1)*rc[0]+cos(-w1)*rc[1];
            rd[2]=rc[2];
            vd[0]=cos(-w1)*vc[0]-sin(-w1)*vc[1];
            vd[1]=sin(-w1)*vc[0]+cos(-w1)*vc[1];
            vd[2]=vc[2];
             
            /* Now so far it is still with respect to the core mass. Now make connections */
            x[count]=rd[0]+r1[0];
            y[count]=rd[1]+r1[1];
            z[count]=rd[2]+r1[2];

            vx[count]=vd[0]+v1[0];
            vy[count]=vd[1]+v1[1];
            vz[count]=vd[2]+v1[2];

            count++;
        }
    }

    /* second galaxy */
    for (int i=0;i<ringn;i++)
    {
        for (int j=0;j<12+3*i;j++)
        {
            /* do similar stuff as before */
            double rb[3]; /* before rotation, the disk particles are on x-y plane */
            double vb[3];
            double vel;

            vel=-sqrt(G*M*r[i]/(pow(r[i],2)+pow(esoft,2)));

            rb[0]=r[i]*cos(j*theta[i]);
            rb[1]=r[i]*sin(j*theta[i]);
            rb[2]=0;
            vb[0]=vel*sin(j*theta[i]);
            vb[1]=-vel*cos(j*theta[i]);
            vb[2]=0;


            /* suppose the axis of intersection is along x-axis */
            /* first rotate around the x-axis by -i2 */

            double rc[3]; /* after the first rotation */
            double vc[3];

            rc[0]=rb[0];
            rc[1]=rb[1]*cos(-i2)-rb[2]*sin(-i2);
            rc[2]=rb[1]*sin(-i2)+rb[2]*cos(-i2);
            vc[0]=vb[0];
            vc[1]=vb[1]*cos(-i2)-vb[2]*sin(-i2);
            vc[2]=vb[1]*sin(-i2)+vb[2]*cos(-i2);

            /* now rotate the disk particles around z-axis by -w1 so that the angle between the pericenter axis and the intersection axis is w1 */
            double rd[3]; /* after the second rotation */
            double vd[3];

            rd[0]=cos(-w2)*rc[0]-sin(-w2)*rc[1];
            rd[1]=sin(-w2)*rc[0]+cos(-w2)*rc[1];
            rd[2]=rc[2];
            vd[0]=cos(-w2)*vc[0]-sin(-w2)*vc[1];
            vd[1]=sin(-w2)*vc[0]+cos(-w2)*vc[1];
            vd[2]=vc[2];
            
            /* Now so far it is still with respect to the core mass. Now make connections */
            x[count]=rd[0]+r2[0];
            y[count]=rd[1]+r2[1];
            z[count]=rd[2]+r2[2];
            vx[count]=vd[0]+v2[0];
            vy[count]=vd[1]+v2[1];
            vz[count]=vd[2]+v2[2];

            count++;
        }
    }

    
    /* check if the number of particles in each disk is equal to 297 */
    if (count!=n)
    {
        printf("There is something wrong with how I compute (r,theta).\n");
    }
    
}

/* give the initial position (x,y,z) of the two core masses at t=-16.4 directly computed from the two ellipitcal orbits */
void corepossep(r1,r2,rp,e,v1,v2)
double r1[];
double r2[];
double rp;
double e;
double v1[];
double v2[];
{
    double a;
    double T; /* period of the kepler orbit */
    double p,b;
    double ra;
    double arealvel;
    double r1i[3];
    double r2i[3];

    p=rp*(1+e);
    b=p/sqrt(1-e*e);
    a=p/(1-e*e);
    ra=a*(1+e);

    /* Period */
    T=sqrt(4*PI*PI*pow(a,3)/(G*2*M)); /* this is in calendar year */

    /* orbital speed */
    arealvel=PI*a*b/T;

    tapo=-T/2.0f;
    tapo=tapo/pow(10,8); /* convert that to readable units */
    
    /* At apocenter the velocity has only tangential component along y-direction */
    v1[0]=0;
    v1[1]=arealvel/ra; /* in kiloparsec/year */
    v1[2]=0;
    for (int j=0;j<3;j++){
        v2[j]=-v1[j];
    }
    r1[0]=ra/2.0f;
    r1[1]=0;
    r1[2]=0;
    for (int j=0;j<3;j++){
        r2[j]=-r1[j];
    }
    printf("the start time for the core masses is %-14.4E\n",tapo);
}

/* Give the initial position (x,y,z) of the two core masses at t=-16.4 of the last apocenter */
/* warning: this needs debugging */
void corepos(r1,r2,Rinit,rp,e,v1,v2)
double r1[];
double r2[];
double Rinit;
double rp;
double e;
double v1[];
double v2[];
{
    double a;
    double phi;
    double relativr[3]; /* relative coordinate vector in r,theta,z*/
    double relativv; /* dtheta/dt of the relative coordinate vector */
    double T; /* period of the kepler orbit */
    double p,b;
    double arealvel;

    /* for the relative coordinate vector, rp is twice the individual one */
    p=rp*(1+e);
    a=p/(1-e*e);
    b=p/sqrt(1-e*e);
    /* Kepler's third law to get period */
    T=sqrt(4*PI*PI*pow(a,3)/(G*pow(1.02201*pow(10,-9.0),2.0)*2*M)); /* this is in calendar year */

    /* Kepler's second law to get orbital speed */
    arealvel=PI*a*b/T;
    relativv=arealvel*2/pow(Rinit,2)*Rinit;

    phi=acos((a*(1-e*e)-Rinit)/(Rinit*e));
    relativr[0]=Rinit;
    relativr[1]=phi;
    relativr[2]=0;

    r1[0]=Rinit/2.0f*cos(phi);
    r1[1]=Rinit/2.0f*sin(phi);
    r1[2]=0;
    v1[0]=-sin(phi)*relativv;
    v1[1]=cos(phi)*relativv;
    v1[2]=0;
         
    for (int i=0;i<3;i++){
        r2[i]=-r1[i];
        v2[i]=-v1[i];
    }
    
}


/* print out documents for initial positions of the disk particles and the core masses in two separate files */
void printdat(r1,r2,v1,v2,x,y,z,vx,vy,vz,fp1,fp2,n)
double r1[];
double r2[];
double v1[];
double v2[];
double x[];
double y[];
double z[];
double vx[];
double vy[];
double vz[];
FILE *fp1;
FILE *fp2;
int n;
{
    
    for (int i=0;i<n;i++)
    {
        fprintf(fp1,"%-14.4E%-14.4E%-14.4E%-14.4E%-14.4E%-14.4E\n",x[i],y[i],z[i],vx[i],vy[i],vz[i]);
    }
    fprintf(fp2,"%-14.4E\n",tapo);
    fprintf(fp2,"%-14.4E%-14.4E%-14.4E%-14.4E%-14.4E%-14.4E\n",r1[0],r1[1],r1[2],v1[0],v1[1],v1[2]);
    fprintf(fp2,"%-14.4E%-14.4E%-14.4E%-14.4E%-14.4E%-14.4E",r2[0],r2[1],r2[2],v2[0],v2[1],v2[2]);
}

void checkinit(r1,r2,x,y,z,fp3,n)
double r1[];
double r2[];
double x[];
double y[];
double z[];
FILE *fp3;
int n;
{
    for (int i=0;i<n;i++)
    {
        fprintf(fp3,"%-14.4E%-14.4E%-14.4E\n",x[i],y[i],z[i]);
    }
    fprintf(fp3,"%-14.4E%-14.4E%-14.4E\n",r1[0],r1[1],r1[2]);
    fprintf(fp3,"%-14.4E%-14.4E%-14.4E\n",r2[0],r2[1],r2[2]);
    /* pick a random point on the outer ring and check its distance from the center mass */
    double dist;
    dist=pow(x[296]-r1[0],2)+pow(y[296]-r1[1],2)+pow(z[296]-r1[2],2);
    dist=sqrt(dist); 
    printf("The distance between an outer radius point and the core mass is %-14.4E\n",dist);
    /* check that the distance between two core masses is Rinit */
    double dist1;
    dist1=pow(r2[0]-r1[0],2)+pow(r2[1]-r1[1],2)+pow(r2[2]-r1[2],2);
    dist1=sqrt(dist1);
    printf("The distance between two core masses is %-14.4E\n",dist1);
}

