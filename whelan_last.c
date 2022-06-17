// This program simulates a three dimensional potential well
// with N (equal to Array_Size) particles randomly
// placed inside the trap. When run
// the program tracks the location of each particle which 
// allows for various sorts of analysis.

//Thomas Francis Whelan May 13, 2019

//THIS CODE MUST BE RUN ON THE CLUSTER with -lm -o command!!!

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define dt 0.0000001 //dt is most accurate a 0.000000001
#define Array_Size 20 //the number of particles in the trap
#define NUMBER_OF_STEPS 100000

double U1(double xi, double yi, int index, int nparticles){
        int i; //finds the potential of a particle in the well
        double r, ro, mu, totalnrgy;
        mu=1;
        r=xi*xi+yi*yi;
        ro=1;
        totalnrgy= (mu*r)/(r+ro);
        return(totalnrgy);
}

double U2(double x[], double y[], int index, int nparticles){
        int i; //finds potential with respect to other particles
        double r, rx,ry,ro, m, mu, xi, xj,yi,yj,F2, totalnrgy;
        mu=1;
        m=1;
        xi=x[index];
        yi=y[index];
        ro=1;
        for(i=0;i<nparticles;i++){
                xj=x[i];
                yi=y[i];
                rx=xj-xi;
                ry=yi-yj;
                r=sqrt(rx*rx+ry*ry);
                F2=0.5*m*ro*ro*(r-ro)*(r-ro)/(r*r);
                totalnrgy+=F2;
        }
        return(totalnrgy);
}
	
	
void Gradient(double x[], double y[], double gradx[], double grady[]){
        int i, nparticles;
        double high, low, x_value,highy,lowy,y_value;
        double U1(),U2();
        double step=0.0001;
	nparticles= Array_Size;
        for(i=0;i<nparticles;i++){ //takes the derivative of our potential
                x_value=x[i];      //to help us calculate RK
                x[i]=x_value+step;
                high= U1(x[i],y[i],i,nparticles) + U2(x,y,i ,nparticles);
                x[i]=x_value-step;
                low= U1(x[i],y[i],i,nparticles) + U2(x,y,i,nparticles);
                x[i]=x_value;
                gradx[i]=(high-low)/(2*step);

                y_value=y[i];     
                y[i]=y_value+step;
                highy= U1(x[i],y[i],i,nparticles) + U2(x,y,i ,nparticles);
                y[i]=y_value-step;
                lowy= U1(x[i],y[i],i,nparticles) + U2(x,y,i ,nparticles);
                y[i]=y_value;
                grady[i]=(highy-lowy)/(2*step);
      }
}

void RK4(double x[], double y[], double vx[], double vy[]){
        int i,N; //uses fourth order runge-kutta integration
	N=Array_Size; //propels simulation forward in time
        double xtemp[N],gradx1[N],gradx2[N],gradx3[N],gradx4[N],x1[N],x2[N],x3[N],x4[N],vx1[N],vx2[N],vx3[N],vx4[N];
        double ytemp[N],grady1[N],grady2[N],grady3[N],grady4[N],y1[N],y2[N],y3[N],y4[N],vy1[N],vy2[N],vy3[N],vy4[N];

        Gradient(x,y,gradx1,grady1);

        for(i=0;i<N;i++){
                x1[i]=vx[i];
                vx1[i]=-gradx1[i];
                xtemp[i]=x[i]+x1[i]*dt/2.;

                y1[i]=vy[i];
                vy1[i]=-grady1[i];
                ytemp[i]=y[i]+y1[i]*dt/2.;
        }
        Gradient(xtemp,ytemp,gradx2,grady2);

        for(i=0;i<N;i++){
                x2[i]=vx[i]+vx1[i]*dt/2.;
                vx2[i]=-gradx2[i];
                xtemp[i]=x[i]+x2[i]*dt/2.;

                y2[i]=vy[i]+vy1[i]*dt/2.;
                vy2[i]=-grady2[i];
                ytemp[i]=y[i]+y2[i]*dt/2.;
        }
        Gradient(xtemp,ytemp,gradx3,grady3);

        for(i=0;i<N;i++){
                x3[i]=vx[i]+vx2[i]*dt/2.;
                vx3[i]=-gradx3[i];
                xtemp[i]=x[i]+x3[i]*dt;

                y3[i]=vy[i]+vy2[i]*dt/2.;
                vy3[i]=-grady3[i];
                ytemp[i]=y[i]+y3[i]*dt;
        }
        Gradient(xtemp,ytemp,gradx4,grady4);

        for(i=0;i<N;i++){
                x4[i]=vx[i]+vx3[i]*dt;
                vx4[i]=-gradx4[i];
                x[i]+=(dt/6.)*(x1[i]+2.*x2[i]+2.*x3[i]+x4[i]);
                vx[i]+=(dt/6.)*(vx1[i]+2.*vx2[i]+2.*vx3[i]+vx4[i]);

                y4[i]=vy[i]+vy3[i]*dt;
                vy4[i]=-grady4[i];
                y[i]+=(dt/6.)*(y1[i]+2.*y2[i]+2.*y3[i]+y4[i]);
                vy[i]+=(dt/6.)*(vy1[i]+2.*vy2[i]+2.*vy3[i]+vy4[i]);
        }
}


double KE(double vx,  double vy){
	double nrgy,m; //finds KE for one particle
	m=1;
	nrgy= 0.5*m*((vx*vx) +(vy*vy));
	return(nrgy);
}

double PE(double x[],  double y[], int index,int nparticles){
        int i; //finds PE for one particle
        double nrgy, U1(), U2();
        nrgy=U1(x[index],y[index],i,nparticles) + U2(x,y,index,nparticles);
	return(nrgy);
}

double TE(double x[], double y[], double vx[], double vy[], int nparticles){
	int i; //finds total energy for system
	double adder,total;
	for(i=0,total=0;i<nparticles;i++){
		total+= KE(vx[i],vy[i]) + PE(x,y,i,nparticles);
	}
	return(total);
}

int checker(double x[], double y[], double vx[], double vy[]){
	int i,escaped,nparticles; //checks if particle has left trap
	double TE, r, KE(), PE();
	nparticles= Array_Size;
	for(i=0,escaped=0;i<nparticles;i++){
		r=sqrt(x[i]*x[i]+y[i]*y[i]);
		TE= KE(vx[i],vy[i]) + PE(x,y,i,nparticles);
		if ((r>50) && (TE>0)){
			escaped+=1;
		}
	}	
	return(escaped);
}

int main(){
        int i,n,j,nparticles,range,steps,goners,checker();
        int vbin[201]={0};
	double val,r(),theta(),TE(),startnrgy,change,endnrgy,v;
        double summer, numb, avgx,help,std_devx,avgy,std_devy;
	nparticles= Array_Size; 
        double x[Array_Size];
	double y[Array_Size];
        double vx[Array_Size]={0};
	double vy[Array_Size]={0};
	int vxbin[121]={0};
	int vybin[121]={0};
	double finalv[Array_Size]={};
        FILE *fp;
	fp=fopen("start.out","w");
	FILE *fo;
	fo=fopen("final.out","w");
	FILE *ft;
        ft=fopen("dist.out","w");
	
	steps=NUMBER_OF_STEPS;
 	srandom(time(0));
        for(j=0;j<nparticles;j++){ //randomizes initial position
        	int i,big,range,ff;
        	double data,spread,decimiler,radius;
        	big= 10000000;
        	data=2*M_PI*(random()% big)/big;
		range=6;
		spread=(random() % range);
        	ff=10000000;
        	decimiler= (random() % ff);
        	radius= spread+decimiler/10000000;

		x[j]= radius*cos(data);
		y[j]= radius*sin(data);
        }
	
	startnrgy= TE(x,y,vx,vy,nparticles); //initial nrgy of system
	for(i=0;i<steps;i++){
              RK4(x,y,vx,vy); //executes RK4 for our arrays
	}	
	endnrgy=TE(x,y,vx,vy,nparticles); //final nrgy of system
	change= ((endnrgy-startnrgy)/startnrgy)*100; //%change in energy
	goners=checker(x,y,vx,vy); //checks how many particles left trap
	for(i=0;i<nparticles;i++){
		finalv[i]=sqrt(vx[i]*vx[i]+vy[i]*vy[i]);
	} 
	for(i=0,summer=0;i<nparticles;i++){
                numb= vx[i];
                summer+= numb;
        }
        avgx= summer/nparticles; //final velocity mean in x direction
        for(i=0,help=0;i<nparticles;i++){
                help+= ((vx[i] - avgx)*(vx[i] - avgx));
        }
        std_devx=sqrt(help/nparticles); //final velocity std dev in x direction

        for(i=0,summer=0,numb=0;i<nparticles;i++){
                numb= vy[i];
                summer+= numb;
        }
        avgy= summer/nparticles; //final velocity mean in y direction
        for(i=0,help=0;i<nparticles;i++){
                help+= ((vy[i] - avgy)*(vy[i] - avgy));
        }
        std_devy=sqrt(help/nparticles); //final velocity std dev in y direction
	for(i=0,v=0;i<nparticles;i++){
                v=vx[i];
		v+=30.5;
		v*=2;
                if(v>120);
                else if(v<0);    //bins our final x velocities
                else vxbin[(int)v]++ ;
	}

	for(i=0,v=0;i<nparticles;i++){
                v=vy[i];
		v+= 30.5;
		v*=2;	
                if(v>120);
                else if(v<0);    //bins our final y velocities
                else vybin[(int)v]++ ;
	}
	
	for(i=0;i<121;i++){
		fprintf(ft,"%d\t%d\t%d\n",i-61,vxbin[i],vybin[i]);
	}
	
	printf("%d\t%f\t%f\t%f\t%f\n",i,x[8],y[8],vx[8],vy[8]);
	
	printf("Mean X Velocity is %f with std dev %f\n",avgx,std_devx);
	printf("Mean Y Velocity is %f with std dev %f\n",avgy,std_devy);
	printf("There are %d particles that escaped the trap\n", goners);
	printf("The percent chance in energy is %f", change);
}
