#include<iostream>
#include<fstream>
#include<cmath>

using namespace std;

const double Re=50.0 ; // length of container= 1unit, rho= 10^3 SI units, viscosity= 8.9x 10^-4
const double length1= 12.0;// so that reynolds number is equal to 100 and peclet number less than two
const double length2= 1.0;
const double dt= 0.0001;// from peclet no. less than 1
const int N= 301;
const int M= 51;
const double rho= 1.0;
double tolerance=pow(10,-6);
const double NT=300000;
// defining variables
double u_infinity= 1.0;
double p[N+1][M+1],p_error[N+1][M+1]={0.0};
double stream_function[N][M], vorticity[N][M]={0.0};
double u[N+1][M+1], un[N+1][M+1],u_n_1[N+1][M+1]={0.0};
double v[N+1][M+1], vn[N][M+1],v_n_1[N+1][M+1]={0.0};
double l[N+1],b[N+1],c[N+1],d[N+1],x[N+1];

double a(int ix, int iy,int z);
double a_w(int ix, int iy,int z);
double a_e(int ix, int iy,int z);
double a_s(int ix, int iy,int z);
double a_n(int ix, int iy,int z);
double p_e(int ix, int iy);
double p_s(int ix, int iy);
double p_n(int ix, int iy);
double p_w(int ix, int iy);
double boundary(int ix, int iy);
double H(int ix, int iy, int z);
double bi(int ix, int iy);
double dx= length1/(N-1);
double dy= length2/(M-1);
//double viscosity= (rho*u_infinity*length2)/rey ;
//double D1= viscosity/dx;
//double D2= viscosity/dy;
double R1=dy/dx;
double R2=dx/dy;
void gs_u(void);
void gs_v(void);
void gs_p(void);
void TDMAm(void);
void TDMAn(void);
void TDMApm(void);


int main(){



    // WRITING VALUES TO TOLERANCE FILE

    ofstream fout15;
    fout15.open("tolerance.dat");

    /*for(int i=0;i<M;i++){
      u[0][i]=u_infinity;// this value from calculation from Reynold's number =100.0
    }*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //                       SIMPLE algorithm
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                   // copying values of U-velocity
        for(int ix=0;ix<N;ix++){
            for(int iy=0;iy<=M;iy++){
                u[ix][iy]=0.0;
            }
       }

       // copying values of V-velocity

       for(int ix=0;ix<=N;ix++){
        for(int iy=0;iy<M;iy++){
            v[ix][iy]=0.0;
        }
       }
        for(int ix=0;ix<=N;ix++){
        for(int iy=0;iy<=M;iy++){
            p[ix][iy]=0.0;
        }
       }

    for(int it=NT;it>0;it--){

        cout<<"Iteration number: "<<NT- it<<endl;


          //Boundary Conditions
    for (int i=1;i<M;i++){
        u[0][i]=1.0;
        u[N-1][i]=u[N-2][i];}
        for (int i=1;i<N;i++)
        {v[i][0]=0.0;
        v[i][M-1]=0.0;}
    for (int i=0;i<=M-1;i++){
        v[0][i]=-v[1][i];
        v[N][i]=-v[N-1][i];}
    for (int i=0;i<=N-1;i++){
        u[i][0]=-u[i][1];
        u[i][M]=-u[i][M-1];

                }

        for (int i=0;i<=N-1;i++)
        //cout<<"hi"<<u[i][M]<<endl;
         for(int ix=0;ix<=M;ix++){
            p[N][ix]=p[N-1][ix];
        }
        for(int ix=0;ix<=N;ix++){
            p[ix][M]=p[ix][M-1];
        }
         for(int ix=0;ix<=M;ix++){
            p[0][ix]=p[1][ix];
        }
        for(int ix=0;ix<=N;ix++){
            p[ix][0]=p[ix][1];
        }

        p[0][0]=0.0;



            // copying values of U-velocity
        for(int ix=0;ix<N;ix++){
            for(int iy=0;iy<=M;iy++){
                u_n_1[ix][iy]=u[ix][iy];
            }
       }

       // copying values of V-velocity

       for(int ix=0;ix<=N;ix++){
        for(int iy=0;iy<M;iy++){
            v_n_1[ix][iy]=v[ix][iy];
        }
       }



        // solving the system of linear equations for u and v momentum equations and pressure correction equation
        gs_u();
        gs_v();

        for(int ix=0;ix<N;ix++){
            for(int iy=0;iy<M;iy++){
                p_error[ix][iy]=0.0;
            }
        }


        gs_p();

          // updating pressure with under-relaxation factor

        for(int ix=1;ix<N;ix++){
            for(int iy=1;iy<M;iy++){

                p[ix][iy]= p[ix][iy]+ 0.1*p_error[ix][iy];// here '0.5' is the under relaxation  factor

               }
        }

        // update velocity values using correction equations

        //updating u velocity





        for(int ix=N-2;ix>0;ix--){
            for(int iy=M-1;iy>0;iy--){


                    u[ix][iy]= u[ix][iy]+ (dy*1.0/a(ix,iy,1))*(p_error[ix][iy]-p_error[ix+1][iy]);//updating velocity values without using under relaxation factor

                    }
            }

        for(int ix=N-1;ix>0;ix--){
            for(int iy=M-2;iy>0;iy--){

            v[ix][iy]= v[ix][iy]+ (dx*1.0/a(ix,iy,2))*(p_error[ix][iy]-p_error[ix][iy+1]);//updating velocity values without using under relaxation factor
            }}
   /*    for(int ix=1;ix<N-1;ix++){
            for(int iy =1;iy<M-1;iy++){

                    u[ix][iy]=1.3*u[ix][iy]- 0.3*u_n_1[ix][iy];// updating velocity values using under relaxation factor

                   // v[ix][iy]=0.5*v[ix][iy]+ 0.5*v_n_1[ix][iy];// updating velocity values using under relaxation factor
            }
        }*/




        cout<<u[N/2][M/2]<<endl;
        cout<<u[N][M/2]<<endl;


        //checking convergence

        double maximum1=0.0;
        double maximum2=0.0;
        double maximum=0.0;

        for(int ix=0;ix<N-1;ix++){
            for(int iy=0;iy<M;iy++){

                if(abs(u[ix][iy]-u_n_1[ix][iy])>maximum1){
                    maximum1= abs(u[ix][iy]-u_n_1[ix][iy]);
                   // cout<<ix<<","<<iy<<endl;
                }
            }
        }
        for(int ix=0;ix<N;ix++){
            for(int iy=0;iy<M-1;iy++){

           if(abs(v[ix][iy]-v_n_1[ix][iy])>maximum1){
                    maximum1= abs(v[ix][iy]-v_n_1[ix][iy]);
                     }
            }
        }

        if(maximum1<maximum2){
            maximum = maximum2;
        }
        else if(maximum2<maximum1){
            maximum = maximum1;
        }
        if(maximum<tolerance)
            break;


        cout<<"Value of error:"<<maximum<<endl;

        fout15<<NT-it<<" "<<maximum<<endl;


}

 cout<<"Writing values to file."<<endl;
            //write values to a file

         ofstream fout12;
         ofstream fout13;
         ofstream fout14;
        ofstream fout10;

            fout12.open("pressure.dat");
            fout13.open("mid_line_x_Velocity.dat");
            fout14.open("mid_line_y_velocity.dat");



            for(int ix=0;ix<N;ix++){
                for(int iy=0;iy<M;iy++){
                    fout12<<ix*dx<<","<<iy*dx<<","<<p[ix][iy]<<endl;// pressure field
                }
            }
           int m;
            for(int i=0;i<M;i++){
                    m=50;
            fout13<<0.5*(u[10][i]+u[10][i+1])<<" ";
            while (m<=N)
           {fout13<<0.5*(u[m][i]+u[m][i+1])<<" ";
                m=m+50;}
            fout13<<(dy*i)/length2<<endl;}
            for(int i=0;i<M;i++){
                //fout13<<u[10][i]/u_infinity<<" "<<u[25][i]/u_infinity<<" "<<u[50][i]/u_infinity<<" "<<u[75][i]/u_infinity<<" "<<u[100][i]/u_infinity<<" "<<u[125][i]/u_infinity<<" "<<u[150][i]/u_infinity<<" "<<u[175][i]/u_infinity<<" "<<u[200][i]/u_infinity<<" "<<u[225][i]/u_infinity<<" "<<u[250][i]/u_infinity<<" "<<u[275][i]/u_infinity<<" "<<u[300][i]/u_infinity<<" "<<(dy*i)/length2<<endl;// mid line u- velocity
                cout<<0.5*(u[i][N-1]+u[i+1][N-1])<<endl<<endl;
            }

            for(int i=0;i<N;i++){

               fout14<<dx*i/length1<<" "<<v[i][65]/u_infinity<<endl;// mid line v- velocity

            }
    fout10.open("u Velocities and pressure.dat");

               for(int i = 0; i < N; i++){
			for(int j = 0; j < M; j++)
			{
				fout10<<i*dx<<" "<<j*dy<<" "<<0.5*(u[i][j]+u[i][j+1])<<" "<<0.5*(v[i][j]+v[i+1][j])<<" "<<0.25*(p[i][j]+p[i+1][j]+p[i+1][j+1]+p[i][j+1])<<endl;
			}
		}
    fout10.close();


            fout12.close();
            fout13.close();
            fout14.close();


fout15.close();

}

// momentum equation coefficients

double a_n(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
     ans= dx*0.25*(v_n_1[ix][iy]+v_n_1[ix+1][iy])-R2/Re;
     else if(z==2)
     ans= dx*0.25*(v_n_1[ix][iy]+v_n_1[ix][iy+1])-R2/Re;

    return ans;
}

double a_s(int ix, int iy,int z){

    double ans=0.0;

        if(z==1)
        ans= -1.0*dx*0.25*(v_n_1[ix][iy-1]+v_n_1[ix+1][iy-1])-R2/Re;
        else if(z==2)
        ans= -1.0*dx*0.25*(v_n_1[ix][iy-1]+v_n_1[ix][iy+1])-R2/Re;

    return ans;
}


double a_e(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
      ans= dy*0.25*(u_n_1[ix][iy]+u_n_1[ix+1][iy])-R1/Re;
      else if(z==2)
      ans= dy*0.25*(u_n_1[ix][iy]+u_n_1[ix][iy+1])-R1/Re;

    return ans;
}

double a_w(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
    ans= -1.0*dy*0.25*(u_n_1[ix][iy]+u_n_1[ix-1][iy])-R1/Re;
    else if(z==2)
    ans= -1.0*dy*0.25*(u_n_1[ix-1][iy]+u_n_1[ix-1][iy+1])-R1/Re;


    return ans;
}


double a(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
     //ans= dx*dy/dt +(a_w(ix,iy,1)+a_e(ix,iy,1)+a_n(ix,iy,1)+a_s(ix,iy,1))+4.0*(R1/Re+R2/Re);
     ans= (a_w(ix,iy,1)+a_e(ix,iy,1)+a_n(ix,iy,1)+a_s(ix,iy,1))+4.0*(R1/Re+R2/Re);
      //ans= (a_w(ix,iy,1)+a_e(ix,iy,1)+a_n(ix,iy,1)+a_s(ix,iy,1));
     else if(z==2)
     //ans= dx*dy/dt +(a_w(ix,iy,2)+a_e(ix,iy,2)+a_n(ix,iy,2)+a_s(ix,iy,2))+4.0*(R1/Re+R2/Re);
    ans= (a_w(ix,iy,2)+a_e(ix,iy,2)+a_n(ix,iy,2)+a_s(ix,iy,2))+4.0*(R1/Re+R2/Re);
    //ans= (a_w(ix,iy,2)+a_e(ix,iy,2)+a_n(ix,iy,2)+a_s(ix,iy,2));
    return ans;
}

double H(int ix, int iy, int z){
    double ans =  0.0;
    if(z==1)
    ans= (p[ix][iy]-p[ix+1][iy])*dy+u_n_1[ix][iy]*dx*dy/dt;

    else if(z==2)
    ans= (p[ix][iy]-p[ix][iy+1])*dx+v_n_1[ix][iy]*dx*dy/dt;
    return ans;
}
// pressure correction equation coefficients

double p_e(int ix, int iy){

    double answer=0.0;

    answer= -1.0*dy*dy*1.0/a(ix,iy,1);

    return answer;

}

double p_n(int ix, int iy){

    double answer=0.0;


    answer= -1.0*dx*dx*1.0/a(ix,iy,2);


    return answer;

}

double p_w(int ix, int iy){

    double answer=0.0;

    answer= -1.0*dy*dy*1.0/a(ix-1,iy,1);

    return answer;

}

double p_s(int ix, int iy){

    double answer=0.0;

    answer= -1.0*dx*dx*1.0/a(ix,iy-1,2);


    return answer;

}

double bi(int ix, int iy){
    double answer=0.0;
    answer= v_n_1[ix][iy-1]*dx-v_n_1[ix][iy]*dx+u_n_1[ix-1][iy]*dy-u_n_1[ix][iy]*dy;
    return answer;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void gs_u(void){
double re;
for(int iy=1;iy<M;iy++){

    for(int ix=1;ix<N-1;ix++){

      if(ix==1)
        {
            d[ix]=(p[ix][iy]-p[ix+1][iy])*dy-u[ix][iy+1]*(a_n(ix,iy,1))-u[ix][iy-1]*(a_s(ix,iy,1))-u[ix-1][iy]*a_w(ix,iy,1);
            l[ix]=0; c[ix]=a_e(ix,iy,1); b[ix]=a(ix,iy,1);
            continue;
        }
        if(ix==N-2)
        {
            d[ix]=(p[ix][iy]-p[ix+1][iy])*dy -u[ix][iy+1]*(a_n(ix,iy,1))-u[ix][iy-1]*(a_s(ix,iy,1))-u[ix+1][iy]*a_e(ix,iy,1);
            l[ix]=a_w(ix,iy,1); b[ix]=a(ix,iy,1);
            continue;
        }
        l[ix]=a_w(ix,iy,1);
        c[ix]=a_e(ix,iy,1);
        b[ix]=a(ix,iy,1);
        d[ix]=(p[ix][iy]-p[ix+1][iy])*dy-u[ix][iy+1]*(a_n(ix,iy,1))-u[ix][iy-1]*(a_s(ix,iy,1));
                  }
    TDMAn();
    for(int ix=1;ix<N-1;ix++)
        u[ix][iy]=x[ix];

}
/*
for(int ix=1;ix<N-1;ix++){

    for(int iy=1;iy<M;iy++){
         if(iy==1)
        {
            d[iy]=(p[ix][iy]-p[ix+1][iy])*dy - u[ix+1][iy]*(a_e(ix,iy,1))-u[ix-1][iy]*(a_w(ix,iy,1))-u[ix][iy-1]*a_s(ix,iy,1);
            l[iy]=0; c[iy]=a_n(ix,iy,1); b[iy]=a(ix,iy,1);
            continue;
        }
        if(iy==M-1)
        {
            d[iy]=(p[ix][iy]-p[ix+1][iy])*dy- u[ix+1][iy]*(a_e(ix,iy,1))-u[ix-1][iy]*(a_w(ix,iy,1))-u[ix][iy+1]*a_n(ix,iy,1);
            l[iy]=a_s(ix,iy,1); b[iy]=a(ix,iy,1);
            continue;
        }
        l[iy]=a_s(ix,iy,1);
        c[iy]=a_n(ix,iy,1);
        b[iy]=a(ix,iy,1);
        d[iy]=(p[ix][iy]-p[ix+1][iy])*dy- u[ix+1][iy]*(a_e(ix,iy,1))-u[ix-1][iy]*(a_w(ix,iy,1));
                  }
    TDMAm();
    for(int iy=1;iy<=M-1;iy++)
        {u[ix][iy]=x[iy];
        //cout<<u[ix][iy]<<endl;
        }
}*/

for(int iy=1;iy<M-1;iy++){

    for(int ix=1;ix<N-1;ix++){

   if(u[ix][iy]!=u[ix][iy])
            cout<<ix<<"..."<<iy<<"....././"<<u_n_1[ix][iy]<<endl;
    }}
}

void gs_p(void){
double qw;
for(int iy=1;iy<M;iy++){

   for(int ix=1;ix<N;ix++){
     if(ix==1)
        {
            d[ix]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy -p_error[ix][iy-1]*p_s(ix,iy) - p_error[ix][iy+1]*p_n(ix,iy)-p_error[ix-1][iy]*p_w(ix,iy);
            l[ix]=0; c[ix]=p_e(ix,iy); b[ix]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
            continue;
        }
        if(ix==N-1)
        {
            d[ix]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy-p_error[ix][iy-1]*p_s(ix,iy) - p_error[ix][iy+1]*p_n(ix,iy)-p_error[ix+1][iy]*p_e(ix,iy);
            l[ix]=p_w(ix,iy); b[ix]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
            continue;
        }
        l[ix]=p_w(ix,iy);
        c[ix]=p_e(ix,iy);
        b[ix]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
        d[ix]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy-p_error[ix][iy-1]*p_s(ix,iy) - p_error[ix][iy+1]*p_n(ix,iy);
                  }
    TDMApm();
    for(int ix=1;ix<N;ix++)
        p_error[ix][iy]=x[ix];
}/*
for(int ix=1;ix<N;ix++){

    for(int iy=1;iy<M;iy++){
         if(iy==1)
        {
            d[iy]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy -p_error[ix-1][iy]*p_w(ix,iy) - p_error[ix+1][iy]*p_e(ix,iy)-p_error[ix][iy-1]*p_s(ix,iy);
            l[iy]=0; c[iy]=p_n(ix,iy); b[iy]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
            continue;
        }
        if(iy==M-1)
        {
            d[iy]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy -p_error[ix-1][iy]*p_w(ix,iy) - p_error[ix+1][iy]*p_e(ix,iy)-p_error[ix][iy+1]*p_n(ix,iy);
            l[iy]=p_s(ix,iy); b[iy]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
            continue;
        }
        l[iy]=p_s(ix,iy);
        c[iy]=p_n(ix,iy);
        b[iy]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
        d[iy]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy -p_error[ix-1][iy]*p_w(ix,iy) - p_error[ix+1][iy]*p_e(ix,iy);
                  }
    TDMAm();
    for(int iy=1;iy<M;iy++)
        p_error[ix][iy]=x[iy];
}*/
}


void gs_v(void){

for(int ix=1;ix<N;ix++){

    for(int iy=1;iy<M-1;iy++){

        v[ix][iy]= (dx*1.0*((p[ix][iy]-p[ix][iy+1]))-1.0*(v[ix+1][iy]*(a_e(ix,iy,2)) + v[ix][iy+1]*(a_n(ix,iy,2)) + v[ix-1][iy]*(a_w(ix,iy,2)) + v[ix][iy-1]*(a_s(ix,iy,2))))/(a(ix,iy,2));

    }
}

}
void TDMAn(void){
 int i;
 double t;
    for (i=2;i<=N-2;i++)
    {
        l[i]=l[i]/(b[i-1]);
        b[i]=b[i]-l[i]*c[i-1];
        d[i]=d[i]-l[i]*d[i-1];
    }
    int n=N-2;
    //double x[n];
    x[n]=d[n]/b[n];
    for (i=n-1;i>=1;i--)
    {
        x[i]=(d[i]-c[i]*x[i+1])/b[i];
        t=x[i];
    }
}
void TDMAm(void){
 int i;
    for (i=2;i<=M-1;i++)
    {
        l[i]=l[i]/(b[i-1]);
        b[i]=b[i]-l[i]*c[i-1];
        d[i]=d[i]-l[i]*d[i-1];
    }
    int n=M-1;
    x[n]=d[n]/b[n];
    for (i=n-1;i>=1;i--)
    {
        x[i]=(d[i]-c[i]*x[i+1])/b[i];
    }
}
void TDMApm(void){
 int i;
    for (i=2;i<=N-1;i++)
    {
        l[i]=l[i]/(b[i-1]);
        b[i]=b[i]-l[i]*c[i-1];
        d[i]=d[i]-l[i]*d[i-1];
    }
    int n=N-1;
    x[n]=d[n]/b[n];
    for (i=n-1;i>=1;i--)
    {
        x[i]=(d[i]-c[i]*x[i+1])/b[i];
    }
}
