 
 
/*
25 August 2026
Problem: Double diffusive Natural convection in an enclosed cavity.
The left wall of the cavity are maintained at uniform high temp. and low concentration.
The right wall of the cavity are maintained at uniform low temp. and high concentration.
The Right wall of the cavity is at ambient and lower concentration.
Top wall is moving with U and thermally insulated.  
Ra=1000-100000,
Pr=0.054,
Le=2-10,
N (buoyancy Ration)=-2-2
Ha=0-50
*/

#include<iostream>
#include<stdio.h>
#include<fstream>
#include<cmath>
//#include<conio.h>
//# define H 1
# define nx 200
# define ny 200
//# define uo 0.0
# define Ra 10000   // Rayleigh number i changed to 10000 from 1000
# define Ha 50    // Hartmann number
# define Le 2     // Lewis number
# define N -2    // Buoyancy Ratio - testing aiding
# define Pr 0.054  // Prandtl number - i changed 
# define Ma 0.1   // Mach number
# define AP 0  // Velocity parameter (A=U*L/alpha) - i changed - working on natural convection solver
#define ci 1.0/sqrt(3)  // lattice speed constant
  void flow_streaming();
  void heat_streaming();
  void heat_boundary_conditions();
  void mass_streaming();
  void mass_boundary_conditions();
using namespace std;
 double f[800][800][9],feq[800][800][9],rho[800][800],u[800][800],v[800][800],w[9],cx[9],cy[9];
 double streamfun[800][800], rhoavg, rhom, old_T[800][800], error_T, error_C, old_u[800][800],old_C[800][800];
 double g[800][800][9], geq[800][800][9],T[800][800],Tw=1.0,Tc=0.0, Tref, Ch=1.0,Cl=0.0,Cref,force[800][800][9];
 double h[800][800][9], heq[800][800][9],C[800][800];
 
    int i=nx+2,j=ny+2,k=9;
    double vor[800][800];
  
  /////////////////////////////////////////
  const double  kinevisc=sqrt((Ma*Ma*ny*ny*Pr*ci*ci)/Ra);
  
  const double omega=1.0/(3.0*kinevisc+0.5);
  ////////////////////////////////////////////
 
  const double Tdiffusivity=kinevisc/Pr;   // Thermal diffusivitty
  
  const double omegaT=1.0/(3.0*Tdiffusivity+0.5);
  /////////////////////////////////////////////
  const double Diffusivity=kinevisc/(Pr*Le);  // Mass Diffusivity 
  
  const double omegaC=1.0/(3.0*Diffusivity+0.5); 
  
  const double omega1=(1.0-omega);  
   
  const double omegaT1=(1.0-omegaT);
    
  const double omegaC1=(1.0-omegaC);
  //*********************************************
  const double gbetaT=(Ra*kinevisc*kinevisc)/((ny*ny*ny)*Pr);
  
  const double betaT=gbetaT/10.0;
  
  const double betaC=N*betaT;
  
 const double A=(Ha*Ha*kinevisc)/(ny*ny); // i changed - added magnetic field effect - MHD force constant
  // charactrisc velocity***********************
  //const double Vref=(kinevisc/ny)*sqrt(Ra/Pr);
//**********************************************
// Lid Velocity/////
 const double ulid=(AP*Tdiffusivity)/ny;
   
int main ()
{
    int t,iterate,maxit=80000000;
  
  // double Nur1[j], Nul1[j]; 
    double error, error_u,error_v, old_v[ny+1][nx+1], max_error=1.e-8;
    double X,Y,p,q;   
   
    double avgNun, avgNus, Nun, Nus, Nun1, Nus1, avgNu; 
    double Nul, Nur, Nul1[1000], Nur1[1000], avgNul,avgNur, Nuw,Nue,Nub,avgNub,Nue1, avgNuw,avgNue, Shl,Shb,avgShb,avgShl,avgShr;
    
   // CONSTANTS
    // weight functions
    w[0]=(4.0/9.0);   // weight function for point 1
    w[1]=(1.0/9.0);   // weight function for point 2
    w[2]=(1.0/9.0);   // weight function for point 3
    w[3]=(1.0/9.0);   // weight function for point 4
    w[4]=(1.0/9.0);   // weight function for point 5
    w[5]=(1.0/36.0);   // weight function for point 6 
    w[6]=(1.0/36.0);   // weight function for point 7 
    w[7]=(1.0/36.0);   // weight function for point 8 
    w[8]=(1.0/36.0);   // weight function for point 9 
    
    // velocity for distribution functions-x
    cx[0]=0.0;   
    cx[1]=1.0;  
    cx[2]=0.0;   
    cx[3]=-1.0;   
    cx[4]=0.0;   
    cx[5]=1.0;   
    cx[6]=-1.0;   
    cx[7]=-1.0;   
    cx[8]=1.0;  
    
     // velocity for distribution functions-y
    cy[0]=0.0;   
    cy[1]=0.0;  
    cy[2]=1.0;   
    cy[3]=0.0;   
    cy[4]=-1.0;   
    cy[5]=1.0;   
    cy[6]=1.0;   
    cy[7]=-1.0;   
    cy[8]=-1.0;  

   //////////////////Mean Values///////////////////
   Tref=(Tw+Tc)/2.0;
   Cref=(Ch+Cl)/2.0;
    
    
    // initialization
    

    for (j=0; j<=ny; j++)
    {
        for (i=0; i<=nx; i++)
        {
            
            rho[j][i]=2.7;
            u[j][i]=0.0;
            v[j][i]=0.0;
           
           
    
    for (k=0; k<=8; k++)
    {
        
                      
    feq[j][i][k]=(rho[j][i]*w[k])*(1.0 +3.0*(u[j][i]*cx[k]+v[j][i]*cy[k])
                                   + (4.50*(u[j][i]*cx[k]+v[j][i]*cy[k])*(u[j][i]*cx[k]+v[j][i]*cy[k]))
                                   - 1.50*(u[j][i]*u[j][i]+v[j][i]*v[j][i])); 
                                   
 
          
          
           // set f as feq
             f[j][i][k]=feq[j][i][k]; 
        
            }
        
        }
    }
    // INITIALIZATION COMPLETED.
    
    // main iteration loop starts
    
    t=0;
    iterate:
            t=t+1;
      {  // main iteration loop starts
            
  
         for (j=0; j<=ny; j++)
            {
                for (i=0; i<=nx; i++)
                {
                    old_u[j][i]=u[j][i];
                    old_v[j][i]=v[j][i];
                    old_T[j][i]=T[j][i];
                    old_C[j][i]=C[j][i];
                    
                }
            }
            
            
       
     
            // collision step
            
     for (j=0; j<=ny; j++)
      {
        for (i=0; i<=nx; i++)
        {  
          
          for (k=0; k<=8; k++)
            {
    feq[j][i][k]=(rho[j][i]*w[k])*(1.0 +3.0*(u[j][i]*cx[k]+v[j][i]*cy[k])
                   + (4.50*(u[j][i]*cx[k]+v[j][i]*cy[k])*(u[j][i]*cx[k]+v[j][i]*cy[k]))
                   - 1.50*(u[j][i]*u[j][i]+v[j][i]*v[j][i]));
    // force evaluation for bouyancy driven convection
   // force[j][i][k]=(3.0*w[k]*gbeta*(T[j][i]-Tref)*cy[k]*rho[j][i]); 
    //  force[j][i][k]=(3.0*w[k]*gbetaT*cy[k]*rho[j][i])*((T[j][i]-Tref)+(C[j][i]-Cref)*N);    
     force[j][i][k]=(3.0*w[k]*cy[k]*rho[j][i])*(gbetaT*((T[j][i]-Tref)+(C[j][i]-Cref)*N) - A*v[j][i]);

 //----------------------------------------------------------- 
  if (i==0){
		    force[j][i][k]=0.0;
	    }
  if (i==nx) {
		force[j][i][k]=0.0;
	     }


   if (j==0)  {
		 force[j][i][k]=0.0;
	      }
  if (j==ny) {
               force[j][i][k]=0.0;
	     }
   //---------------------------------------------------------       
        // collision equation : RHS of LBE
        
        f[j][i][k]=(omega1*f[j][i][k])+(omega*feq[j][i][k])+force[j][i][k];
        
            }
            
        }
      }
      
       
      
      // streaming step: LHS of LBE

  flow_streaming();
 
 //----------------------------------------------------------

      // Flow boundary conditions
  
    for (j=0; j<=ny; j++)
    {
        
        // west wall: bounceback bounadry conditions
             f[j][0][1]=f[j][0][3];        
              f[j][0][5]=f[j][0][7];         
               f[j][0][8]=f[j][0][6];        
      
         
          // east wall: bounceback bounadry conditions:Open ended conditions
               f[j][nx][3]=f[j][nx][1];       
                  f[j][nx][7]=f[j][nx][5];        
                    f[j][nx][6]=f[j][nx][8];        
           //------------------------------------------------------------
      }
      
      
    for (i=0; i<=nx; i++)
    {
       
	// North wall
       // rho[ny][i]=f[ny][i][0]+f[ny][i][1]+f[ny][i][3]+2.0*(f[ny][i][2]+f[ny][i][5]+f[ny][i][6]);
	  f[ny][i][4]=f[ny][i][2];  
	   f[ny][i][7]=f[ny][i][5];
	    f[ny][i][8]=f[ny][i][6];  

        // south wall: bounceback bounadry conditions
        
          f[0][i][2]=f[0][i][4];  
            f[0][i][6]=f[0][i][8];      
               f[0][i][5]=f[0][i][7];  
          
    }
     

    // post streaming density and velocity calculations
    
      for (j=0; j<=ny; j++)
      {
        for (i=0; i<=nx; i++)
        {
            rho[j][i]=0.0;
            u[j][i]=0.0;
            v[j][i]=0.0;
          
            for (k=0; k<=8; k++)
            {
                rho[j][i]+=f[j][i][k];
                u[j][i] +=(cx[k]*f[j][i][k]);
                v[j][i] +=(cy[k]*f[j][i][k]);
            
            }
            
            u[j][i]/=rho[j][i];
            v[j][i]/=rho[j][i];
         
        }
        
      }
      
 // scalar field----------------------------------------------------------------
 //****************************************************************************
 
   for (j=0; j<=ny; j++)
    {
        for (i=0; i<=nx; i++)
        {
            T[j][i]=0.0;
            for (k=0; k<=8; k++)
            {
                T[j][i]+=(g[j][i][k]);
                
            }
            
        }
    }     
   // collision: scalar 
  

for (j=0; j<=ny; j++)
{
    for (i=0; i<=nx; i++)
    {
        for (k=0; k<=8; k++)
            {
                  
          //equlibrium distribution fun
            
            geq[j][i][k]=(w[k]*T[j][i])*(1.0+ 3.0*(cx[k]*u[j][i]+cy[k]*v[j][i]));           
            
           
            // collision equation
            g[j][i][k]=(omegaT1*g[j][i][k]+omegaT*geq[j][i][k]);
            
           }
     }
 }    
 
 // streaming
 
  heat_streaming();
 
  // tempearture boundary conditionsS
  
 heat_boundary_conditions();
 
        // POST STREAMING CALCULATION OF TEMP
    
    for (j=0; j<=ny; j++)
    {
        for (i=0; i<=nx; i++)
        {
            T[j][i]=0.0;
            for (k=0; k<=8; k++)
            {
                T[j][i]+=(g[j][i][k]);
                
            }
            
        }
    }     
 //***************************************************************************
  // Concentration field----------------------------------------------------------------
 
   for (j=0; j<=ny; j++)
    {
        for (i=0; i<=nx; i++)
        {
            C[j][i]=0.0;
            for (k=0; k<=8; k++)
            {
                C[j][i]+=(h[j][i][k]);
                
            }
            
        }
    }     
   // collision: scalar 
  

for (j=0; j<=ny; j++)
{
    for (i=0; i<=nx; i++)
    {
        for (k=0; k<=8; k++)
            {
                  
          //equlibrium distribution fun
            
            heq[j][i][k]=(w[k]*C[j][i])*(1.0+ 3.0*(cx[k]*u[j][i]+cy[k]*v[j][i]));           
            
           
            // collision equation
            h[j][i][k]=(omegaC1*h[j][i][k]+omegaC*heq[j][i][k]);
            
           }
     }
 }    
 
 // streaming //////////////////////////////////////
 
  mass_streaming();
 
  // tempearture boundary conditionsS //////////////
  
 mass_boundary_conditions();
 
        // POST STREAMING CALCULATION OF TEMP
    
    for (j=0; j<=ny; j++)
    {
        for (i=0; i<=nx; i++)
        {
            C[j][i]=0.0;
            for (k=0; k<=8; k++)
            {
                C[j][i]+=(h[j][i][k]);
                
            }
            
        }
    }     
 //-----------------------------------------------------------------------------
      
 }  // end of main iterattion loop
 //--------------------------------------------------------------

    // nusslet number calculation
// nusslet number calculation: west and East walls


 //error - absolute
  
  error_u =0.0;
  error_v = 0.0; 
   error_T =0.0;
   error_C = 0.0; 
  error = 0.0; 
 
 
  
  for (j=0; j<=ny; j++)
       {
      for (i=0; i<=nx; i++)
           { 
  
   error_u = error_u + fabs(fabs(u[j][i])-fabs(old_u[j][i]));
  // error_v = error_v + fabs(fabs(v[j][i])-fabs(old_v[j][i])); 
   error_T = error_T + fabs(fabs(T[j][i])-fabs(old_T[j][i]));
   error_C = error_C + fabs(fabs(C[j][i])-fabs(old_C[j][i])); 
 
   
           }
      }           

  error=max(max(error_u,error_T),error_C)/(nx*ny);
  
//   cout<<"t="<<t<<"\t Error="<<error<<endl; - i changed - it will only print at 1k intervals to avoid too much output
  if (t % 1000 == 0) {
    cout << "t=" << t << "\t Error=" << error << endl;
}
    
  if (error >= max_error &&  t<=maxit)            goto iterate;  
  /*
   // stream function calculation//--------------------------------------------------------
    streamfun[0][0]=0.0;
     
      for (i=1; i<=nx; i++)
       { 
      rhoavg=0.5*(rho[0][i-1]+rho[0][i]); 
      if (i!=0)   {         
     streamfun[0][i]=streamfun[0][i-1]-(rhoavg*0.5*(v[0][i-1]+v[0][i])); 
		  }
    
     for (j=0; j<=ny; j++)
     {
      
     rhom=0.5*(rho[j][i]+rho[j][i+1]);
     streamfun[j][i]=streamfun[j-1][i]+(rhom*0.5*(u[j-1][i]+u[j][i]));
	      
      } 
      
       } 
*/
// stream function calculation//--------------------------------------------------------
    streamfun[0][0]=0.0;
     
      for (i=0; i<=nx; i++)
       { 
     
      if (i!=0)   {
     rhoavg=0.5*(rho[0][i-1]+rho[0][i]);      
     streamfun[0][i]=streamfun[0][i-1]-(rhoavg*0.5*(v[0][i-1]+v[0][i])); 
		  }
               
	
	    for (j=0; j<=ny; j++)
	    {
            //  for (i=0; i<=nx; i++)
            //   { 
	    
     rhom=0.5*(rho[j][i]+rho[j][i+1]);
     streamfun[j][i]=streamfun[j-1][i]+(rhom*0.5*(u[j-1][i]+u[j][i]));
	      
          } 
      
       } 
//-------------------------------------------------------------------------------------------
// vorticity           

  /* for (i=1; i<nx; i++)
   {
        for (j=1; j<ny; j++)
     {
     
   vor[j][i]=(1.0/8.0)*((u[j-1][i-1]+2.0*u[j-1][i]+u[j-1][i+1])+(v[j-1][i+1]+2.0*v[j][i+1]-v[j+1][i+1])-
              (u[j+1][i+1]+2.0*u[j+1][i]+u[j+1][i-1])-(v[j+1][i-1]+2.0*v[j][i-1]+v[j-1][i-1]));
    }
 }
*/
/*
for (i=0; i<nx; i++)
   {
     for (j=0; j<ny; j++)
    {
     
   vor[j][i]=v[j][i+1]-v[j][i]-u[j+1][i]+u[j][i];
              
    }
 }
 */
// Estimation of local and average Nusselt number//----------------
Nul=0.0;
//Nur=0.0;
for (j=0; j<=ny; j++)
{

Nul+=(T[j][0]-T[j][1])*double (nx);  // Left wall: Heated
//Nur+=(T[j][nx-1]-T[j][nx])*double (nx);
avgNul=Nul/double(ny);
//avgNur=Nur/double(ny);
}
//////////////////////////
Nub=0.0;
//Nur=0.0;
for (i=0; i<=nx; i++)
{

Nub+=(T[0][i]-T[1][i])*double (nx);  // Bottom wall: Heated
//Nur+=(T[j][nx-1]-T[j][nx])*double (nx);
avgNub=Nub/double(ny);
//avgNur=Nur/double(ny);
}

// Estimation of local and average Sherwood number//----------------
Shl=0.0;
//Nur=0.0;
// Estimation of local and average Sherwood number//----------------
for (j=0; j<=ny; j++)
{

Shl+=(C[j][0]-C[j][1])*double (nx);
//Shr+=(C[j][nx-1]-C[j][nx])*double (nx);
avgShl=Shl/double(ny);
//avgShr=Shr/double(ny);
}
//////////
Shb=0.0;
//Nur=0.0;
// Estimation of local and average Sherwood number//----------------
for (i=0; i<=nx; i++)
{

Shb+=(C[0][i]-C[1][i])*double (nx);
//Shr+=(C[j][nx-1]-C[j][nx])*double (nx);
avgShb=Shb/double(ny);
//avgShr=Shr/double(ny);
///////////////////////////////////////////////////
}
/*
// nusslet number calculation: west and East walls

//Nul1[j]=0.0;
//Nur1[j]=0.0;
for (j=0; j<=ny; j++)
{

Nul1[j]=(T[j][0]-T[j][1]);
Nur1[j]=(T[j][nx-1]-T[j][nx]);

}
*/

//-----------------------------------------------------------------------------------------

 
//---------------------------------------------------------------------------------------------
  cout<<"\n "<<endl;
  
 /*
// ofstream outfile("AR1.0localNusseltnumberRa10e5.xls");
 ofstream outfile("Y0.5UVTNuAR1.0halfhotandcoldRa1e4.xls");
 //ofstream outfile("Y0.5TempNaturalconvPr10AR1Ra10e4.xls");
  
  
  
  
  
  //   for (j=0; j<=ny; j++)
   // {
       for (i=0; i<=nx; i++)
      {
      //   p=double (j)/double (ny);
         q=double (i)/double (ny);
 // cout<<"Y="<<p<<"\t X="<<q<<"\t j="<<j<<"\t i="<<i<<"\t ux="<<u[j][i]/uo<<"\t uy="<<v[j][i]/uo<<"\t T="<<T[j][i]<<endl;
   //  
   // outfile<<" "<<p<<"\t "<<q<<"\t "<<T[j][i]<<"\t "<<streamfun[j][i]<<"\t "<<vor[j][i]<<"\t "<<u[j][i]/V<<"\t "<<v[j][i]/V<<endl;
 // outfile<<" "<<j<<"\t "<<p<<"\t "<<i<<"\t "<<q<<"\t "<<u[j][i]/V<<"\t "<<v[j][i]/V<<endl;
   // outfile<<" "<<j<<"\t "<<p<<"\t "<<Nul1[j]<<"\t "<<Nur1[j]<<endl;
  //   outfile<<"Title=Natural Convection in cavity"<<"Variables= Y, X, Temperature"<<endl<<"\n";
    //  cout<<" "<<q<<"\t "<<p<<"\t "<<vor[j][i]<<endl;
       // outfile<<" "<<p<<"\t "<<q<<"\t "<<T[j][i]<<endl;
      //  outfile<<" "<<j<<"\t "<<p<<"\t "<<T[j][nx/2]<<"\t "<<u[j][nx/2]/V<<"\t "<<v[j][nx/2]/V<<endl;
         outfile<<" "<<i<<"\t "<<q<<"\t "<<T[ny/2][i]<<"\t "<<u[ny/2][i]/V<<"\t "<<v[ny/2][i]/V<<endl;
   // }
  }
  

  outfile.close();
*/
  
  
   /*
   ofstream outfile("contoursRa10e3.xls");
  // ofstream outfile("bottomheatedcontoursRa10e5.xls");
  
  for (j=0; j<=ny; j++)
  {
 for (i=0; i<=nx; i++)
   {
        p=double (j)/double (ny);
        q=double (i)/double (ny);
 // cout<<"Y="<<p<<"\t X="<<q<<"\t j="<<j<<"\t i="<<i<<"\t ux="<<u[j][i]/uo<<"\t uy="<<v[j][i]/uo<<"\t T="<<T[j][i]<<endl;
   //   outfile<<" "<<p<<"\t X="<<q<<"\t "<<streamfun[j][i]<<endl;
   outfile<<" "<<p<<"\t "<<q<<"\t "<<T[j][i]<<"\t "<<streamfun[j][i]<<"\t "<<vor[j][i]<<endl;
 // outfile<<" "<<j<<"\t "<<p<<"\t "<<i<<"\t "<<q<<"\t "<<u[j][i]/V<<"\t "<<v[j][i]/V<<endl;
   
  //   outfile<<"Title=Natural Convection in cavity"<<"Variables= Y, X, Temperature"<<endl<<"\n";
      // outfile<<" "<<j<<"\t "<<p<<"\t "<<Nul1[j]<<endl;
       // outfile<<" "<<p<<"\t "<<q<<"\t "<<T[j][i]<<endl;
       // outfile<<" "<<j<<"\t "<<p<<"\t "<<T[j][nx/2]<<"\t "<<u[j][nx/2]/V<<"\t "<<v[j][nx/2]/V<<endl;
      //   outfile<<" "<<i<<"\t "<<q<<"\t "<<T[ny/2][i]<<"\t "<<u[ny/2][i]/V<<"\t "<<v[ny/2][i]/V<<endl;
    }
  }
  

  outfile.close();
  */

 // cout<<"Richardson Number="<<Ri<<endl;
//   cout<<"Raylegh number="<<Ra<<"\t AR= "<<double(nx)/double(ny)<<"\t Prandtl number="<<Pr<<endl<<"\t Nx="<<nx<<"\n"; 
  cout<<"Raylegh number="<<Ra<<"\t Ha="<<Ha<<"\t AR= "<<double(nx)/double(ny)<<"\t Prandtl number="<<Pr<<endl<<"\t Nx="<<nx<<"\n"; // - i changed - added Hartmann number
 cout<<"Velocity of wall="<<ulid<<endl<<"\n";
 cout<<"avgNu left="<<avgNul<<endl<<"\n";
 cout<<"avg Sh left="<<avgShl<<endl<<"\n";
 cout<<"avgNu bottom="<<avgNub<<endl<<"\n";
 cout<<"avg Sh bottom="<<avgShb<<endl<<"\n";
 // cout<<"avgNu north="<<avgNun<<"\t avgNu south="<<avgNus<<"\t Local Nu at north="<<Nun<<"\t Local Nu at east="<<Nus<<endl<<"\n";
  //cout<<"Avg Nu"<<avgNu<<endl;  
 // cout<<" Kinematic Viscocity="<<kinevisc<<"\t thermal Diffusivity="<<Tdiffusivity<<"\t Charct. V="<<Vref<<endl<<"\n";
 cout<<"**************Simulation is Complete***************"<<endl<<"\n";

      
 }

void flow_streaming()
{
// STREAMING PROCESS-------------------------
 
for (j=0; j<=ny; j++)
{
    for (i=nx; i>=1; i--)
    {
        // f1 streaming
        f[j][i][1]=f[j][i-1][1];
        
    }
    
    for (i=0; i<=nx-1; i++)
    {
        //f3 streaming
        f[j][i][3]=f[j][i+1][3];
        
    }
 }
 
 for (j=ny; j>=1; j--)
 {
     for (i=0; i<=nx; i++)
     {
         // f2 streaming
         f[j][i][2]=f[j-1][i][2];
         
     }
     
     for (i=nx; i>=1; i--)
     {
         //f5 streaming
         f[j][i][5]=f[j-1][i-1][5];
         
     }
     
     for (i=0; i<=nx-1; i++)
     {
         f[j][i][6]=f[j-1][i+1][6];
         
     }
 }
 
 for (j=0; j<=ny-1; j++)
 {
     for (i=0; i<=nx; i++)
     {
         f[j][i][4]=f[j+1][i][4];
         
     }
     
     for (i=0; i<=nx-1; i++)
     {
         f[j][i][7]=f[j+1][i+1][7];
         
     }
     
     for (i=nx; i>=1; i--)
     {
         f[j][i][8]=f[j+1][i-1][8];
         
     }
 }
 
}

void heat_streaming()
{
 // streaming of g
  for (j=0; j<=ny; j++)
{
    for (i=nx; i>=1; i--)
    {
        // f1 streaming
        g[j][i][1]=g[j][i-1][1];
        
    }
    
    for (i=0; i<=nx-1; i++)
    {
        //f3 streaming
        g[j][i][3]=g[j][i+1][3];
        
    }
 }
 
 for (j=ny; j>=1; j--)
 {
     for (i=0; i<=nx; i++)
     {
         // f2 streaming
         g[j][i][2]=g[j-1][i][2];
         
     }
     
     for (i=nx; i>=1; i--)
     {
         //f5 streaming
         g[j][i][5]=g[j-1][i-1][5];
         
     }
     
     for (i=0; i<=nx-1; i++)
     {
         g[j][i][6]=g[j-1][i+1][6];
         
     }
 }
 
 for (j=0; j<=ny-1; j++)
 {
     for (i=0; i<=nx; i++)
     {
       g[j][i][4]=g[j+1][i][4];
         
     }
     
     for (i=0; i<=nx-1; i++)
     {
         g[j][i][7]=g[j+1][i+1][7];
         
     }
     
     for (i=nx; i>=1; i--)
     {
         g[j][i][8]=g[j+1][i-1][8];
         
     }
 }

 
}

void heat_boundary_conditions()
{

 // west Boundary, T=1; in the middle of domain
 
  for (j=0; j<=ny; j++)
    {
        
        // Left wall: bounceback bounadry conditions:Tw=1.0 : high temp
           
              g[j][0][1]=Tw*(w[1]+w[3])-g[j][0][3];              
                g[j][0][5]=Tw*(w[5]+w[7])-g[j][0][7];              
                  g[j][0][8]=Tw*(w[8]+w[6])-g[j][0][6];     
    
	
          // east wall: cold wall
	  
	         g[j][nx][3]=-g[j][nx][1];              
                  g[j][nx][7]=-g[j][nx][5];              
                   g[j][nx][6]=-g[j][nx][8];     
      }
      
      
     // adiabatic part of top AND bottom walls: adiabgic 
     
   /* 
     for (j=h1; j<=ny; j++)
     {
       g[j][0][k]=g[j][1][k];
     }
    */
  // heated bottom wall
   /* 
    for (i=0; i<=nx; i++)
    {
     g[0][i][2]=Tw*(w[2]+w[4])-g[0][i][4];
      g[0][i][5]=Tw*(w[5]+w[7])-g[0][i][7];
       g[0][i][6]=Tw*(w[6]+w[8])-g[0][i][8];
     }*/
    // adiabatic: north and bottom wall BC  
    for (i=0; i<=nx; i++)
    {
        for (k=0; k<=8; k++)
        {
      
        g[ny][i][k]=g[ny-1][i][k];    
          g[0][i][k]=g[1][i][k]; 
        }       
          
    }   
}
//////////////////////////////////////////////////////////////
	
	void mass_streaming()
{
 // streaming of g
  for (j=0; j<=ny; j++)
{
    for (i=nx; i>=1; i--)
    {
        // f1 streaming
        h[j][i][1]=h[j][i-1][1];
        
    }
    
    for (i=0; i<=nx-1; i++)
    {
        //f3 streaming
        h[j][i][3]=h[j][i+1][3];
        
    }
 }
 
 for (j=ny; j>=1; j--)
 {
     for (i=0; i<=nx; i++)
     {
         // f2 streaming
         h[j][i][2]=h[j-1][i][2];
         
     }
     
     for (i=nx; i>=1; i--)
     {
         //f5 streaming
         h[j][i][5]=h[j-1][i-1][5];
         
     }
     
     for (i=0; i<=nx-1; i++)
     {
         h[j][i][6]=h[j-1][i+1][6];
         
     }
 }
 
 for (j=0; j<=ny-1; j++)
 {
     for (i=0; i<=nx; i++)
     {
       h[j][i][4]=h[j+1][i][4];
         
     }
     
     for (i=0; i<=nx-1; i++)
     {
         h[j][i][7]=h[j+1][i+1][7];
         
     }
     
     for (i=nx; i>=1; i--)
     {
         h[j][i][8]=h[j+1][i-1][8];
         
     }
 }

 
}

void mass_boundary_conditions()
{

 
  for (j=0; j<=ny; j++)
    {
      /*
      // west wall: Unifromly Concentrated
           
              h[j][0][1]=Ch*(w[1]+w[3])-h[j][0][3];              
                h[j][0][5]=Ch*(w[5]+w[7])-h[j][0][7];              
                  h[j][0][8]=Ch*(w[8]+w[6])-h[j][0][6];     
    
      */  
        // Left wall: Low Concentrated
           
              h[j][0][1]=-h[j][0][3];              
                h[j][0][5]=-h[j][0][7];              
                  h[j][0][8]=-h[j][0][6];     
    
	
          // Right wall: Higher concentration wall
	  
	           h[j][nx][3]=Ch*(w[1]+w[3])-h[j][nx][1];              
                  h[j][nx][7]=Ch*(w[7]+w[5])-h[j][nx][5];              
                   h[j][nx][6]=Ch*(w[6]+w[8])-h[j][nx][8];     
      }
      
      for (i=0; i<=nx; i++)
       {
         }
     // No mass transfer from north wall
    for (i=0; i<=nx; i++)
    {
        for (k=0; k<=8; k++)
        {
              
        h[ny][i][k]=h[ny-1][i][k];
              
        h[0][i][k]=h[1][i][k];    
        }       
          
    }       
   
}