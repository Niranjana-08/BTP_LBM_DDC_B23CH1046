/*
 * Double-Diffusive MHD Natural Convection in an Enclosed Cavity
 * Focused Parametric Study: Varying Rayleigh (Ra) and Hartmann (Ha) Numbers Only
 * Fixed Constants: Pr = 0.054, Le = 1.0, N = 0.0
 * Saves all output CSV files into an 'output_csv_data' subfolder.
 */

#include<iostream>
#include<stdio.h>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<sstream>

#ifdef _WIN32
    #include <direct.h>
#else
    #include <sys/stat.h>
#endif

# define nx 200
# define ny 200
# define Ma 0.1   // Mach number
# define AP 0     // Velocity parameter (A=U*L/alpha)
# define ci 1.0/sqrt(3)  // lattice speed constant

void flow_streaming();
void heat_streaming();
void heat_boundary_conditions();
void mass_streaming();
void mass_boundary_conditions();

using namespace std;

double f[800][800][9], feq[800][800][9], rho[800][800], u[800][800], v[800][800], w[9], cx[9], cy[9];
double streamfun[800][800], rhoavg, rhom, old_T[800][800], error_T, error_C, old_u[800][800], old_v[800][800], old_C[800][800];
double g[800][800][9], geq[800][800][9], T[800][800], Tw=1.0, Tc=0.0, Tref, Ch=1.0, Cl=0.0, Cref, force[800][800][9];
double h[800][800][9], heq[800][800][9], C[800][800];

int i = nx + 2, j = ny + 2, k = 9;
double vor[800][800];

// Dynamic Parameter Holders
double Ra, Ha, Le, N, Pr;

// Derived Physical Constants (recalculated per run)
double kinevisc, omega, Tdiffusivity, omegaT, Diffusivity, omegaC;
double omega1, omegaT1, omegaC1, gbetaT, betaT, betaC, A, ulid;

void recalculate_physics() {
    kinevisc = sqrt((Ma * Ma * ny * ny * Pr * ci * ci) / Ra);
    omega = 1.0 / (3.0 * kinevisc + 0.5);
    Tdiffusivity = kinevisc / Pr;
    omegaT = 1.0 / (3.0 * Tdiffusivity + 0.5);
    Diffusivity = kinevisc / (Pr * Le);
    omegaC = 1.0 / (3.0 * Diffusivity + 0.5);
    omega1 = (1.0 - omega);
    omegaT1 = (1.0 - omegaT);
    omegaC1 = (1.0 - omegaC);
    gbetaT = (Ra * kinevisc * kinevisc) / ((ny * ny * ny) * Pr);
    betaT = gbetaT / 10.0;
    betaC = N * betaT;
    A = (Ha * Ha * kinevisc) / (ny * ny);
    ulid = (AP * Tdiffusivity) / ny;
}

void print_progress_bar(int current, int total, double error, int iteration) {
    int barWidth = 30;
    float progress = (float)current / total;
    cout << "\r[Run " << setw(2) << current << "/" << total << "] [";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) cout << "=";
        else if (i == pos) cout << ">";
        else cout << " ";
    }
    cout << "] " << int(progress * 100.0) << "% | Iter: " << iteration << " | Err: " << scientific << setprecision(2) << error << flush;
}

int main ()
{
    // Create output folder automatically if it doesn't exist
    #ifdef _WIN32
        _mkdir("output_csv_data");
    #else
        mkdir("output_csv_data", 0777);
    #endif

    // Targeted Parameter Arrays (3 x 3 = 9 Total Runs)
    double Ra_vals[] = {1000.0, 10000.0, 100000.0};
    double Ha_vals[] = {0.0, 50.0, 100.0};

    // Fixed Baseline Constants
    double N_fixed  = 0.0;
    double Pr_fixed = 0.054;
    double Le_fixed = 1.0;

    int total_runs = 3 * 3;
    int run_counter = 0;

    // Save master log inside 'output_csv_data' folder
    ofstream master_csv("output_csv_data/summary_results_Ra_Ha_study.csv");
    master_csv << "Run,Ra,Ha,N,Pr,Le,avgNu_left,avgSh_left,avgNu_bottom,avgSh_bottom,Converged_Iter\n";

    cout << "========================================================================\n";
    cout << "   STARTING TARGETED 9-RUN Ra vs Ha PARAMETRIC LBM STUDY               \n";
    cout << "========================================================================\n\n";

    for (int r = 0; r < 3; r++) {
        for (int h_idx = 0; h_idx < 3; h_idx++) {

            run_counter++;
            Ra = Ra_vals[r];
            Ha = Ha_vals[h_idx];
            N  = N_fixed;
            Pr = Pr_fixed;
            Le = Le_fixed;

            recalculate_physics();

            int t, maxit = 80000000;
            double error, error_u, error_T, error_C, max_error = 1.e-8;
            double Nul, Nub, Shl, Shb, avgNul, avgNub, avgShl, avgShb;

            // LBM Weight and Velocity Vectors
            w[0]=(4.0/9.0); w[1]=(1.0/9.0); w[2]=(1.0/9.0); w[3]=(1.0/9.0); w[4]=(1.0/9.0);
            w[5]=(1.0/36.0); w[6]=(1.0/36.0); w[7]=(1.0/36.0); w[8]=(1.0/36.0);

            cx[0]=0.0; cx[1]=1.0; cx[2]=0.0; cx[3]=-1.0; cx[4]=0.0;
            cx[5]=1.0; cx[6]=-1.0; cx[7]=-1.0; cx[8]=1.0;

            cy[0]=0.0; cy[1]=0.0; cy[2]=1.0; cy[3]=0.0; cy[4]=-1.0;
            cy[5]=1.0; cy[6]=1.0; cy[7]=-1.0; cy[8]=-1.0;

            Tref = (Tw + Tc) / 2.0;
            Cref = (Ch + Cl) / 2.0;

            // Grid initialization
            for (j=0; j<=ny; j++) {
                for (i=0; i<=nx; i++) {
                    rho[j][i] = 2.7;
                    u[j][i] = 0.0;
                    v[j][i] = 0.0;
                    T[j][i] = 0.0;
                    C[j][i] = 0.0;
                    for (k=0; k<=8; k++) {
                        feq[j][i][k] = (rho[j][i]*w[k]);
                        f[j][i][k]   = feq[j][i][k];
                        g[j][i][k]   = 0.0;
                        h[j][i][k]   = 0.0;
                    }
                }
            }

            t = 0;
            iterate:
            t = t + 1;
            {
                for (j=0; j<=ny; j++) {
                    for (i=0; i<=nx; i++) {
                        old_u[j][i] = u[j][i];
                        old_v[j][i] = v[j][i];
                        old_T[j][i] = T[j][i];
                        old_C[j][i] = C[j][i];
                    }
                }

                // Hydrodynamic Collision
                for (j=0; j<=ny; j++) {
                    for (i=0; i<=nx; i++) {
                        for (k=0; k<=8; k++) {
                            feq[j][i][k]=(rho[j][i]*w[k])*(1.0 + 3.0*(u[j][i]*cx[k]+v[j][i]*cy[k])
                                           + (4.50*(u[j][i]*cx[k]+v[j][i]*cy[k])*(u[j][i]*cx[k]+v[j][i]*cy[k]))
                                           - 1.50*(u[j][i]*u[j][i]+v[j][i]*v[j][i]));

                            force[j][i][k]=(3.0*w[k]*cy[k]*rho[j][i])*(gbetaT*((T[j][i]-Tref)+(C[j][i]-Cref)*N) - A*v[j][i]);

                            if (i==0 || i==nx || j==0 || j==ny) force[j][i][k] = 0.0;

                            f[j][i][k]=(omega1*f[j][i][k])+(omega*feq[j][i][k])+force[j][i][k];
                        }
                    }
                }

                flow_streaming();

                // Flow boundary conditions
                for (j=0; j<=ny; j++) {
                    f[j][0][1]=f[j][0][3]; f[j][0][5]=f[j][0][7]; f[j][0][8]=f[j][0][6];
                    f[j][nx][3]=f[j][nx][1]; f[j][nx][7]=f[j][nx][5]; f[j][nx][6]=f[j][nx][8];
                }
                for (i=0; i<=nx; i++) {
                    f[ny][i][4]=f[ny][i][2]; f[ny][i][7]=f[ny][i][5]; f[ny][i][8]=f[ny][i][6];
                    f[0][i][2]=f[0][i][4]; f[0][i][6]=f[0][i][8]; f[0][i][5]=f[0][i][7];
                }

                // Density and velocity update
                for (j=0; j<=ny; j++) {
                    for (i=0; i<=nx; i++) {
                        rho[j][i]=0.0; u[j][i]=0.0; v[j][i]=0.0;
                        for (k=0; k<=8; k++) {
                            rho[j][i]+=f[j][i][k];
                            u[j][i]+=(cx[k]*f[j][i][k]);
                            v[j][i]+=(cy[k]*f[j][i][k]);
                        }
                        u[j][i]/=rho[j][i];
                        v[j][i]/=rho[j][i];
                    }
                }

                // Thermal collision & streaming
                for (j=0; j<=ny; j++) {
                    for (i=0; i<=nx; i++) {
                        T[j][i]=0.0;
                        for (k=0; k<=8; k++) T[j][i]+=(g[j][i][k]);
                    }
                }
                for (j=0; j<=ny; j++) {
                    for (i=0; i<=nx; i++) {
                        for (k=0; k<=8; k++) {
                            geq[j][i][k]=(w[k]*T[j][i])*(1.0+ 3.0*(cx[k]*u[j][i]+cy[k]*v[j][i]));
                            g[j][i][k]=(omegaT1*g[j][i][k]+omegaT*geq[j][i][k]);
                        }
                    }
                }
                heat_streaming();
                heat_boundary_conditions();
                for (j=0; j<=ny; j++) {
                    for (i=0; i<=nx; i++) {
                        T[j][i]=0.0;
                        for (k=0; k<=8; k++) T[j][i]+=(g[j][i][k]);
                    }
                }

                // Concentration collision & streaming
                for (j=0; j<=ny; j++) {
                    for (i=0; i<=nx; i++) {
                        C[j][i]=0.0;
                        for (k=0; k<=8; k++) C[j][i]+=(h[j][i][k]);
                    }
                }
                for (j=0; j<=ny; j++) {
                    for (i=0; i<=nx; i++) {
                        for (k=0; k<=8; k++) {
                            heq[j][i][k]=(w[k]*C[j][i])*(1.0+ 3.0*(cx[k]*u[j][i]+cy[k]*v[j][i]));
                            h[j][i][k]=(omegaC1*h[j][i][k]+omegaC*heq[j][i][k]);
                        }
                    }
                }
                mass_streaming();
                mass_boundary_conditions();
                for (j=0; j<=ny; j++) {
                    for (i=0; i<=nx; i++) {
                        C[j][i]=0.0;
                        for (k=0; k<=8; k++) C[j][i]+=(h[j][i][k]);
                    }
                }
            }

            // Error evaluation
            error_u = 0.0; error_T = 0.0; error_C = 0.0;
            for (j=0; j<=ny; j++) {
                for (i=0; i<=nx; i++) {
                    error_u += fabs(fabs(u[j][i])-fabs(old_u[j][i]));
                    error_T += fabs(fabs(T[j][i])-fabs(old_T[j][i]));
                    error_C += fabs(fabs(C[j][i])-fabs(old_C[j][i]));
                }
            }
            error = max(max(error_u, error_T), error_C) / (nx * ny);

            if (t % 2000 == 0) {
                print_progress_bar(run_counter, total_runs, error, t);
            }

            if (error >= max_error && t <= maxit) goto iterate;

            // Final Metrics Calculation
            Nul = 0.0; Nub = 0.0; Shl = 0.0; Shb = 0.0;
            for (j=0; j<=ny; j++) { Nul += (T[j][0] - T[j][1]) * double(nx); }
            avgNul = Nul / double(ny);

            for (i=0; i<=nx; i++) { Nub += (T[0][i] - T[1][i]) * double(nx); }
            avgNub = Nub / double(ny);

            for (j=0; j<=ny; j++) { Shl += (C[j][0] - C[j][1]) * double(nx); }
            avgShl = Shl / double(ny);

            for (i=0; i<=nx; i++) { Shb += (C[0][i] - C[1][i]) * double(nx); }
            avgShb = Shb / double(ny);

            cout << "\n -> Completed Run " << run_counter << "/" << total_runs 
                 << " [Ra=" << Ra << ", Ha=" << Ha << ", N=" << N << ", Pr=" << Pr << "]"
                 << " | avgNu_L=" << fixed << setprecision(4) << avgNul 
                 << " | avgSh_L=" << avgShl << "\n";

            master_csv << run_counter << "," << Ra << "," << Ha << "," << N << "," << Pr << "," << Le << ","
                       << avgNul << "," << avgShl << "," << avgNub << "," << avgShb << "," << t << "\n";

            // Save grid file inside 'output_csv_data' folder
            stringstream grid_filename;
            grid_filename << "output_csv_data/grid_Ra" << int(Ra) << "_Ha" << int(Ha) << "_N" << int(N) << "_Pr" << Pr << ".csv";
            ofstream grid_file(grid_filename.str().c_str());
            grid_file << "i,j,T,C,u,v\n";
            for (j=0; j<=ny; j++) {
                for (i=0; i<=nx; i++) {
                    grid_file << i << "," << j << "," << T[j][i] << "," << C[j][i] << "," << u[j][i] << "," << v[j][i] << "\n";
                }
            }
            grid_file.close();

        }
    }

    master_csv.close();
    cout << "\n========================================================================\n";
    cout << "   9-RUN TARGETED STUDY COMPLETE! All data saved in 'output_csv_data/'   \n";
    cout << "========================================================================\n";

    return 0;
}

void flow_streaming() {
    for (j=0; j<=ny; j++) {
        for (i=nx; i>=1; i--) f[j][i][1]=f[j][i-1][1];
        for (i=0; i<=nx-1; i++) f[j][i][3]=f[j][i+1][3];
    }
    for (j=ny; j>=1; j--) {
        for (i=0; i<=nx; i++) f[j][i][2]=f[j-1][i][2];
        for (i=nx; i>=1; i--) f[j][i][5]=f[j-1][i-1][5];
        for (i=0; i<=nx-1; i++) f[j][i][6]=f[j-1][i+1][6];
    }
    for (j=0; j<=ny-1; j++) {
        for (i=0; i<=nx; i++) f[j][i][4]=f[j+1][i][4];
        for (i=0; i<=nx-1; i++) f[j][i][7]=f[j+1][i+1][7];
        for (i=nx; i>=1; i--) f[j][i][8]=f[j+1][i-1][8];
    }
}

void heat_streaming() {
    for (j=0; j<=ny; j++) {
        for (i=nx; i>=1; i--) g[j][i][1]=g[j][i-1][1];
        for (i=0; i<=nx-1; i++) g[j][i][3]=g[j][i+1][3];
    }
    for (j=ny; j>=1; j--) {
        for (i=0; i<=nx; i++) g[j][i][2]=g[j-1][i][2];
        for (i=nx; i>=1; i--) g[j][i][5]=g[j-1][i-1][5];
        for (i=0; i<=nx-1; i++) g[j][i][6]=g[j-1][i+1][6];
    }
    for (j=0; j<=ny-1; j++) {
        for (i=0; i<=nx; i++) g[j][i][4]=g[j+1][i][4];
        for (i=0; i<=nx-1; i++) g[j][i][7]=g[j+1][i+1][7];
        for (i=nx; i>=1; i--) g[j][i][8]=g[j+1][i-1][8];
    }
}

void heat_boundary_conditions() {
    for (j=0; j<=ny; j++) {
        g[j][0][1]=Tw*(w[1]+w[3])-g[j][0][3];
        g[j][0][5]=Tw*(w[5]+w[7])-g[j][0][7];
        g[j][0][8]=Tw*(w[8]+w[6])-g[j][0][6];
        g[j][nx][3]=-g[j][nx][1];
        g[j][nx][7]=-g[j][nx][5];
        g[j][nx][6]=-g[j][nx][8];
    }
    for (i=0; i<=nx; i++) {
        for (k=0; k<=8; k++) {
            g[ny][i][k]=g[ny-1][i][k];
            g[0][i][k]=g[1][i][k];
        }
    }
}

void mass_streaming() {
    for (j=0; j<=ny; j++) {
        for (i=nx; i>=1; i--) h[j][i][1]=h[j][i-1][1];
        for (i=0; i<=nx-1; i++) h[j][i][3]=h[j][i+1][3];
    }
    for (j=ny; j>=1; j--) {
        for (i=0; i<=nx; i++) h[j][i][2]=h[j-1][i][2];
        for (i=nx; i>=1; i--) h[j][i][5]=h[j-1][i-1][5];
        for (i=0; i<=nx-1; i++) h[j][i][6]=h[j-1][i+1][6];
    }
    for (j=0; j<=ny-1; j++) {
        for (i=0; i<=nx; i++) h[j][i][4]=h[j+1][i][4];
        for (i=0; i<=nx-1; i++) h[j][i][7]=h[j+1][i+1][7];
        for (i=nx; i>=1; i--) h[j][i][8]=h[j+1][i-1][8];
    }
}

void mass_boundary_conditions() {
    for (j=0; j<=ny; j++) {
        h[j][0][1]=-h[j][0][3];
        h[j][0][5]=-h[j][0][7];
        h[j][0][8]=-h[j][0][6];
        h[j][nx][3]=Ch*(w[1]+w[3])-h[j][nx][1];
        h[j][nx][7]=Ch*(w[7]+w[5])-h[j][nx][5];
        h[j][nx][6]=Ch*(w[6]+w[8])-h[j][nx][8];
    }
    for (i=0; i<=nx; i++) {
    }
    for (i=0; i<=nx; i++) {
        for (k=0; k<=8; k++) {
            h[ny][i][k]=h[ny-1][i][k];
            h[0][i][k]=h[1][i][k];
        }
    }
}