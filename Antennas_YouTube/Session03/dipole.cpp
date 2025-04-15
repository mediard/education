// This code was written with the help of ChatGPT as part of the Antenna Course at @ComputationalElectromagnetics YouTube channel.
// The mathematical equations were broken into smaller pieces and function, clear algorithms were provided.
// Types of inputs and outputs were defined. ChatGPT provided initial C++ codes.
// Codes were modified and debugged.
// by: Mehdi Ardavan
#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <algorithm>
#include <tuple>
#include <fstream>




// Constants
const double c_0 = 299792458.0;  // Speed of light in m/s
const double mu_0 = 4.0 * M_PI * 1e-7;  // Permeability of free space in H/m
const double epsilon_0 = 8.854187817e-12;  // Permittivity of free space in F/m
const double TOL = 1e-10; // A small threshold for checking whether a number is close to zero

std::tuple<double, double, double> computeValues(double l, double a, double phi_M, int i, int j, int N){//, double l ,double N) {
    // Calculating x_m, y_m, and z_m
    double x_m = -a * std::cos(phi_M);
    double y_m = a - a * std::sin(phi_M);
    double z_m = (i - j) * ((l) / (N)); // using M_PI as value for pi

    return std::make_tuple(x_m, y_m, z_m);
}

std::tuple<double, double, double> cartesianToSpherical(std::tuple<double,double,double> xyz_in) {
    double r, theta, phi;
    double x = std::get<0>(xyz_in);  
    double y = std::get<1>(xyz_in);
    double z = std::get<2>(xyz_in);

    // Handling the special conditions first
    if (x == 0 && y == 0 && z == 0) {
        return std::make_tuple(0, 0, 0);
    }

    if (x == 0 && y == 0 && z != 0) {
        return std::make_tuple(z, 0, 0);
    }

    // Computing r
    r = std::sqrt(x * x + y * y + z * z);
    // Computing theta
    theta = std::acos(z / r);

    if (y == 0 && x != 0) {
        phi = std::acos(x / std::sqrt(x * x + y * y));
    } else {
        // Computing phi with sign of y considered
        int sgn_y = (y > 0) - (y < 0); // This will return 1 for positive y, -1 for negative y, and 0 for y=0
        phi = -M_PI * (sgn_y - 1) + sgn_y * std::acos(x / std::sqrt(x * x + y * y));
    }

    return std::make_tuple(r, theta, phi);
}


std::pair<std::complex<double>, std::complex<double>> computeFields(double r, double theta, double phi, double f,   double a2z) {
    std::complex<double> cj(0, 1);  // Imaginary unit

    // Calculate necessary constants
    double omega = 2 * M_PI * f;
    double lambda = c_0 / f;
    double k = omega * std::sqrt(mu_0 * epsilon_0);
    double eta = std::sqrt(mu_0 / epsilon_0);
    
    // Calculate E_r and E_theta
    std::complex<double> E_r = (eta * std::cos(theta)  / (2. )) * (a2z/(r*r) + a2z / (cj * k * r*r*r)) * std::exp(-cj * k * r);
    std::complex<double> E_theta = (cj * eta * k * std::sin(theta) / (4. )) * (a2z/r + a2z / (cj * k * r*r) - a2z / (k * r * k * r*r)) * std::exp(-cj * k * r);

    return {E_r, E_theta};
}

std::complex<double> computeEz(double r, double theta, double phi, const std::complex<double>& Er, const std::complex<double>& Etheta) {
    return std::cos(theta) * Er -  std::sin(theta) * Etheta;
}

std::complex<double> computeSumEzm(double z, double z_prime,  double a, double a2z) {
    std::complex<double> sumEzm(0, 0);
    std::tuple<double, double, double> xyz;
    std::tuple<double, double, double> rtp;
    std::pair<std::complex<double>, std::complex<double>> Etp;
    std::complex<double> cj(0, 1);  // Imaginary unit
    std::complex<double> Ez,Er,Etheta;
    double r, theta,phi;
    double lambda =1.;
    double f = c_0/lambda;
    //for (int m = 1; m <= M; ++m) {
        //double phi_m = m * 2.0 * M_PI / M;
        //xyz = computeValues(l,a, phi_m,  i,  j,N);
        xyz = std::make_tuple(0., a, z-z_prime);
        rtp = cartesianToSpherical( xyz );
        r = std::get<0>(rtp);  
        theta = std::get<1>(rtp);
        phi = std::get<2>(rtp);
        Etp = computeFields( r,  theta, phi, f, a2z);
        Er = Etp.first;
        Etheta = Etp.second;
        Ez = computeEz(r, theta, phi, Er, Etheta);
        sumEzm = Ez;// computeEz(phi_m, f, I_0, a, l, i, j, M, N);
        //computeSumEzm(double z, double z_prime,  double a, double a2z)
        //std::cout<<"z = "<< z << "   ,  z_prime = " << z_prime << "   , Ez = " <<sumEzm<<std::endl;
    //}



    return sumEzm;
}




std::complex<double> compositeSimpson(double z_a, double z_b, double z_m, int n_slices) {
    /*if (N % 2 != 0) {
        std::cerr << "N must be even!" << std::endl;
        return -1;
    }*/
    if (n_slices % 2 != 0) {
        std::cout << "n_slices is odd. Making it even..." << std::endl;
        n_slices++;
    }
    //int N = 38;
    double a=0.005*1.;
    double z,  z_prime,   a2z;
    double h = (z_b - z_a) / n_slices;
    //computeSumEzm(double z, double z_prime,  double a, double a2z) 
    std::complex<double> integral = computeSumEzm(z_m, z_a,  a, 1./M_PI) + computeSumEzm(z_m, z_b,  a, 1./M_PI);
    //std::complex<double> integral = computeSumEzm(z_a, f, I_0, a, l,N) + computeSumEzm(z_b, f, I_0, a, l,N);

    for (int i = 1; i <= n_slices-1; i++) {
        double z = z_a + i * h;
        if (i % 2 == 0) {
            integral += 2. * computeSumEzm(z_m, z,  a, 1./M_PI);
        } else {
            integral += 4. * computeSumEzm(z_m, z,  a, 1./M_PI);
        }
    }

    integral *= h/3.0;

    return integral;
}



std::complex<double> g(double z, double zprime, double a) {
    const double lambda = 1.0;
    
    double eta0 = std::sqrt(mu_0 / epsilon_0);
    double k = 2.0 * M_PI / lambda;
    double RR = std::sqrt(std::pow(z - zprime, 2) + std::pow(a, 2));

    std::complex<double> cj(0.0, 1.0);  // Complex unit j
    std::complex<double> g_value = 
        -cj * eta0 / k * std::exp(-cj * k * RR) / (4.0 * M_PI * std::pow(RR, 5)) *
        ((1.0 + cj * k * RR) * (2.0 * std::pow(RR, 2) - 3.0 * std::pow(a, 2)) + 
         std::pow(k * a * RR, 2));

    return g_value;
}




std::vector<std::vector<std::complex<double> > > computeG(double a, double l, double f, int k_prime) {
    int N = 2 * k_prime + 1;
    std::vector<std::vector<std::complex<double> > > G(N, std::vector<std::complex<double> >(N, 0));
    double zm,zs, z_a, z_b, dZ2;
    dZ2 = l/N/2.;
    int n_slices = 100 * (0.005/a) * (l/N*80);
    std::cout<<"n_slices = "<< n_slices<<std::endl;
    for (int i = 0; i < N; ++i) {
        zm = -l/2.+i*l/N+dZ2;
        for (int j = 0; j < N; ++j) {
            //if (i != j) {
                zs = -l/2.+j*l/N+dZ2;
                z_a = zs - dZ2;
                z_b = zs + dZ2;
                //std::cout<<"hello"<<std::endl;
                G[i][j] = - compositeSimpson( z_a,  z_b, zm,n_slices);
                //std::cout<<"G = "<<G[i][j]<<",   z_a = "<<z_a<<",   z_b = "<<z_b<<",   zm = "<<zm<<std::endl;
                //G[i][j] = compositeSimpson(z_a, z_b, f, I_0, a, l, N+1);
                //G[i][j] = computeSumEzm(f, I_0, a, l, i, j, M, N);  
            //}
        }
    }
    return G;
}



std::vector<std::complex<double>> computeEtans(double l, int k_prime) {
    int N = 2 * k_prime + 1;
    std::vector<std::complex<double>> Etans(N, 0);
    double delta_z = l / N;
    std::cout<<" !!!!!!!!!!!!!!!! delta_z = "<<delta_z<<std::endl;
    Etans[k_prime] = -1.0 / delta_z;

    return Etans;
}


void forwardElimination(std::vector<std::vector<std::complex<double>>>& G, std::vector<std::complex<double>>& Etans, int N) {
    for (int k = 0; k < N-1; k++) {
        int maxRow = k;
        for (int i = k+1; i < N; i++) {
            if (std::abs(G[i][k]) > std::abs(G[maxRow][k])) {
                maxRow = i;
            }
        }
        std::swap(G[k], G[maxRow]);
        std::swap(Etans[k], Etans[maxRow]);

        for (int i = k+1; i < N; i++) {
            std::complex<double> factor = G[i][k] / G[k][k];
            for (int j = k; j < N; j++) {
                G[i][j] -= factor * G[k][j];
            }
            Etans[i] -= factor * Etans[k];
        }
    }
}

std::vector<std::complex<double>> backSubstitution(const std::vector<std::vector<std::complex<double>>>& G, const std::vector<std::complex<double>>& Etans, int N) {
    std::vector<std::complex<double>> I(N);
    for (int i = N-1; i >= 0; i--) {
        I[i] = Etans[i];
        for (int j = i+1; j < N; j++) {
            I[i] -= G[i][j] * I[j];
        }
        I[i] /= G[i][i];
    }
    return I;
}

std::vector<std::complex<double>> solveUsingElimination(const std::vector<std::vector<std::complex<double>>>& G, const std::vector<std::complex<double>>& Etans) {
    int N = G.size();
    std::vector<std::vector<std::complex<double>>> G_copy = G;
    std::vector<std::complex<double>> Etans_copy = Etans;

    forwardElimination(G_copy, Etans_copy, N);
    return backSubstitution(G_copy, Etans_copy, N);
}


int main(){
    // z_a = -0.235,   z_b = -0.0783333,   zm = -0.156667
    
    std::cout<<computeSumEzm(0, 0,  0.005, 1./M_PI )<<std::endl;
    std::cout<<"hello\n";
    std::cout<<g(0, 0,  0.005)<<std::endl;

    double lambda = 1.;
    double f = c_0/lambda;
    double l = 0.47 * lambda;
    double a = 0.005 * lambda;
    int k_prime = 20;
    int N=2*k_prime+1;
    
    std::vector<std::complex<double>> Etans;
    Etans = computeEtans(l, k_prime);

    std::vector<std::vector<std::complex<double> > > G;
    G=computeG(a, l, f, k_prime);
    //std::cout<<-compositeSimpson(-0.234999999, -0.230346546, -0.232673272)<<std::endl;
     std::cout << "\n G :" << std::endl;
    for (int i = 0; i < N; i++) {
        for (int j=0; j<N; j++) {
        //std::cout << "G( " << i<<","<<j<<  ")= " << G[i][j] << std::endl;
        }
        /*
        std::cout << "Press any key to continue . . .";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // clear any existing input
        std::cin.get();  // wait for user to press any key
        */
    }
    std::cout << "G( " << N-1<<","<<N-1<<  ")= " << G[N-1][N-1] << std::endl;
    std::vector<std::complex<double>> Currents;
    Currents = solveUsingElimination( G, Etans);

    std::ofstream outFile("CPP_MoF.txt");
    for (int i = 0; i < Currents.size(); i++) {
        std::cout << "Magnitude of Current " << i + 1 << ": " << std::abs(Currents[i]) << std::endl;
        double zp = -l/2.0 + l/(2.0*N) + i*l/N;
        outFile << zp << "\t" << std::abs(Currents[i]) << std::endl;
    }
    outFile.close();


    std::vector<std::complex<double>> results(N, 0);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            results[i] += G[i][j] * Currents[j];
        }
    }

    std::cout << "\\nResults of G x Currents:" << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << "Result " << i + 1 << ": " << results[i] << std::endl;
    }



    std::cout<<"Z0 = -1/cuurent(k_prime) =  "<<std::complex<double> (-1.,0.)/Currents[k_prime]<<std::endl;
    //std::cout<<"Z0 = Etan/Current(k_prime)= "<<Etans[k_prime]/Currents[k_prime]<<std::endl;
    //std::cout<<"Current(k_prime+1)="<<Currents[k_prime]<<std::endl;

    return 0;


}
