#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <complex>
using namespace std;

#define PI 3.14159265358979323846
//k_B = 1;
const int q = 3; //number of spin states; we have 3: 0, 1 and 2
const int L = 16; //number of sites; L = N 
const double T = 0.25; //temperature in units J 
//changing temperature has to be done manually, a bit hacky but it will have to do for now  

int dim = 1; //dimension of the lattice; 1D chain case

const int N = L; 
const double pconnect = 1. - exp(-2. / T); //probability of connecting two sites; this expression works because our temperature is in units of J 
//might need to make a function for connecting probability so that it can be computed for different temperatures 

const int Nclusters = 1; //number of cluster builds in one MC step 
const int NEsteps = 100000; //number of equilibrium MC steps - for equilibration/burn-in i assume 
const int NMCsteps = 100000; //number of MC steps for measurement
const int Nbins = 10; //number of measurement bins/experiments 

vector<int> Spins(N); //array of spins 
vector<int> M(q); //number of spins in different states - will be used in flipping the spins 
vector<complex<double>> W(q); //oder parameter weights 
vector<complex<double>> m0r(N); //order parameter for each site
vector<complex<double>> mr(N); //order parameter for each site
vector<complex<double>> Cr(N); //order parameter for each site


//function to initialize the lattice with the spins 
enum directions {right, left}; //1D chain case
int indx(int x) {return x;}; //make an indx on every site 
int xpos(int i) {return i%L;} 

int Nbr(int i, int dir)
{
    int x = xpos(i);

    switch(dir)
    {
        case ::right: return indx((x+1)%L); //return index of the right neighbour 
        case ::left: return indx((x-1+L)%L); //return index of the left neighbour
        default: return -1; //if the direction is not right or left, return -1
    }

}

//implementing cluster-flipping algorithm 
void FlipAndBuild_from(int s) //might need to implement temperature as an argument later as well when looking for phase transitions 
{
    int oldstate(Spins[s]); 
    int newstate((Spins[s] + 1) % q); //save the old state of the site
    //newstate is modulated so that it is a state between 0,1 and 2 in our case of Potts model with q = 3
    Spins[s] = newstate; //flip the spin of the site
    M[oldstate]--; //decrease the number of spins in the old state by 1
    M[newstate]++; //increase the number of spins in the new state by 1

    for (int dir = 0; dir < 2; dir++) //loop over all neighbours s
    //2 neighbours in 1D case, 4 neighbours in 2D case
        {
            int j = Nbr(s, dir); //get the index of the neighbour
            if (Spins[j] == oldstate && std::rand() / (RAND_MAX + 1.) < pconnect) //if the neighbour has the same spin as the site and the random number is less than the probability of connecting the two sites 
            {
                FlipAndBuild_from(j); //call the function again 
            }
        }
}


int main()
{
    for(int s = 0; s < q; s++)
    {
        W[s] = complex<double>(cos(2 * PI * s / q), sin(2 * PI * s / q));
    } 
    
  
    for(int i = 0; i < N; i++)
    {
        Spins[i] = 0; //initialise all spins to be in the 0 state 
        m0r[i] = 0; //initialise all m0r to be 0
        mr[i] = 0; //initialise all mr to be 0
    }

    for(int s = 1; s < q; s++)
    {
        M[s] = 0; //initialise counters 
    }
    M[0] = N;


    srand((unsigned) time(0)); //seed the random number generator


    //equilibration/burn-in
    for(int t = 0; t < NEsteps; t++)
    {
        for(int c = 0; c < Nclusters; c++)
        {
            FlipAndBuild_from(rand() % N); //flip the spin of a random site and build the cluster 
        }
    }

    //do your measurements
    for(int n = 0; n < Nbins; n++)
    {
        complex<double> m(0., 0.);
        double m1 = 0; double m2 = 0; double m4 = 0;
        complex<double> m0(0., 0.);

        for(int t = 0; t < NMCsteps; t++)
        {
            for(int c = 0; c < Nclusters; c++)
            {
                FlipAndBuild_from(rand() % N); //flip the spin of the site and build the cluster 
            }

            complex<double> tm(0., 0.);
            for(int s = 0; s < q; s++)
            {
                tm += W[s] * double(M[s]); //compute the order parameter
            }

            tm/=N; //normalise the order parameter
            double tm1 = abs(tm); //compute the first moment of the order parameter
            double tm2 = tm1 * tm1; //compute the second moment of the order parameter

            m += tm; m1 += tm1; m2 += tm2; m4 += tm2 * tm2; //compute the moments of the order parameter

            complex<double> m0conj = conj(W[Spins[0]]); //compute the order parameter for each site
            m0 += m0conj; //compute the order parameter for each site
            
            for(int j = 0; j < N; j++)
            {
                m0r[j] += m0conj * W[Spins[j]]; //compute the order parameter for each site
                mr[j] += W[Spins[j]]; //compute the order parameter for each site
            }
            m+= tm; m1 += tm1, m2 += tm2, m4 += tm2 * tm2; //compute the magnetic moments
        }
        
        //normalise everything by number of MC steps 
        m/=NMCsteps; m1/=NMCsteps; m2/=NMCsteps; m4/=NMCsteps; 
        m0 /= NMCsteps; //normalise the order parameter for each site

        for(int j = 0; j < N; j++)
        {
            m0r[j] /= NMCsteps;
            mr[j] /= NMCsteps;

        //print length of m0r for a sanity check

        // cout << m0r.size() << endl;

        //compute the correlation function: i think that this is computing the 'experimental' average 
        for(int j = 0; j < N; j++)
        {
            Cr[j] += m0r[j] - m0 * mr[j]; //do not understand why there is a summation here: it could be for the experimental average 
            // cout << Cr[j] << endl;
        }

        // cout << Cr.size() << endl;

    // cout << endl;

    }

    //  //print size of Cr for each j 

    // for(int j = 0; j < N; j++)
    // {
    //     cout << Cr[j] << endl;
    // } 

    for (int j = 0; j < N; j++)
    {
        Cr[j] /= Nbins;
        std::cout << Cr[j].real() << endl;
    }

    
    ofstream file("correlation_0.25.txt");

        if (file.is_open()) {
        // Write vector values to the file

        for (int j = 0; j < N; j++)
        {
            file << Cr[j].real() << " ";
        }

        // Close the file
        file.close();

        cout << "Vector values saved to file successfully." << endl;
    } else {
        cerr << "Failed to open file for writing." << endl;
    }


    //not ready yet 

    // ofstream file("magnetic_moments.txt");

    //     if (file.is_open()) {
    //     // Write vector values to the file

    //     for (int j = 0; j < N; j++)
    //     {
    //         file << Cr[j].real() << " ";
    //     }

    //     // Close the file
    //     file.close();

    //     cout << "Magnetic moments saved to file successfully." << endl;
    // } else {
    //     cerr << "Failed to open file for writing." << endl;
    // }

    }
}


















