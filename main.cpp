/*
Program to solve the two-dimensional Ising model
The coupling constant J = 1
Boltzmann’s constant = 1, temperature has thus dimension energy
Metropolis sampling is used. Periodic boundary conditions.
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
using namespace std;
ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}

// Function to read in data from screen
void read_input(int& n_spins, int& mcs, double& initial_temp, double& final_temp, double& temp_step){
    cout << "n_spins:";
    cin >> n_spins;

    cout << "mcs:";
    cin >> mcs;

    cout << "initial_temp:";
    cin >> initial_temp;

    cout << "final_temp:";
    cin >> final_temp;

    cout << "temp_step:";
    cin >> temp_step;
};


// Function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, double temp, int **spin_matrix, double& E, double& M)
{
    // setup spin matrix and intial magnetization
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
        if (temp < 1.5) spin_matrix[y][x] = 1; // spin orientation for the ground state
        M += (double) spin_matrix[y][x];
        }
    }
    // setup initial energy
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
        E -= (double) spin_matrix[y][x]*
        (spin_matrix[periodic(y,n_spins,-1)][x] +
        spin_matrix[y][periodic(x,n_spins,-1)]);
        }
    }
}// end function initialise

// The Metropolis algorithm
void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w)
{
    // loop over all spins
        for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            // Find random position
            int ix = (int) (ran1(&idum)*(double)n_spins);
            int iy = (int) (ran1(&idum)*(double)n_spins);

            int deltaE = 2*spin_matrix[iy][ix]*
            (spin_matrix[iy][periodic(ix,n_spins,-1)]+
            spin_matrix[periodic(iy,n_spins,-1)][ix] +
            spin_matrix[iy][periodic(ix,n_spins,1)] +
            spin_matrix[periodic(iy,n_spins,1)][ix]);

            // Here we perform the Metropolis test
                if ( ran1(&idum) <= w[deltaE+8] ) {
                spin_matrix[iy][ix] *= -1;
                // flip one spin and accept new spin config
                // update energy and magnetization
                M += (double) 2*spin_matrix[iy][ix];
                E += (double) deltaE;
            }
        }
    }
} // end of Metropolis sampling over spins

void output(int n_spins, int mcs, double temp, double *average){
    cout << "n_spins:" << n_spins << endl;
    cout << "mcs:" << mcs << endl;
    cout << "temp:" << temp << endl;
    cout << "average:" << average << endl;
};



// main program
int main(int argc, char* argv[])
{
    char *outfilename;
    long idum;
    int **spin_matrix, n_spins, mcs;
    double w[17], average[5], initial_temp, final_temp, E, M, temp_step;
    // Read in output file, abort if there are too few command-line arguments
    if( argc <= 1 ){
        cout << "Bad Usage: " << argv[0] <<
        " read also output file on same line" << endl;
        exit(1);
    }
    else{
        outfilename=argv[1];
    }
ofile.open(outfilename);

// Read in initial values such as size of lattice, temp and cycles
read_input(n_spins, mcs, initial_temp, final_temp, temp_step);

spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));

idum = -1; // random starting point

for (double temp = initial_temp; temp <= final_temp; temp+=temp_step){
    // initialise energy and magnetization
    E = M = 0.;

    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temp);

    // initialise array for expectation values
    for( int i = 0; i < 5; i++) average[i] = 0.;
    initialize(n_spins, temp, spin_matrix, E, M);

    // start Monte Carlo computation
    for (int cycles = 1; cycles <= mcs; cycles++){
        Metropolis(n_spins, idum, spin_matrix, E, M, w);
        // update expectation values
        average[0] += E; average[1] += E*E;
        average[2] += M; average[3] += M*M; average[4] += fabs(M);
    }
// print results
    output(n_spins, mcs, temp, average);
}

free_matrix((void **) spin_matrix); // free memory
ofile.close(); // close output file
return 0;
}
