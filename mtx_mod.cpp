#include "El.hpp"   // Elemental library
#include <cstdlib>  // Standard import
#include <string>   // Text for file names, etc.

#include <ctime>

#define NOINVERT
// Redefine namespace and working set for ease of use
using namespace std;
using namespace El;

typedef double Real;

/**
 * mtx_mod performs the compound matrix operation (A^T * A)^-1 * A^T.
 * Params: --filein  - Input file name, default 'mtx.txt'
 *         --fileout - Output file name, defaults to the input filename
 *         --display - 'true' allows intermediate output; 'false' 
 *                     supresses it.  Default 'false'.
 */
int main(int argc, char* argv[])
{
    // Initialize Elemental and MPI
    Initialize(argc, argv);
 
    // We have an even number of arguments iff there is a lambda argument
    float lambda = 0.0f;
    if (argc % 2 == 0) {
        lambda = stof(argv[1]);
    }
    try
    {
        // Command line arguments are optional here, but highly recommended for clarity.
        const string input  = Input("--filein",    "Input file name", "mtx.mat");
        const string output = Input("--fileout",   "Output file name", input);
        const int height    = Input("--h",         "Input matrix height", 1);
        const int width     = Input("--w",         "Input matrix width",  1);
        bool display  = Input("--display",   "Display temporary matrices", false);
        ProcessInput();
        display = display && (mpi::WorldRank() == 0);
        if (display) PrintInputReport();

        if (display)
        {
            cout << "Grid is " << DefaultGrid().Height() << " x " << DefaultGrid().Width() << ", "
                 << DefaultGrid().Size() << " active processes\n";
            cout << "Matrix is " << height << " x " << width << endl;
            cout << "Lambda is " << lambda << endl;
        }
        
        // Initialize the matrices over the processor grid.
        DistMatrix<Real> mtx(height, width, DefaultGrid());
        DistMatrix<Real> identity(DefaultGrid());
        DistMatrix<Real> modified(DefaultGrid());

        // Read the initial matrix from the file.
        if (display) cout << "Reading...\n";
        Read(mtx, input, BINARY_FLAT);
        if (display) Print(mtx, "A");

        time_t start;
        if (mpi::WorldRank() == 0) start = time(0);

#ifndef INVERT
        // Transpose the original matrix
        Transpose(mtx, modified);
        if (display) Print(modified, "A^T");
#endif

        // Perform A^T * A, transposing A implicitly.
        // (Note: gemm is the Level 3 BLAS matrix multiplication 
        // function.)
        if (lambda) {
            Identity(identity, width, width);
            Gemm(TRANSPOSE, NORMAL, Real(1.0), mtx, mtx, Real(lambda), identity);
        } else {
            Gemm(TRANSPOSE, NORMAL, Real(1.0), mtx, mtx, identity);
        }
        if (display) Print(identity, "A^T*A");

// If options specify to use a full inversion rather than a linear solution, do that
#ifdef INVERT

        // Invert the result in place.  Note that 'product' is 
        // guaranteed to be symmetric.
        SymmetricInverse(LOWER, identity);
        if (display) Print(identity, "(A^T*A)^-1");

        // Perform the final product, again transposing A implicitly.
        Gemm(NORMAL, TRANSPOSE, Real(1.0), identity, mtx, modified);
        if (display) Print(modified, "(A^T*A)^-1*A^T");

// Otherwise, solve a linear system
#else
        // Solve the system of linear equations (A^T*A)X = A^T  Note
        // that 'product' is guaranteed to be symmetric.
        SymmetricSolve(LOWER, NORMAL, mtx, modified);
        if (display) Print(modified, "(A^T*A)^-1*A^T");

#endif

        if (mpi::WorldRank() == 0) cout << "Time: " <<
            (time(0) - start) << endl;

        // Write the result to a file
        Write(modified, output, BINARY_FLAT);
    }
        catch (exception &e)
    { 
        // The result of an error is that a file is not written.
    }
    // End Elemental and MPI
    Finalize();

    return 0;
}
