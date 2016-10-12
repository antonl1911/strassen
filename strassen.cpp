#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;
#include "matrix.hpp"

typedef MatrixVA<double> Matrix;

// Set LEAF_SIZE to 1 if you want to the pure strassen algorithm
// otherwise, the ikj-algorithm will be applied when the split
// matrices are as small as LEAF_SIZE x LEAF_SIZE
int leafsize = 64;


/*
 * Implementation of the strassen algorithm, similar to 
 * http://en.wikipedia.org/w/index.php?title=Strassen_algorithm&oldid=498910018#Source_code_of_the_Strassen_algorithm_in_C_language
 */

void printMatrix(Matrix matrix, int n);
void read(string filename, Matrix &A, Matrix &B);

Matrix strassenR(Matrix &A, Matrix &B, const int tam)
{
    if (tam <= leafsize)
        return A*B;

    // other cases are treated here:
    int newTam = tam/2;
    Matrix a11(newTam, A), a12(newTam), a21(newTam), a22(newTam),
	   b11(newTam, B), b12(newTam), b21(newTam), b22(newTam),
	    c12(newTam), c21(newTam), c22(newTam);

    // dividing the matrices in 4 sub-matrices:
    for (int i = 0; i < newTam; i++) {
	for (int j = 0; j < newTam; j++) {
	    a12(i, j) = A(i, j + newTam);
	    a21(i, j) = A(i + newTam, j);
	    a22(i, j) = A(i + newTam, j + newTam);

	    b12(i, j) = B(i, j + newTam);
	    b21(i, j) = B(i + newTam, j);
	    b22(i, j) = B(i + newTam, j + newTam);
	}
    }

    // Calculating p1 to p7:
    Matrix aResult = Matrix(a11 + a22);
    Matrix bResult = Matrix(b11 + b22);
    auto p1 = strassenR(aResult, bResult, newTam); // p1 = (a11+a22) * (b11+b22)
    aResult = a21 + a22;
    auto p2 = strassenR(aResult, b11, newTam);	// p2 = (a21+a22) * (b11)
    bResult = b12 - b22;
    auto p3 = strassenR(a11, bResult, newTam);	// p3 = (a11) * (b12 - b22)
    bResult = b21 - b11;
    auto p4 = strassenR(a22, bResult, newTam);	// p4 = (a22) * (b21 - b11)
    aResult = a11 + a12;
    auto p5 = strassenR(aResult, b22, newTam);	// p5 = (a11+a12) * (b22)
    aResult = a21 - a11;
    bResult = b11 + b12;
    auto p6 = strassenR(aResult, bResult, newTam); // p6 = (a21-a11) * (b11+b12)
    aResult = a12 - a22;
    bResult = b21 + b22;
    auto p7 = strassenR(aResult, bResult, newTam); // p7 = (a12-a22) * (b21+b22)

    // calculating c21, c21, c11 e c22:

    c12 = p3 + p5; // c12 = p3 + p5
    c21 = p2 + p4; // c21 = p2 + p4
    c22 = p1 + p6 + (p3 - p2); // c22 = p1 + p3 - p2 + p6
    // Grouping the results obtained in a single matrix:
    Matrix C = Matrix(tam, p1 + p7 + p4 - p5); // c11 = p1 + p4 - p5 + p7
    for (int i = 0; i < newTam ; i++) {
	for (int j = 0 ; j < newTam ; j++) {
	    C(i, j + newTam) = c12(i, j);
	    C(i + newTam, j) = c21(i, j);
	    C(i + newTam, j + newTam) = c22(i, j);
	}
    }
    return C;
}

unsigned int nextPowerOfTwo(int n) {
    return pow(2, int(ceil(log2(n))));
}

Matrix strassen(Matrix &A, Matrix &B, int n)
{
    unsigned int m = nextPowerOfTwo(n);
    Matrix APrep(m, A), BPrep(m, B);

    return Matrix(n, strassenR(APrep, BPrep, m));
}

int getMatrixSize(string filename)
{
    string line;
    ifstream infile;
    infile.open (filename.c_str());
    getline(infile, line);
    return count(line.begin(), line.end(), '\t') + 1;
}

void read(string filename, Matrix &A, Matrix &B)
{
    string line;
    FILE* matrixfile = freopen(filename.c_str(), "r", stdin);

    if (matrixfile == 0) {
        cerr << "Could not read file " << filename << endl;
        return;
    }

    int i = 0, j;
    double a;
    while (getline(cin, line) && !line.empty()) {
        istringstream iss(line);
        j = 0;
        while (iss >> a) {
            A(i, j) = a;
            j++;
        }
        i++;
    }

    i = 0;
    while (getline(cin, line)) {
        istringstream iss(line);
        j = 0;
        while (iss >> a) {
            B(i, j) = a;
            j++;
        }
        i++;
    }

    fclose (matrixfile);
}

void printMatrix(Matrix matrix, int n)
{
    for (int i = 0; i < n; i++) {
        for (int j=0; j < n; j++) {
            if (j != 0) {
                cout << "\t";
            }
            cout << matrix(i,j);
        }
        cout << endl;
    }
}

int main (int argc, char* argv[])
{
    string filename;
    if (argc < 3) {
        filename = "2000.in";
    } else {
        filename = argv[2];
    }

    if (argc < 5) {
        leafsize = 64;
    } else {
        leafsize = atoi(argv[4]);
    }

    int n = getMatrixSize(filename);
    Matrix A(n), B(n);
    read (filename, A, B);
    Matrix C = strassen(A, B, n);
    printMatrix(C, n);
    return 0;
}
