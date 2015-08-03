/*
 Trouble Shooting
 - unknown file: error: SEH exception with code 0xc00000005
                        thrown in the test body
    -> Assign NULL to your variable as a initial value.
*/
#include "data.h"

Data::Data() {
}

Data::~Data() {   
}

int* Data::Allocate1Di(int nx) {
    int *array = NULL;

    array = (int *)malloc(nx * sizeof(int));

    return array;
}

float* Data::Allocate1Df(int nx) {
    float *array = NULL;

    array = (float *)malloc(nx * sizeof(float));

    return array;
}

double* Data::Allocate1Dd(int nx) {
    double *array = NULL;

    array = (double *)malloc(nx * sizeof(double));

    return array;
}

long double* Data::Allocate1Dld(int nx) {
    long double *array = NULL;

    array = (long double *)malloc(nx * sizeof(long double));

    return array;
}

int Data::Deallocate1Di(int *to_deallocate) {
    free(to_deallocate);

    return 0;
}

int Data::Deallocate1Df(float *to_deallocate) {
    free(to_deallocate);

    return 0;
}

int Data::Deallocate1Dd(double *to_deallocate) {
    free(to_deallocate);

    return 0;
}

int Data::Deallocate1Dld(long double *to_deallocate) {
    free(to_deallocate);

    return 0;
}

int** Data::Allocate2Di(int nx, int ny) {
    int **array = NULL;
    
    array = (int **) malloc(nx * sizeof(int *));

    array[0] = (int *) malloc(nx * ny * sizeof(int));

    // Make contiguous array due to HDF5
    for (int i = 1; i < nx; ++i) {
        array[i] = array[i-1] + ny;
    }

    // Initialize
    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j) {
        array[i][j] = 0;
    }
    
    return array;
}

float** Data::Allocate2Df(int nx, int ny) {
    float **array = NULL;
    array = (float **) malloc(nx * sizeof(float *));

    array[0] = (float *) malloc(nx * ny * sizeof(float));

    // Make contiguous array due to HDF5
    for (int i = 1;  i < nx;  ++i) {
        array[i] = array[i-1] + ny;
    }

    // Initialize
    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j) {
        array[i][j] = 0.0;
    }

    return array;
}

double** Data::Allocate2Dd(int nx, int ny) {
    double **array = NULL;
    array = (double **) malloc(nx * sizeof(double *));

    array[0] = (double *) malloc(nx * ny * sizeof(double));

    // Make contiguous array due to HDF5
    for (int i = 1; i < nx; ++i ) {
        array[i] = array[i-1] + ny;
    }
    
    // Initialize
    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j) {
        array[i][j] = 0.0;
    }

    return array;
}

long double** Data::Allocate2Dld(int nx, int ny) {
    long double **array = NULL;
    array = (long double **) malloc(nx * sizeof(long double *));

    array[0] = (long double *) malloc(nx * ny * sizeof(long double));

    // Make contiguous array due to HDF5
    for (int i = 1; i < nx; ++i) {
        array[i] = array[i-1] + ny;
    }
    
    // Initialize
    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j) {
        array[i][j] = 0.0;
    }
    return array;
}

int Data::Deallocate2Di(int **array) {
    // allocated twice, free twice
    free(array[0]);
    free(array);

    return 0;
}

int Data::Deallocate2Df(float **array) {
    // allocated twice, free twice
    free(array[0]);
    free(array);

    return 0;
}

int Data::Deallocate2Dd(double **array) {
    // allocated twice, free twice
    free(array[0]);
    free(array);

    return 0;
}

int Data::Deallocate2Dld(long double **array) {
    // allocated twice, free twice
    free(array[0]);
    free(array);

    return 0;
}

int*** Data::Allocate3Di(int nx, int ny, int nz) {
    int ***array = NULL;
    array = (int ***) malloc(nx * sizeof(int **));

    array[0] = (int **) malloc(nx * ny * sizeof(int *));

    array[0][0] = (int *) malloc(nx * ny * nz * sizeof(int));

    // Make contiguous array due to HDF5
    for (int j = 1; j < ny; ++j ) {
        array[0][j] = array[0][j-1] + nz;
    }
    for (int i = 1; i < nx; ++i ) {
        array[i] = array[i-1] + ny;

        array[i][0] = array[i-1][ny-1] + nz;
        for (int j = 1; j < ny; ++j) {
            array[i][j] = array[i][j-1] + nz;
        }
    }
    
    // Initialize
    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    for (int k = 0; k < nz; ++k) {
        array[i][j][k] = 0.0;
    }

    return array;
}

float*** Data::Allocate3Df(int nx, int ny, int nz) {
    float ***array = NULL;
    array = (float ***) malloc(nx * sizeof(float **));

    array[0] = (float **) malloc(nx * ny * sizeof(float *));

    array[0][0] = (float *) malloc(nx * ny * nz * sizeof(float));

    // Make contiguous array due to HDF5
    for (int j = 1; j < ny; ++j ) {
        array[0][j] = array[0][j-1] + nz;
    }
    for (int i = 1; i < nx; ++i ) {
        array[i] = array[i-1] + ny;

        array[i][0] = array[i-1][ny-1] + nz;
        for (int j = 1; j < ny; ++j) {
            array[i][j] = array[i][j-1] + nz;
        }
    }

    // Initialize
    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    for (int k = 0; k < nz; ++k) {
        array[i][j][k] = 0;
    }

    return array;
}

double*** Data::Allocate3Dd(int nx, int ny, int nz) {
    double ***array = NULL;
    array = (double ***) malloc(nx * sizeof(double **));

    array[0] = (double **) malloc(nx * ny * sizeof(double *));

    array[0][0] = (double *) malloc(nx * ny * nz * sizeof(double));

    // Make contiguous array due to HDF5
    for (int j = 1; j < ny; ++j ) {
        array[0][j] = array[0][j-1] + nz;
    }
    for (int i = 1; i < nx; ++i ) {
        array[i] = array[i-1] + ny;

        array[i][0] = array[i-1][ny-1] + nz;
        for (int j = 1; j < ny; ++j) {
            array[i][j] = array[i][j-1] + nz;
        }
    }

    // Initialize
    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    for (int k = 0; k < nz; ++k) {
        array[i][j][k] = 0.0;
    }

    return array;
}

long double*** Data::Allocate3Dld(int nx, int ny, int nz) {
    long double ***array = NULL;
    array = (long double ***) malloc(nx * sizeof(long double **));

    array[0] = (long double **) malloc(nx * ny * sizeof(long double *));

    array[0][0] = (long double *) malloc(nx * ny * nz * sizeof(long double));

    // Make contiguous array due to HDF5
    for (int j = 1;  j < ny ;  ++j) {
        array[0][j] = array[0][j-1] + nz;
    }
    for (int i = 1;  i < nx;  ++i) {
        array[i] = array[i-1] + ny;

        array[i][0] = array[i-1][ny-1] + nz;
        for (int j = 1; j < ny; ++j) {
            array[i][j] = array[i][j-1] + nz;
        }
    }

    // Initialize
    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    for (int k = 0; k < nz; ++k) {
        array[i][j][k] = 0.0;
    }
    return array;
}

int Data::Deallocate3Di(int ***array) {
    // allocated thrice, free thrice
    free(array[0][0]);
    free(array[0]);
    free(array);

    return 0;
}

int Data::Deallocate3Df(float ***array) {
    // allocated thrice, free thrice
    free(array[0][0]);
    free(array[0]);
    free(array);

    return 0;
}

int Data::Deallocate3Dd(double ***array) {
    // allocated thrice, free thrice
    free(array[0][0]);
    free(array[0]);
    free(array);
    
    return 0;
}

int Data::Deallocate3Dld(long double ***array) {
    free(array[0][0]);
    free(array[0]);
    free(array);
    
    return 0;
}

// see http://stackoverflow.com/questions/12951162/error-stdios-baseios-baseconst-stdios-base-is-private
// for 'std::ios_base &std::ios_base::operator=(const std::ios_base &)" is inaccessible' error
int Data::OpenPLTforRead(string file_name, std::ifstream& in_file) {

    in_file.open(file_name.c_str());

    return 0;
}

int Data::OpenPLTforWrite(string file_name, std::ofstream& out_file) {

    out_file.open(file_name.c_str(), std::ios::out);

    return 0;
}

int Data::OpenPLTforAppend(string file_name, std::ofstream& out_file) {
    
    out_file.open(file_name.c_str(), std::ios::app);

    return 0; 
}

int Data::WritePLTBasic1D(std::ofstream& of_plt, int nx) {
    of_plt << string("TITLE = VELOCITY") << endl;
    of_plt << string("VARIABLES = \"X\" \"U\" \n") << endl;

    return 0;
}

int Data::WritePLTBasic2D(std::ofstream& of_plt, int nx, int ny) {
    of_plt << string("TITLE = VELOCITY") << endl;
    of_plt << string("VARIABLES = \"X\" \"Y\" \"U\" \"V\" \n") << endl;

    return 0;
}

int Data::WritePLTBasic3D(std::ofstream& of_plt, int nx, int ny, int nz) {
    of_plt << string("TITLE = VELOCITY") << endl;
    of_plt << string("VARIABLES = \"X\" \"Y\" \"Z\"\"U\" \"V\" \"W\"\n") << endl;

    return 0;
}


int Data::ClosePLTWrite(ofstream& file) {
    file.close();
    
    return 0;
}

int Data::ClosePLTRead(ifstream& file) {
    file.close();
    
    return 0;
}

/*
template <typename T>
T Data::getNumericValFromJson(jsonxx::Object io, std::string name) {
    T val;
    if (io.has<Number>(name.c_str()))
        val = io.get<Number>(name.c_str());
    else {
        perror(std::string(std::string("JSON Parse Error;") + name).c_str());
        exit(1);
    }

    return val;
}

bool Data::getBooleanValFromJson(jsonxx::Object io, std::string name) {
    bool val;
    if (io.has<Boolean>(name.c_str()))
        val = io.get<Boolean>(name.c_str());
    else {
        perror(std::string(std::string("JSON Parse Error;") + name).c_str());
        exit(1);
    }

    return val;
}

std::string Data::getStringValFromJson(jsonxx::Object io, std::string name) {
    std::string val;
    if (io.has<String>(name.c_str()))
        val = io.get<String>(name.c_str());
    else {
        perror(std::string(std::string("JSON Parse Error;") + name).c_str());
        exit(1);
    }

    return val;
}
*/