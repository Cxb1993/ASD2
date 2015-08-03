#include "mac2d.h"

int GNx = 0, GNy = 0, GNumBCGrid = 0;
int idx3(int nx, int i, int j);

int main(int argc, char *argv[]) {
    
    return 0;
}


inline int idx3(int nx, int i, int j) {
    return i + (GNx + 2 * GNumBCGrid) * j;
}