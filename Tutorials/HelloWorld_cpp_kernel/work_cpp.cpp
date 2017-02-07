#include <Box.H>
#include <REAL.H>
#include <cmath>

void work_cpp(const Box* bx,
          const int ng,
          Real *dataptr) {

    int ijk;

    // When we flatten the (i,j,k) indices into a single integer, we need to
    // know how far to stride when increment j or k (the i-stride is 1). We
    // must pad these stride lengths by the number of ghost cells on each end
    // of the axis.
    const int BL_jStride = bx->length(0) + 2*ng;
    const int BL_kStride = BL_jStride * (bx->length(1) + 2*ng);

    // The "absolute" coordinates of the box on the problem domain.
    const int *lo = bx->loVect();
    const int *hi = bx->hiVect();

    // Used to convert relative C-style grid indices (which always start at 0)
    // to the absolute grid coordinates.
    int abs_i, abs_j, abs_k;

    for (int k = 0; k < bx->length(2); ++k) {
        for (int j = 0; j < bx->length(1); ++j) {
            for (int i = 0; i < bx->length(0); ++i) {

                abs_i = lo[0] + i;
                abs_j = lo[1] + j;
                abs_k = lo[2] + k;

                // Compute the flattened 1-D index from the 3-D indices.
                ijk = (i + ng) + (j + ng)*BL_jStride + (k + ng)*BL_kStride;

                dataptr[ijk] = exp(-(double)(abs_i*abs_i + abs_j*abs_j + abs_k*abs_k)/512.0);
            }
        }
    }
}
