#include "edge.hpp"
#include "lagrange_polynomials.hpp"
#include <cmath>

static real_t TOL = 1e-12; // Tolerance for saying two points are identical

void Edge::setup(const int edge_id, const int Nf[3], const VectorField xf)
{
    id = edge_id;

    offset[0] = 0.0;
    offset[1] = 0.0;
    offset[2] = 0.0;

    switch (id)
    {
        case 0:
            dir = 0;
            break;
        case 1:
            dir = 1;
            offset[0] = 1.0;
            break;
        case 2:
            dir = 0;
            offset[1] = 1.0;
            break;
        case 3:
            dir = 1;
            break;
        case 4:
            dir = 0;
            offset[2] = 1.0;
            break;
        case 5:
            dir = 1;
            offset[0] = 1.0;
            offset[2] = 1.0;
            break;
        case 6:
            dir = 0;
            offset[1] = 1.0;
            offset[2] = 1.0;
            break;
        case 7:
            dir = 1;
            offset[2] = 1.0;
            break;
        case 8:
            dir = 2;
            break;
        case 9:
            dir = 2;
            offset[0] = 1.0;
            break;
        case 10:
            dir = 2;
            offset[0] = 1.0;
            offset[1] = 1.0;
            break;
        case 11:
            dir = 2;
            offset[1] = 1.0;
            break;
    }

    N = Nf[dir];

    x = new real_t [N];

    for (int i: dirs)
        r(i) = new real_t [N];

    for (int i = 0; i < N; ++i)
        x[i] = xf(dir, i);

    w = lagrange::barycentric_weights(x, N);

    return;
}


void Edge::free()
{
    delete[] x;
    delete[] w;

    for (int d: dirs)
        delete[] r(d);

    return;
}


/* Evaluate the coordinate vector r at location s, store in rs */
void Edge::eval(const real_t s, real_t* const __restrict__ rs) const
{
    /* Maybe slightly faster? */
    //lagrange::interpolate_vector(x, w, r, s, N, rs);

    /* Check to see if at an edge point. If so, don't need to interpolate,
     * and the interpolation function will give a NaN. */
    for (int i = 0; i < N; ++i)
        if (std::abs(s - x[i]) < TOL)
        {
            for (int d: dirs)
                rs[d] = r(d,i);
            return;
        }

    for (int d: dirs)
        rs[d] = lagrange::interpolate(x, w, r(d), s, N);

    return;
}


void Edge::diff(const real_t s, real_t* const __restrict__ drdx) const
{
    /* If at an interpolation node, use the at-node version of the
     * differentiation algorithm */
    for (int i = 0; i < N; ++i)
        if (std::abs(s - x[i]) < TOL)
        {
            for (int d: dirs)
                drdx[d] = lagrange::diff_at_node(x, w, r(d), i, N);
            return;
        }

    /* If didn't hit a node, use the in-between-nodes version */
    for (int d: dirs)
        drdx[d] = lagrange::differentiate(x, w, r(d), s, N);

    return;
}
