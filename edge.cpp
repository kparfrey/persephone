#include "edge.hpp"
#include "lagrange_polynomials.hpp"


void Edge::setup(const int Ns, const real_t* const __restrict__ xs)
{
    N = Ns;

    x = new real_t [N];

    for (int i: dirs)
        r(i) = new real_t [N];

    for (int i = 0; i < N; ++i)
        x[i] = xs[i];

    w = lagrange::barycentric_weights(x, N);

    return;
}


/* Evaluate the coordinate vector r at location s, store in rs */
void Edge::eval(const real_t s, real_t* const __restrict__ rs)
{
    /* Maybe slightly faster? */
    //lagrange::interpolate_vector(x, w, r, s, N, rs);

    for (int d: dirs)
        rs[d] = lagrange::interpolate(x, w, r(d), s, N);

    return;
}


void Edge::diff(const real_t s, real_t* const __restrict__ drdx)
{
    for (int d: dirs)
        drdx[d] = lagrange::differentiate(x, w, r(d), s, N);

    return;
}
