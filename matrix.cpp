#include "matrix.hpp"

Matrix::Matrix(real_t (&m)[3][3])
    : m(m)
{
    a = m[0][0];
    b = m[0][1];
    c = m[0][2];
    d = m[1][0];
    e = m[1][1];
    f = m[1][2];
    g = m[2][0];
    h = m[2][1];
    i = m[2][2];

    cofactors_found   = false;
    determinant_found = false;
}


void Matrix::find_cofactors()
{
    A = e*i - f*h;
    B = f*g - d*i;
    C = d*h - e*g;
    D = c*h - b*i;
    E = a*i - c*g;
    F = b*g - a*h;
    G = b*f - c*e;
    H = c*d - a*f;
    I = a*e - b*d;

    cofactors_found = true;

    return;
}


void Matrix::find_determinant()
{
    if (!cofactors_found)
    {
        A = e*i - f*h;
        B = f*g - d*i;
        C = d*h - e*g;
    }

    det = a*A + b*B + c*C;

    determinant_found = true;

    return;
}


void Matrix::find_inverse()
{
    if (!cofactors_found)
        find_cofactors();

    if (!determinant_found)
        find_determinant();

    minv[0][0] = A;
    minv[0][1] = D;
    minv[0][2] = G;
    minv[1][0] = B;
    minv[1][1] = E;
    minv[1][2] = H;
    minv[2][0] = C;
    minv[2][1] = F;
    minv[2][2] = I;

    for (int i: dirs)
        for (int j: dirs)
            minv[i][j] /= det;

    return;
}

void Matrix::test()
{
    find_inverse();

    std::cout << "Determinant: " << det << std::endl;

    real_t ident[3][3] = {}; // Zero-initialize

    for (int a: dirs)
    for (int b: dirs)
        for (int i: dirs)
            ident[a][b] += minv[a][i] * m[i][b];
        
    for (int a: dirs)
    for (int b: dirs)
        std::cout<<a<<" "<<b<<"   "<<ident[a][b]<<std::endl;

    return;
}
