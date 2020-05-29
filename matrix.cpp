#include "matrix.hpp"
#include "write_screen.hpp"


Matrix::Matrix()
{
}


/* Need to create a new object like this and set storeArray=true
 * if you want to use the m[3][3] array, eg in test(). */
Matrix::Matrix(real_t (&m)[3][3], bool storeArray)
{
    fill(m);

    if (storeArray)
    {
        for (int i: dirs)
            for (int j: dirs)
                arr[i][j] = m[i][j];
    }
}


void Matrix::fill(real_t (&m)[3][3])
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

    inv[0][0] = A;
    inv[0][1] = D;
    inv[0][2] = G;
    inv[1][0] = B;
    inv[1][1] = E;
    inv[1][2] = H;
    inv[2][0] = C;
    inv[2][1] = F;
    inv[2][2] = I;

    for (int i: dirs)
        for (int j: dirs)
            inv[i][j] /= det;

    return;
}


void Matrix::test()
{
    find_inverse();

    write::variable<real_t>("Determinant", det);

    real_t ident[3][3] = {}; // Zero-initialize

    for (int a: dirs)
    for (int b: dirs)
        for (int i: dirs)
            ident[a][b] += inv[a][i] * arr[i][b];
        
    for (int a: dirs)
    for (int b: dirs)
        write::variable<real_t>(std::to_string(a)+" "+std::to_string(b), ident[a][b]);

    return;
}
