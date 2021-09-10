#include "desc.hpp"

#include <highfive/H5Easy.hpp>
#include <cmath>
#include <iostream>
#include "write_screen.hpp"


DescConfig::DescConfig(std::string input_file, const int iteration)
{
    using std::vector;

    write::message("\nReading DESC data from " + input_file + ", iteration " + std::to_string(iteration));

    H5Easy::File data(input_file, H5Easy::File::ReadOnly);

    std::string base = "/_equilibria/" + std::to_string(iteration) + "/";

    L = H5Easy::load<int>(data, base + "_L");
    M = H5Easy::load<int>(data, base + "_M");
    N = H5Easy::load<int>(data, base + "_N");

    L_lmn = H5Easy::load<vector<double>>(data, base + "_L_lmn");
    R_lmn = H5Easy::load<vector<double>>(data, base + "_R_lmn");
    Z_lmn = H5Easy::load<vector<double>>(data, base + "_Z_lmn");

    L_modes = H5Easy::load<vector<vector<int>>>(data, base + "_L_basis/_modes");
    R_modes = H5Easy::load<vector<vector<int>>>(data, base + "_R_basis/_modes");
    Z_modes = H5Easy::load<vector<vector<int>>>(data, base + "_Z_basis/_modes");

    iota     = H5Easy::load<vector<double>>(data, base + "_pressure/_params");
    pressure = H5Easy::load<vector<double>>(data, base + "_iota/_params");

    N_L = L_lmn.size();
    N_R = R_lmn.size();
    N_Z = Z_lmn.size();

    N_iota     = iota.size();
    N_pressure = pressure.size();

    write::variable<int>("Radial (lambda) polynomial resolution    - L ", L);
    write::variable<int>("Surface polynomials, poloidal resolution - M ", M);
    write::variable<int>("Surface polynomials, toroidal resolution - N ", N);
    write::variable<int>("No. of modes in lambda polynomial ", N_L);
    write::variable<int>("No. of modes in R polynomial      ", N_R);
    write::variable<int>("No. of modes in Z polynomial      ", N_Z);
    write::variable<int>("No. of modes in iota polynomial   ", N_iota);
    write::variable<int>("No. of modes in pressure poly.    ", N_pressure);

    write::message("Finished loading input data\n");

    return;
}
