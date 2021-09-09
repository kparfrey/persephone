#include <highfive/H5Easy.hpp>
#include <iostream>

int main()
{
    using std::vector;
    using std::cout;
    using std::endl;


    H5Easy::File dfile("../../DESC/examples/DESC/DSHAPE_output.h5", H5Easy::File::ReadOnly);

    vector<double> R_lmn;

    R_lmn = H5Easy::load<vector<double>>(dfile, "/_equilibria/2/_R_lmn");

    const int num_R_modes = R_lmn.size(); 
    cout << "No. of modes in R function: " << num_R_modes << endl;

    for (int i = 0; i < 10; ++i)
        cout << "Coeff " << i << " :  " << R_lmn[i] << endl;


    cout << "\n" << endl;

    vector<vector<int>> R_modes;

    R_modes = H5Easy::load<vector<vector<int>>>(dfile, "/_equilibria/2/_R_basis/_modes");

    const int n0 = R_modes.size();
    const int n1 = R_modes[0].size();

    cout << "Mode array dimensions: " << n0 << " " << n1 << endl;
    for (int i = 0; i < n0; ++i)
        cout << "Coeff " << i << " :  " << R_modes[i][0] << " " << R_modes[i][1] << " " <<  R_modes[i][2] << endl;



    /* Count the number of existing equilibria in the file */
    int n = 0;
    bool found = false;
    while(!found)
    {
        try
        {
            H5Easy::load<int>(dfile,"/_equilibria/"+std::to_string(n)+"/_L");
        }
        catch(...)
        {
            break;
        }
        n++;
    }

    cout << "Max index of equilibria: " << n-1 << endl;



    return 1;
}
