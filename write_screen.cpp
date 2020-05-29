#include "write_screen.hpp"

namespace write
{
    using std::cout;
    using std::endl;


    void message(string message, string prefix, bool set_root, int rank)
    {
        static bool am_root;

        if (set_root)
        {
            if (rank < 0)
            {
                cout << "Error: must set process rank in write::message()." << endl;
                exit(1);
            }

            if (rank == 0)
                am_root = true;
            else
                am_root = false;
        }

        if (am_root)
            cout << prefix << message << endl;

        return;
    }


    void error(string error, Prognosis prognosis, bool set_rank, int rank_)
    {
        static int rank;

        if (set_rank)
        {
            if (rank < 0)
            {
                cout << "Error: must set process rank in write::error()." << endl;
                exit(1);
            }
            
            rank = rank_;

            if (rank == 0)
                cout << error << endl;
            return;
        }

        cout << "Error --- rank: " << rank << ": " << error << endl;

        if (prognosis == destroy)
            exit(1);

        return;
    }


    template<typename type>
    void variable(string message, type variable, WhoWrites who_writes, bool set_rank, int rank_)
    {
        static bool am_root;
        static int rank;

        if (set_rank)
        {
            if (rank < 0)
            {
                cout << "Error: must set process rank in write::variable()." << endl;
                exit(1);
            }
            
            rank = rank_;

            if (rank == 0)
                am_root = true;
            else
                am_root = false;
            
            if (am_root)
                cout << message << endl;

            return;
        }


        if ((who_writes == root_only) && am_root)
            cout << "Data on root: " << message << ": " << variable << endl;

        if (who_writes == all_ranks)
            cout << "Data on rank " << rank << ": " << message << ": " << variable << endl;

        return;
    }
    
    /* Explicit instantiation of the template for the types we'll use it for */
    template void variable<int>(string message, int variable, WhoWrites who_writes,
                                                     bool set_rank, int rank_);
    template void variable<real_t>(string message, real_t variable, WhoWrites who_writes,
                                                     bool set_rank, int rank_);


    void basic(string message)
    {
        cout << message << endl;
        return;
    }



    void store_rank(int rank)
    {
        write::message("Setting root process in write::message()", "Progress: ",
                                                                      true, rank);
        write::error("Progress: Setting process rank in write::error()", survive, true, rank);

        real_t dummy_real = 0.0;
        int    dummy_int  = 0;
        write::variable<real_t>("Progress: Setting rank for real data in write::variable()", 
                                                            dummy_real,root_only, true, rank);
        write::variable<int>("Progress: Setting rank for int data in write::variable()", 
                                                            dummy_int, root_only, true, rank);

        return;
    }

}

