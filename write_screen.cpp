#include "write_screen.hpp"

/* This way of using static variables in functions is a bit ugly, since it
 * requires these function arguments that are only used once, at setup.
 * Should convert these functions into functors, with "set_rank()" etc.
 * member functions. Would be cleaner. 
 *
 * On the other hand, this approach means the functions are immediately 
 * available in any file which includes this one --- it doesn't need to have
 * access to an instance of a functor class.... 
 *
 * Maybe this can survive as an example of using function-static variables
 * and default arguments... */

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
            cout << message << ": " << variable << endl;
            //cout << "Data on root: " << message << ": " << variable << endl;

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
        write::message("Setting root process in write::message()", "", true, rank);
        write::error("Setting process rank in write::error()", survive, true, rank);

        real_t dummy_real = 0.0;
        int    dummy_int  = 0;
        write::variable<real_t>("Setting rank for real data in write::variable()", 
                                                            dummy_real,root_only, true, rank);
        write::variable<int>("Setting rank for int data in write::variable()", 
                                                            dummy_int, root_only, true, rank);

        return;
    }


    string str(int n)
    {
        return std::to_string(n);
    }

    string str(double n)
    {
        return std::to_string(n);
    }

    string str(float n)
    {
        return std::to_string(n);
    }
}

