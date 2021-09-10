#ifndef WRITE_SCREEN_HPP
#define WRITE_SCREEN_HPP

#include <string>
#include "common.hpp"

/* Move iostream to the cpp file. Currently it is needed to supply cout etc.
 * to a few other files, but this is sloppy. */
#include <iostream> 

namespace write
{
    using std::string;

    void message(string message, string prefix="", 
                                        bool set_root=false, int rank=-1);

    void error(string error, Prognosis prognosis=destroy, 
                                        bool set_rank=false, int rank_=-1);

    template<typename type>
    void variable(string message, type variable, WhoWrites who_writes=root_only,
                                        bool set_rank=false, int rank_=-1);

    void basic(string message);

    void store_rank(int rank);

    string str(int n);
    string str(double n);
    string str(float n);
}
#endif
