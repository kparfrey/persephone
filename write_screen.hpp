#ifndef WRITE_SCREEN_HPP
#define WRITE_SCREEN_HPP

#include <iostream>
#include <string>
#include "common.hpp"

namespace write
{
    using std::string;

    void message(string message, string prefix="Progress: ", 
                                        bool set_root=false, int rank=-1);

    void error(string error, Prognosis prognosis=destroy, 
                                        bool set_rank=false, int rank_=-1);

    template<typename type>
    void variable(string message, type variable, WhoWrites who_writes=root_only,
                                        bool set_rank=false, int rank_=-1);

    void basic(string message);

    void store_rank(int rank);
}
#endif
