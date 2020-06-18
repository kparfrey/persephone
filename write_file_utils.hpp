#ifndef WRITE_FILE_UTILS_HPP
#define WRITE_FILE_UTILS_HPP

#include <highfive/H5File.hpp>
#include "common.hpp"
#include "element_block.hpp"

using std::string;

void repack(real_t* original, real_t* eb_logical, ElementBlock& eb);

void vector_organization(HighFive::File file,
                         const int Nvec, int& count, string group,
                         string vecnames[], VectorField veclist[],
                         string names[], real_t* datalist[]);
#endif
