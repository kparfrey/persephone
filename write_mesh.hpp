#ifndef WRITE_MESH_HPP
#define WRITE_MESH_HPP

#include <vector>
#include <mpi.h>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

#include "common.hpp"
#include "process.hpp"
#include "element_block.hpp"


void write_mesh(Process &proc);

#endif
