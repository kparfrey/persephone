#ifndef INITIAL_STATE_TORUS_HPP
#define INITIAL_STATE_TORUS_HPP

#include "common.hpp"

class ElementBlock;
class CerfonFreidbergConfig;

void set_euler_torus(ElementBlock& eb);
void set_CerfonFreidberg(ElementBlock& eb, CerfonFreidbergConfig& cf_config);
void set_uniform(ElementBlock& eb);

#endif
