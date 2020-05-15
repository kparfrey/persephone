#ifndef EDGE_HPP
#define EDGE_HPP

#include "common.hpp"

/* Just a placeholder class for now for generally curved edge object */
/* Think you should only need real Edge objects around processes/element blocks
 * --- within a block you can find all the physical locations and derivatives
 *     by linear blending inwards from the outer edges. */
class Edge
{
    public:

    int temp_data;
};

#endif
