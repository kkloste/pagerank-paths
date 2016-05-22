/**
 * @file sparselist.hpp
 * A set of functions to maintain a fast hash_table for unsigned integers
 */

#ifndef _SPARSELIST_H_
#define _SPARSELIST_H_

#include <vector>
#include <stdlib.h>
#include <limits>
#include "sparsehash/dense_hash_map.h"

#include "mex.h"

struct sparselist {
    typedef google::dense_hash_map<mwIndex,size_t> map_type;
    map_type map;
    mwIndex lastval; // = std::numeric_limits<local_index_type>::max();      
    sparselist()
    {
    	lastval = std::numeric_limits<mwIndex>::max();
	    map.set_empty_key(std::numeric_limits<mwIndex>::max()); 
    }
};

#endif /*  _SPARSELIST_H_  */