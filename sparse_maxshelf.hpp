/** 
  * A soft max-bucket sort -- replaces a heap when
  * it is expensive to make lots of updates to the heap,
  * and you don't always need access to max.
  *
  * Each shelf holds values in a different interval;
  * this shelf-system has log-spaced intervalues:
  * values in shelf[j] are divided
  * roughly as follows:
  *     s^j <= shelf[j] < s^(j-1)
  *
  * where s is an input parameter between 0 and 1.
  */
  
#ifndef _SPARSE_MAXSHELF_H_
#define _SPARSE_MAXSHELF_H_

#include "mex.h"

#include <vector>
#include <stdlib.h>
#include <limits>
#include <math.h>
#include <utility> // for std::pair

#include "sparsehash/dense_hash_map.h"


#define DEBUGPRINT3(x) do { if (debugflag>=2) { \
mexPrintf x; mexEvalString("drawnow"); } \
} while (0)


class sparse_max_shelf {
    
    public:
    const double theta;  // `decay` factor of the shelves
    const double eps_min; // smallest value of eps
    const double eps_max; // starting value of eps
    const double log2_theta;

    mwIndex max_size;
    
    google::dense_hash_map<mwIndex, double> values;
    google::dense_hash_map<mwIndex, std::pair<mwIndex,mwIndex> > index2shelf;
    
    std::vector< std::vector<mwIndex> > shelves; // shelves[i][j] gives index in values of j^th entry on i^th shelf

    mwIndex num_shelves;
    mwIndex lastval;
    std::pair<mwIndex,mwIndex> empty_pair;
    
    int debugflag;
    
    sparse_max_shelf(double _theta, double _eps_min, double _eps_max, mwIndex _max_size, int _debugflag)
    : theta(_theta), eps_min(_eps_min), eps_max(_eps_max), log2_theta(log2(_theta)), max_size(_max_size),
      shelves( 1 + ceil(log2(_eps_min/_eps_max)/log2(_theta)) ), debugflag(_debugflag)
    {
        num_shelves = shelves.size();
    	lastval = std::numeric_limits<mwIndex>::max();
    	empty_pair = std::make_pair (lastval, lastval);
        values.set_empty_key(_max_size);
        index2shelf.set_empty_key(_max_size);
    }


    /**
    *   Given a node_ID, return what shelf it is on. Return lastval if not on a shelf
    */
    mwIndex get_shelf(mwIndex node_ID){
        return std::get<0>(index2shelf[node_ID]);
    }

   /**
    * Determines highest-up non-empty shelf, starting from input top_shelf.
    * Returns 0 if all empty, 1 otherwise.
    */
    bool find_top_shelf( mwIndex& top_shelf ) {  
        for ( ; top_shelf < num_shelves; top_shelf++){
            if ( shelves[top_shelf].empty() == 0 ){                 
            return 1; }
        }
        top_shelf = 0;
        return 0;
    }

   /**
    * Returns max entry on top_shelf, looks for top non-empty shelf starting from top_shelf
    */
    double look_max( mwIndex& top_shelf ) {
        if ( find_top_shelf(top_shelf) == 0 ){ return 0.0;}
        double max = 0.0;
        for (mwIndex j = 0; j < shelves[top_shelf].size(); j++){
            mwIndex ind = shelves[top_shelf][j];
            if (ind == lastval ) { continue;}
            if ( values[ind] > max ){ max = values[ind]; }
        }       
        return max;
    }

   /**
    * Returns index of the shelf 'value' belongs in.
    * Returning num_shelves means the value doesn't go inside the shelf.
    * Returns 0 for all entries `too big` to go into the shelf.
    */
    mwIndex which_shelf(double value) {
        if (value < eps_min) {
            return num_shelves;
        }
        return std::max(0.0,ceil(log2(value/eps_max)/log2_theta));
    }


    /**
     * Put new entry on a shelf determined by 'which_shelf', return that shelf_ID.
     * The only way this is called is if update() determines an entry is not
     * in the shelf.
     */
    mwIndex insert(mwIndex node_ID, double value){
        mwIndex shelf_ID = which_shelf(value);
        if (shelf_ID >= num_shelves){
            return num_shelves;
        }       
        shelves[shelf_ID].push_back(node_ID);
        mwIndex ind2 = shelves[shelf_ID].size()-1;
        shelves[shelf_ID][ind2] = node_ID;
        index2shelf[node_ID] = std::make_pair(shelf_ID,ind2);
        return shelf_ID;
    }

    /**
     * put entry on appropriate shelf, returns ID of shelf where it goes
     */
    mwIndex update(mwIndex node_ID, double value){
        // update val, then check if it's in shelf
        double newval = values[node_ID] + value;
        values[node_ID] = newval;

        // if not in shelf, insert() it.
        if ( index2shelf.count(node_ID) == 0 ) { 
            return insert(node_ID, newval); 
        }
        if ( index2shelf[node_ID] == empty_pair ) {
            return insert(node_ID, newval);
        }
        
        // ... then determine if shelf location needs to change
        mwIndex old_shelf_ID = std::get<0>(index2shelf[node_ID]);
        mwIndex old_shelf_spot = std::get<1>(index2shelf[node_ID]);

        mwIndex shelf_ID = which_shelf(fabs(newval));
        if (shelf_ID == old_shelf_ID){ return shelf_ID; } // already where it belongs!

        // delete the obsolete entry by swapping it to back and pop_back() ing
        mwIndex ind_end = swap_to_end(old_shelf_spot, old_shelf_ID);
        shelves[old_shelf_ID][ind_end] = lastval; // clear from shelf
        shelves[old_shelf_ID].pop_back();

        if (shelf_ID >= num_shelves){ // too small for shelf, make sure it's removed
            index2shelf[node_ID] = empty_pair;        
            return num_shelves;
        }
    
        // then put it on its new shelf
        shelves[shelf_ID].push_back(node_ID);
        mwIndex ind2 = shelves[shelf_ID].size()-1;
        index2shelf[node_ID] = std::make_pair(shelf_ID,ind2);
        return shelf_ID;
    }
    
    
    /**
    *   Swaps entry on shelf `shelf_ID` with shelf_spot ind1 to end of the shelf
    */
    mwIndex swap_to_end(mwIndex ind1, mwIndex shelf_ID){
        int ind_end = shelves[shelf_ID].size()-1;
        assert( ind_end >= 0 );
        if (ind1 == (mwIndex)ind_end){return ind1;} // already at back
        mwIndex temp_ind = shelves[shelf_ID][ind1];
        mwIndex temp_ind_end = shelves[shelf_ID][ind_end];
        shelves[shelf_ID][ind1] = temp_ind_end;
        shelves[shelf_ID][ind_end] = temp_ind;
        index2shelf[temp_ind] = std::make_pair(shelf_ID,(mwIndex)ind_end);
        index2shelf[temp_ind_end] = std::make_pair(shelf_ID,ind1);
        return (mwIndex)ind_end;
    }

    
    /**
     * pops back on specified shelf, sets node_ID to the index of that entry,
     * sets value to values[node_ID], returns 1 if possible, 0 if empty
     */
    bool pop_from_shelf( mwIndex current_shelf, mwIndex& node_ID, double& value ){
        mwIndex shelf_size = shelves[current_shelf].size();
        if (shelf_size == 0){ return 0;} // if shelf is empty, do nothing
        node_ID = shelves[current_shelf].back();

        // clear entry in shelf
        shelves[current_shelf][shelf_size-1] = lastval;
        shelves[current_shelf].pop_back();
        if (node_ID == lastval){ return 0; } // 0 signals the entry was a blank
        index2shelf[node_ID] = empty_pair; // entry erased from shelf and index2shelf

        // pass residual out and clear inside shelf
        value = values[node_ID];
        values[node_ID] = 0.0;
        return 1;
    }
};

#endif  /*  _SPARSE_MAXSHELF_H_  */