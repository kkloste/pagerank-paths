import sys
import collections
import math

class BucketShelf:
    def __init__(self,theta,minval,vals,location=None):
        nshelves = int(math.ceil(math.log(minval)/math.log(theta)))
        self._shelves = [ list() for i in xrange(nshelves) ]
        if location is None:
            self._location = {}
        else:
            self._location = location
        self._vals = vals
        self._minval = minval
        self._theta = theta
        
        
    def update(self,item):
        """ Update the position in the shelf for item. 
        Note that update will ignore any items with values less than minval.
        So you shouldn't use this to decrease items. Remove those first.
        """
        
        val = self._vals[item]
        if val < self._minval:
            return
            
        newshelfid = int(math.floor(math.log(val)/math.log(self._theta)))
        if item not in self._location or self._location[item][0] is None:
            # then the item isn't in the shelf yet
            self._shelves[newshelfid].append(item)
            self._location[item] = (len(self._shelves[newshelfid]) - 1,newshelfid)
        else:
            # then the item is in the shelf
            (pos, oldshelf) = self._location[item]
            if oldshelf is not newshelfid:
                # then we need to move it.
                self._remove(item, oldshelf, pos)
                self._shelves[newshelfid].append(item)
                self._location[item] = (len(self._shelves[newshelfid]) - 1,newshelfid)
                
    def remove(self, item):
        #assert( self._vals[item] < self._minval )
        (pos,shelfid) = self._location[item]
        if pos is not None:
            self._remove(item, shelfid, pos)
            
    def _remove(self, item, shelf, pos):
        """ Remove an item that we know the location """
        s = self._shelves[shelf] # get a ref to the shelf
        assert(len(s) > 0)
        if len(s) > 1:
            # then there will be something left, so we are going to 
            # pivot the element to the end, and then pop it later   
            endpos = len(s)-1
            lastitem = s[endpos]
            
            s[pos] = lastitem
            s[endpos] = item    
            self._location[lastitem] = (pos, shelf)
        else:
            assert( pos == 0 )
        
        # last step, remove the item
        s.pop() # remove the last thing
        self._location[item] = (None, 0) # and update 
        
    def find_max(self):
        """ 
        Return, (-1, None) if the shelf is empty
        otherwise, return (max, item) if the shelf is not empty
        
        This will also fix up the first row of the shelf if it isn't valid.
        """
        smin = 1.
        curmax = -1 # this is guaranteed to be smaller than everything
        maxitem = None
        validshelf = False
        for s in self._shelves:
            smin = smin*self._theta
            for item in s:
                val = self._vals[item]
                if val >= smin:
                    validshelf = True
                if val > curmax:
                    maxitem = item
                    curmax = val
            if validshelf:
                break # we an stop looking at later shelves
                
        return (curmax, maxitem)
        
    def iter_above(self, limit):
        """ Return a list of all items above a threshold. """
        smin = 1.
        for s in self._shelves:
            smin = smin*self._theta
            for item in s:
                val = self._vals[item]
                if val >= limit:
                    yield item
            if smin < limit:
                break
                
    def iter_above_and_clear(self, limit):
        smin = 1.
        itemlist = list(self.iter_above(limit))
        for i in itemlist:
            self.remove(i)
        return itemlist
                

class StaleShelf:
    """
    We assume that there is a finite universe and so we can use
    a lookup table instead of a hash table for true O(1) lookups
    
    We assume that values are in [0,1).
    
    We assume that values can only increase. If they can decrease, then 
    you must manually remove them, and re-update them.
    
    """
    
    def __init__(self,minval,nuniverse,valdict):
        """ 
        minval - the smallest value that should be tracked
        nuniverse -
        valdict - the values that determine behavior in the shelf
        """
        self._nuniv = nuniverse
        self._location = [ (nuniverse,0) for i in xrange(nuniverse) ]
        nshelves = -math.frexp(minval)[1]+1
        self._shelves = [ list() for i in xrange(nshelves) ]
        self._shelves_sorted = [ False for i in xrange(nshelves) ]
        self._vals = valdict
        self._minval = minval
        
    def update(self,item):
        """ Update the position in the shelf for item. 
        Note that update will ignore any items with values less than minval.
        So you shouldn't use this to decrease items. Remove those first.
        """
        
        val = self._vals[item]
        if val < self._minval:
            return
            
        newshelfid = -math.frexp(val)[1]
        if self._location[item][0] == self._nuniv:
            # then the item isn't in the shelf yet
            self._shelves[newshelfid].append(item)
            self._location[item] = (len(self._shelves[newshelfid]) - 1,newshelfid)
            self._shelves_sorted[newshelfid] = len(self._shelves[newshelfid]) > 1
        else:
            # then the item is in the shelf
            (pos, oldshelf) = self._location[item]
            #print (pos, oldshelf)
            if oldshelf is not newshelfid:
                # then we need to move it.
                self._remove(item, oldshelf, pos)
                self._shelves[newshelfid].append(item)
                self._shelves_sorted[newshelfid] = len(self._shelves[newshelfid]) > 1
                self._location[item] = (len(self._shelves[newshelfid]) - 1,newshelfid)
                
    def remove(self, item):
        #assert( self._vals[item] < self._minval )
        (pos,shelfid) = self._location[item]
        if pos is not self._nuniv:
            self._remove(item, shelfid, pos)
            
    def _sort_shelf(self, shelf):
        return sorted( shelf, reverse=True, key=lambda item: self._vals[item] )

    def _remove(self, item, shelf, pos):
        """ Remove an item that we know the location """
        s = self._shelves[shelf] # get a ref to the shelf
        assert(len(s) > 0)
        if len(s) > 1:
            # then there will be something left, so we are going to 
            # pivot the element to the end, and then pop it later   
            endpos = len(s)-1
            lastitem = s[endpos]
            
            s[pos] = lastitem
            s[endpos] = item    
            self._location[lastitem] = (pos, shelf)
        else:
            assert( pos == 0 )
        
        # last step, remove the item
        s.pop() # remove the last thing
        self._location[item] = (self._nuniv, 0) # and update 
        
    def find_max(self):
        """ 
        Return, (-1, None) if the shelf is empty
        otherwise, return (max, item) if the shelf is not empty
        
        This will also fix up the first row of the shelf if it isn't valid.
        """
        smin = 1.
        curmax = -1 # this is guaranteed to be smaller than everything
        maxitem = None
        validshelf = False
        for s in self._shelves:
            smin = smin/2.
            for item in s:
                val = self._vals[item]
                if val >= smin:
                    validshelf = True
                if val > curmax:
                    maxitem = item
                    curmax = val
            if validshelf:
                break # we an stop looking at later shelves
                
        return (curmax, maxitem)
        
    def iter_above(self, limit):
        """ Return a list of all items above a threshold. """
        smin = 1.
        for s in self._shelves:
            smin = smin/2.
            for item in s:
                val = self._vals[item]
                if val >= limit:
                    yield item
            if smin < limit:
                break
                
    
def test1():
    
    x = [3.,5.,1.,2.,4.]
    x = [i/(max(x) + 1.) for i in x]
    print x
    nuniv = len(x)
    shelf = StaleShelf(2.5/6., nuniv, x)
    for i in xrange(nuniv):
        shelf.update(i)
    print shelf.find_max()
    assert( shelf.find_max() == (5./(6.), 1 ) )
    x[1] = 2.6/6.
    shelf.update(1)
    print shelf.find_max()
    assert( shelf.find_max() == (4./(6.), 4 ) )
    
    x[1] = 0.
    print shelf.find_max()
    assert( shelf.find_max() == (4./(6.), 4 ) )
    
    x[1] = 0.9
    shelf.update(1)
    print shelf.find_max()
    assert( shelf.find_max() == (0.9, 1 ) )
    
    
    shelf.remove(1)
    print shelf.find_max()
    assert( shelf.find_max() == (4./6., 4) )
    
    shelf.remove(4)
    print shelf.find_max()    
    assert( shelf.find_max() == (3./6., 0) )
    
    shelf.remove(0)
    print shelf.find_max()    
    assert( shelf.find_max() == (-1, None) )
    
if __name__=='__main__':
    # run our tests
    test1()    
        
        