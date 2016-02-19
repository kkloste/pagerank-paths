#!/usr/bin/env python

"""
This Python script has a set of pseudo-code for the pprpath or pprall 
algorithms that Kyle Kloster and I are developing.

We begin with a highly simplified implementation and proceed to more 
intricate implementations. We have three test cases:

netscience-cc
pgp-cc
itdk-cc

These will be test cases for the staleshelf structure
"""

import sys
import math
import time
import collections

from heapq import heapify, heappush, heappop

from staleshelf import StaleShelf, BucketShelf

class SmatGraph(dict):
    def __init__(self, filename):
        dict.__init__(self)
        with open(filename) as file:
            header = file.readline().split()
            self.nverts = int(header[0])
            nlines = int(header[2])
            nedges = 0
            nself = 0
            for line in file:
                parts = line.split()
                src = int(parts[0])
                dst = int(parts[1])
                nedges += 1
                if src == dst:
                    nself += 1
                    continue # don't add these
                
                if src in self:
                    self[src].add(dst)
                else:
                    self[src] = set([dst])
            
        if nlines != nedges:
            print >> sys.stderr, "Expecting %i lines, got %i edges"%(nlines, nedges)
            assert(nlines == nedges)
                
        self.nedges = nedges
        
        # check symmetric
        for u,neigh in self.iteritems():
            for v in neigh:
                if u not in self[v]:
                    print >>sys.stderr, "Expecting an undirected graph" 
                    print >>sys.stderr, "Found (%i,%i) but not (%i,%i)"%(
                        u,v, v,u)
                    assert( False )
        
        
    def out_degree(self,v):
        return len(self[v])
        
class priority_dict(dict):
    """Dictionary that can be used as a priority queue.

    Keys of the dictionary are items to be put into the queue, and values
    are their respective priorities. All dictionary methods work as expected.
    The advantage over a standard heapq-based priority queue is
    that priorities of items can be efficiently updated (amortized O(1))
    using code as 'thedict[item] = new_priority.'

    The 'smallest' method can be used to return the object with lowest
    priority, and 'pop_smallest' also removes it.

    The 'sorted_iter' method provides a destructive sorted iterator.
    """
    
    def __init__(self, *args, **kwargs):
        super(priority_dict, self).__init__(*args, **kwargs)
        self._rebuild_heap()

    def _rebuild_heap(self):
        self._heap = [(v, k) for k, v in self.iteritems()]
        heapify(self._heap)

    def smallest(self):
        """Return the item with the lowest priority.

        Raises IndexError if the object is empty.
        """
        
        heap = self._heap
        v, k = heap[0]
        while k not in self or self[k] != v:
            heappop(heap)
            v, k = heap[0]
        return k

    def pop_smallest(self):
        """Return the item with the lowest priority and remove it.

        Raises IndexError if the object is empty.
        """
        
        heap = self._heap
        v, k = heappop(heap)
        while k not in self or self[k] != v:
            v, k = heappop(heap)
        del self[k]
        return k

    def __setitem__(self, key, val):
        # We are not going to remove the previous value from the heap,
        # since this would have a cost O(n).
        
        super(priority_dict, self).__setitem__(key, val)
        
        if len(self._heap) < 2 * len(self):
            heappush(self._heap, (val, key))
        else:
            # When the heap grows larger than 2 * len(self), we rebuild it
            # from scratch to avoid wasting too much memory.
            self._rebuild_heap()

    def setdefault(self, key, val):
        if key not in self:
            self[key] = val
            return val
        return self[key]

    def update(self, *args, **kwargs):
        # Reimplementing dict.update is tricky -- see e.g.
        # http://mail.python.org/pipermail/python-ideas/2007-May/000744.html
        # We just rebuild the heap from scratch after passing to super.
        
        super(priority_dict, self).update(*args, **kwargs)
        self._rebuild_heap()

    def sorted_iter(self):
        """Sorted iterator of the priority dictionary items.

        Beware: this will destroy elements as they are returned.
        """
        
        while self:
            yield self.pop_smallest()        
        
def ppr_naive(G, s, epsmin):
    """ Run a naive PPR to a tolerance of epsilon. """
    
    alpha = 0.99
    epscur = 1.
    epsm = (1.-sys.float_info.epsilon)
    x = [0. for i in xrange(G.nverts)]
    r = [0. for i in xrange(G.nverts)]
    rpri = priority_dict()
    npush = 0
    pushvol = 0
    
    r[s] = 1./G.out_degree(s)
    rpri[s] = -1./G.out_degree(s)
    npusheps = -1
    
    while len(rpri) > 0:
        v = rpri.pop_smallest()
        curr = r[v] 
        
        if curr < epscur:
            epscur = curr*epsm
            #print "%.10f  %10i  %10i"%(epscur, npush, pushvol)

        #print "push on %i", v

        pushamount = alpha*curr
        x[v] += curr
        r[v] = 0.
        
        for u in G[v]:
            rchange = pushamount/G.out_degree(u)
            rucur = r[u]
            runew = rucur + rchange
            r[u] = runew
            if runew > epsmin:
                rpri[u] = -runew
            pushvol += 1
        npush += 1
        
        
def ppr_shelf(G,s,epsmin):
    alpha = 0.99
    onesmall = (1.-sys.float_info.epsilon)
    x = [0. for i in xrange(G.nverts)]
    r = [0. for i in xrange(G.nverts)]
    npush = 0
    pushvol = 0
    
    rho = 0.
    rhom1 = 1.-rho
    
    r[s] = 1./G.out_degree(s)
    shelf = StaleShelf(epsmin,G.nverts,r)
    shelf.update(s)
    
    epscur = r[s]*onesmall
    
    while epscur >= epsmin:
        Q = collections.deque(shelf.iter_above(epscur))
        changedset = []
        changedsetind = 0
        
        #print "%.10f  %10i  %10i"%(epscur, npush, pushvol)
        
        while len(Q) > 0:
            v = Q.popleft()
            curr = r[v] 
            
            #print "push on %i"%(v)
            
            pushamount = curr - epscur*rho
            x[v] += pushamount
            pushamount *= alpha # fix the push amount
            r[v] = epscur*rho
            
            if r[v] > epsmin:
                shelf.update(v)
            else:
                shelf.remove(v)
            
            for u in G[v]:
                rchange = pushamount/G.out_degree(u)
                rucur = r[u]
                runew = rucur + rchange
                r[u] = runew
                
                if rucur >= epscur:
                    pass # must be in the queue already
                elif rucur < epscur and runew >= epscur:
                    Q.append(u)
                elif runew >= epsmin:
                    shelf.update(u)
                    
                pushvol += 1
            npush += 1
            
        epscur = onesmall*shelf.find_max()[0]
            
def ppr_grid(G,s,epsmin):
    alpha = 0.99
    theta = 0.95
    x = [0. for i in xrange(G.nverts)]
    r = [0. for i in xrange(G.nverts)]
    npush = 0
    pushvol = 0
    
    rho = 0.
    rhom1 = 1.-rho
    
    r[s] = 1./G.out_degree(s)
    shelf = BucketShelf(theta,epsmin,r)
    shelf.update(s)
    
    epsbound = r[s]
    epscur = 1.
    
    while epsbound >= epsmin:
        # find the next value of eps to solve
        
        while epscur > epsbound:
            epscur *= theta
            
        Q = collections.deque(shelf.iter_above_and_clear(epscur))
        changedset = []
        changedsetind = 0
        
        while len(Q) > 0:
            v = Q.popleft()
            curr = r[v] 
            
            #print "push on %i"%(v)
            
            pushamount = curr - epscur*rho
            x[v] += pushamount
            pushamount *= alpha # fix the push amount
            r[v] = epscur*rho
            
            if r[v] > epsmin:
                shelf.update(v)
            
            for u in G[v]:
                rchange = pushamount/G.out_degree(u)
                rucur = r[u]
                runew = rucur + rchange
                r[u] = runew
                
                if rucur >= epscur:
                    pass # must be in the queue already
                elif rucur < epscur and runew >= epscur:
                    Q.append(u)
                elif runew >= epsmin:
                    shelf.update(u)
                    
                pushvol += 1
            npush += 1
            
        epsbound = shelf.find_max()[0]
        #print "%.10f  %10i  %10i"%(epsbound, npush, pushvol)
        

def perfrun():
    G = SmatGraph('itdk0304-cc.smat')
    for epsmin in [ 1.e-4, 1.e-5, 1.e-7 ]:
        v = 1
        t0 = time.time() 
        ppr_grid(G, v, epsmin)
        t1 = time.time()
        print "%8s  %5.1e  %6.3f sec"%("grid", epsmin, t1-t0)
        
        t0 = time.time() 
        ppr_naive(G, 1, epsmin)
        t1 = time.time()
        print "%8s  %5.1e  %6.3f sec"%("naive", epsmin, t1-t0)

        t0 = time.time() 
        ppr_shelf(G, 1, epsmin)
        t1 = time.time()
        print "%8s  %5.1e  %6.3f sec"%("shelf", epsmin, t1-t0)
        
if __name__ == '__main__':
    if sys.argv[1] == 'perf':
        perfrun()
        sys.exit(0)
        
    G = SmatGraph('netscience-cc.smat')       
    epsmin = 0.00326814570000 
    ppr_grid(G, 1, 1.e-7)
    #ppr_grid(G, 1, 1.e-7)
    #ppr_naive(G, 1, epsmin)
    #ppr_shelf(G, 1, epsmin)
    
""" New idea for resweep algorithm.
We are looking to stop resweeping unless there is a large gap in the tail 
of the vector.
There is always a gap of at least epsmin. A gap beyond     
