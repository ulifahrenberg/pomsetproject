# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 09:22:54 2018

@author: uli
"""

import poset
import itertools

class GluingInterfacesDontMatchError(Exception):
    def __init__(self, ipo_object1, ipo_object2):
        self.ipo_object1 = ipo_object1
        self.ipo_object2 = ipo_object2

class StartingInterfaceNotMinimalError(Exception):
    def __init__(self, ipo_object, vertex):
        self.ipo_object = ipo_object
        self.vertex = vertex
        #raise

class TerminatingInterfaceNotMaximalError(Exception):
    def __init__(self, ipo_object, vertex):
        self.ipo_object = ipo_object
        self.vertex = vertex

class InterfacesNotUnique(Exception):
    def __init__(self, ipo_object, intf):
        self.ipo_object = ipo_object
        self.intf = intf

class Iposet:
    """Data attributes:
    _poset -- underlying poset
    _s, _t -- tuples of vertices in the starting resp. terminating interfaces
    """
    def __init__(self, _poset, s, t):
        self._poset = _poset
        self._s = s
        self._t = t
        # are interfaces minimal / maximal?
        for x in s:
            #print(x, s, _poset)
            if not _poset.is_minimal(x):
                raise StartingInterfaceNotMinimalError(self, x) from Exception
        for x in t:
            if not _poset.is_maximal(x):
                raise TerminatingInterfaceNotMaximalError(self, x) from Exception
        # are interfaces unique?
        if len(s) != len(set(s)):
            raise InterfacesNotUnique(self, s) from Exception
        if len(t) != len(set(t)):
            raise InterfacesNotUnique(self, t) from Exception
        self._s = s
        self._t = t

    def __eq__(self, other):
        return self._poset == other._poset and self._s == other._s and self._t == other._t

    def __hash__(self):
        return hash(hash(self._poset) + hash(self._s) + hash(self._t))

    def __repr__(self):
        return "Iposet(Poset({}, {}), {}, {})".format(
                self._poset._vertices, self._poset.minrep(), self._s, self._t)

    @classmethod
    def from_str(cls, s):
        """Return Iposet generated from serial representation."""
        return eval(s)

    def to_poset(self):
        """Image of forgetful functor."""
        return self._poset

    @classmethod
    def from_poset(cls, _poset):
        """Image of cofree functor."""
        return cls(_poset, (), ())

    def is_poset(self):
        """Is in the image of the cofree functor."""
        return self._s == () and self._t == ()

    def is_connected(self, print_counterexample=False):
        """Return True iff self is (weakly) connected."""
        return self._poset.is_connected(print_counterexample)

    def to_tikz(self, vertex_placement=None):
        """Convert self to tikz code for drawing.
        vertex_placement -- matrix of vertices in xy style (NOT IMPLEMENTED)
        BAILS OUT if a vertex is both source and target
        """
        retl = []
        successors = self._poset.hasse_diagram()
        if vertex_placement:
            pass #TODO
        else:
            col = list(self._poset.minimals())
            xpos = 0
            while col != []:
                nextcol = set()
                for x in col:
                    nextcol |= successors[x]
                    #print(x, nextcol)
                #print(col, nextcol)
                for i in range(len(col)):
                    if col[i] not in nextcol: # otherwise, don't print it now!
                        if col[i] in self._s and col[i] in self._t:
                            #FIXME not implemented
                            raise ValueError('not implemented')
                            #retl.append('\\node ({}) at ({},{}) {{\\inoutpt}};\n'.format(col[i], xpos, -i))
                        elif col[i] in self._s:
                            retl.append('\\node [label=left:{{\\tiny {}}}] ({}) at ({},{}) {{\\inpt}};\n'.format(
                                    self._s.index(col[i]) + 1, col[i], xpos, -i))
                        elif col[i] in self._t:
                            retl.append('\\node [label=right:{{\\tiny {}}}] ({}) at ({},{}) {{\\outpt}};\n'.format(
                                    self._t.index(col[i]) + 1, col[i], xpos, -i))
                        else:
                            retl.append('\\node ({}) at ({},{}) {{\\intpt}};\n'.format(col[i], xpos, -i))
                col = list(nextcol)
                xpos += 1
        for x in self._poset._vertices:
            for y in successors[x]:
                retl.append('\\path ({}) edge ({});\n'.format(x, y))
        with open('test.txt', 'w') as f: f.write(''.join(retl))
        return ''.join(retl)

    def is_iso(self, other):
        """Return True iff self is isomorphic to other."""
        # obvious things first
        if len(self._poset._vertices) != len(other._poset._vertices):
            return False
        if len(self._poset._edges) != len(other._poset._edges):
            return False
        if len(self._s) != len(other._s):
            return False
        if len(self._t) != len(other._t):
            return False
        selfv = tuple(self._poset._vertices)
        for otherv in itertools.permutations(other._poset._vertices):
            permgood = True
            # check interfaces first
            for i in range(len(selfv)):
                if selfv[i] in self._s:
                    if otherv[i] not in other._s:
                        permgood = False
                        break
                    if self._s.index(selfv[i]) != other._s.index(otherv[i]):
                        permgood = False
                        break
                if selfv[i] in self._t:
                    if otherv[i] not in other._t:
                        permgood = False
                        break
                    if self._t.index(selfv[i]) != other._t.index(otherv[i]):
                        permgood = False
                        break
                if otherv[i] in other._s:
                    if selfv[i] not in self._s:
                        permgood = False
                        break
                    if self._s.index(selfv[i]) != other._s.index(otherv[i]):
                        permgood = False
                        break
                if otherv[i] in other._t:
                    if selfv[i] not in self._t:
                        permgood = False
                        break
                    if self._t.index(selfv[i]) != other._t.index(otherv[i]):
                        permgood = False
                        break
            # check edges
            if permgood:
                for i in range(len(selfv)):
                    for j in range(len(selfv)):
                        if i != j:
                            if (selfv[i], selfv[j]) in self._poset._edges and (otherv[i], otherv[j]) not in other._poset._edges:
                                permgood = False
                                break
                            if (selfv[i], selfv[j]) not in self._poset._edges and (otherv[i], otherv[j]) in other._poset._edges:
                                permgood = False
                                break
                    if permgood == False:
                        break
            if permgood: return True
        return False

    def parallel(self, other):
        """Return parallel product of self and other."""
        pposet, oldnew1, oldnew2 = self._poset.parallel(other._poset, _return_mappings=True)
        ps, pt = [], []
        for i in self._s:
            ps.append(oldnew1[i])
        for i in other._s:
            ps.append(oldnew2[i])
        for i in self._t:
            pt.append(oldnew1[i])
        for i in other._t:
            pt.append(oldnew2[i])
        return Iposet(pposet, tuple(ps), tuple(pt))

    def gluing(self, other):
        """Return gluing product of self below other."""
        if len(self._t) != len(other._s): #interfaces don't match!
            raise GluingInterfacesDontMatchError(self, other) from Exception
        # create vertices; keep track where they come from
        vertices, count, oldnew1, oldnew2 = set(), 0, dict(), dict()
        for v in self._poset._vertices:
            count += 1
            vertices.add(count)
            oldnew1[v] = count
        for v in other._poset._vertices:
            if v in other._s:
                pos = other._s.index(v)
                corrtarget = self._t[pos] # the self vertex with which v has been identified
                oldnew2[v] = oldnew1[corrtarget]
            else:
                count += 1
                vertices.add(count)
                oldnew2[v] = count
        #print(vertices, oldnew1, oldnew2)
        # create edges
        edges = set()
        for (x, y) in self._poset._edges:
            edges.add((oldnew1[x], oldnew1[y]))
        for (x, y) in other._poset._edges:
            edges.add((oldnew2[x], oldnew2[y]))
        for x in self._poset._vertices:
            if x not in self._t:
                for y in other._poset._vertices:
                    if y not in other._s:
                        edges.add((oldnew1[x], oldnew2[y]))
        #print(edges)
        ps, pt = [], []
        for i in self._s:
            ps.append(oldnew1[i])
        for i in other._t:
            pt.append(oldnew2[i])
        return Iposet(poset.Poset(vertices, edges), tuple(ps), tuple(pt))


def testiposet():
    t = poset.Poset({1,2,3,4}, {(1,2), (3,4)})
    ti = Iposet(t, (1,), (2,))
    print('ti:', ti)
    print(Iposet.from_poset(t))
    print(ti.to_poset())
    print('ti.parallel(ti):', ti.parallel(ti))
    print('ti.gluing(ti):', ti.gluing(ti))
    N1 = Iposet(poset.Poset({1,2,3}, {(2,3)}), (), (3,)).gluing(
                Iposet(poset.Poset({1,2}, set()), (2,), ()))
    print('N1:', N1)
    N2 = Iposet(poset.Poset({1,2}, set()), (), (1,)).gluing(
                Iposet(poset.Poset({1,2,3}, {(1,2)}), (1,), ()))
    print('N2:', N2)
    N = Iposet(poset.Poset({1,2,3,4}, {(1,2), (3,2), (3,4)}), (), ())
    print('N:', N)
    print('N.gluing(N):', N.gluing(N))
    Ni = Iposet(poset.Poset({1,2,3,4}, {(1,2), (3,2), (3,4)}), (), (2,))
    iN = Iposet(poset.Poset({1,2,3,4}, {(1,2), (3,2), (3,4)}), (1,), ())
    print('Ni.gluing(iN):', Ni.gluing(iN))
    Nij = Iposet(poset.Poset({1,2,3,4}, {(1,2), (3,2), (3,4)}), (), (2,4))
    ijN = Iposet(poset.Poset({1,2,3,4}, {(1,2), (3,2), (3,4)}), (1,3), ())
    print('Nij.gluing(ijN):', Nij.gluing(ijN))
    print(Nij)
    #return Nij.to_tikz()
    p = Iposet(poset.Poset({1,2,3,4}, {(1,3), (3,4), (2,4)}), (), ())
    print(p)
    p.to_tikz()

