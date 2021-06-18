# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 09:22:54 2018

@author: uli
"""

import itertools
import more_itertools as mitertools

class PartialOrderNotAntisymmetricError(Exception):
    def __init__(self, po_object):
        self.po_object = po_object

class Poset:
    """Data attributes:
    vertices -- finite set
    edges -- finite set of tuples of vertices which is guaranteed to be
    reflexive, transitive and antisymmetric
    """
    def __init__(self, vertices, edges):
        self._vertices = vertices
        # transitive closure of edges
        done = False
        while not done:
            done = True
            e2 = edges.copy()
            for (x, y1) in e2:
                for (y2, z) in e2:
                    if y1 == y2:
                        if (x, z) not in edges:
                            edges.add((x, z))
                            done = False
        self._edges = set(edges)
        # check for antisymmetry
        for x in self._vertices:
            for y in self._vertices:
                if x != y:
                    if (x, y) in edges and (y, x) in edges:
                        raise PartialOrderNotAntisymmetricError(self) from Exception
        # reflexive closure
        edges |= set((x, x) for x in self._vertices)
        self._edges = edges

    def __eq__(self, other):
        if self._vertices != other._vertices:
            return False
        return self._edges == other._edges

    def __hash__(self):
        return hash(frozenset(self._edges))

    def __repr__(self):
        return "Poset({}, {})".format(self._vertices, self.minrep())

    def add_edge(self, e):
        # does not check for asymmetry!
        self._edges.add(e)

    def add_vertex(self, v):
        self._vertices.add(v)

    def remove_vertex(self, v):
        self._vertices.remove(v)
        for e in self._edges.copy():
            if e[0] == v or e[1] == v:
                self._edges.remove(e)

    def deepcopy(self):
        return Poset(self._vertices.copy(), self._edges.copy())

    @classmethod
    def from_str(cls, s):
        """Return Poset generated from serial representation."""
        return eval(s)

    def hasse_diagram(self):
        """Return dict from vertices to sets of immediate successors."""
        return {x: self.successors(x) for x in self._vertices}

    def filtration(self):
        """Return poset filtration of self as list of nonempty disjoint sets of vertices."""
        pc = self.deepcopy()
        retl = []
        mi = pc.minimals()
        while mi:
            retl.append(mi)
            for v in mi:
                pc.remove_vertex(v)
            mi = pc.minimals()
        return retl

    def to_tikz(self, vertex_placement=None, vertex_labels=False):
        """Convert self to tikz code for drawing.
        vertex_placement -- matrix of vertices in xy style (NOT IMPLEMENTED)
        """
        retl = []
        successors = self.hasse_diagram()
        if vertex_placement:
            pass #TODO
        else:
            #TODO should use self.filtration()
            col = list(self.minimals())
            xpos = 0
            while col != []:
                for i in range(len(col)):
                    if vertex_labels:
                        retl.append('\\node ({}) at ({},{}) {{${}$}};\n'.format(col[i], xpos, -i, col[i]))
                    else:
                        retl.append('\\node ({}) at ({},{}) {{\\intpt}};\n'.format(col[i], xpos, -i))
                nextcol = set()
                for x in col:
                    nextcol |= successors[x]
                col = list(nextcol)
                xpos += 1
        for x in self._vertices:
            for y in successors[x]:
                retl.append('\\path ({}) edge ({});\n'.format(x, y))
        #with open('test.txt', 'w') as f: f.write(''.join(retl))
        return ''.join(retl)

    def to_tikz_olavi(self, rem_point, forbidden, vertex_labels=False, forb_right=True):
        """Convert self to tikz code for drawing.
        - rem_point: point which, when removed, gives forbidden
        - forbidden: name of forbidden substructure
        - vertex_labels: print numbers; else bullets
        - forb_right: print forbidden right of fig; else below
        """
        retl = []
        successors = self.hasse_diagram()
        filt = self.filtration()
        xrange = len(filt)
        yrange = 0 # greatest y coordinate
        for xpos in range(len(filt)):
            col = list(filt[xpos])
            yrange = min(yrange, 1-len(col))
            for i in range(len(col)):
                ypos = -i
                if xpos > 1:
                    ypos -= .3
                if col[i] == rem_point:
                    if vertex_labels:
                        retl.append('\\node[red] ({}) at ({},{}) {{${}$}};\n'.format(
                                col[i], xpos, ypos, col[i]))
                    else:
                        retl.append('\\node[red] ({}) at ({},{}) {{$\\bullet$}};\n'.format(
                                col[i], xpos, ypos))
                else:
                    if vertex_labels:
                        retl.append('\\node ({}) at ({},{}) {{${}$}};\n'.format(
                                col[i], xpos, ypos, col[i]))
                    else:
                        retl.append('\\node ({}) at ({},{}) {{\\intpt}};\n'.format(
                                col[i], xpos, ypos))
        for x in self._vertices:
            for y in successors[x]:
                if x == rem_point or y == rem_point:
                    retl.append('\\path ({}) edge[red] ({});\n'.format(x, y))
                else:
                    retl.append('\\path ({}) edge ({});\n'.format(x, y))
        if forb_right:
            retl.append('\\node[align=left] at ({},{}) {{{}}};\n'.format(
                    xrange+1, yrange/2, '-- '+forbidden))
        else:
            retl.append('\\node at ({},{}) {{{}}};\n'.format(
                    xrange/2, yrange-1, '-- '+forbidden))
        #with open('test.txt', 'w') as f: f.write(''.join(retl))
        return ''.join(retl)

    def precedes(self, x, y):
        """Return True iff x precedes y."""
        return (x, y) in self._edges

    def concurrent(self, x, y):
        """Return True iff x and y are concurrent / independent."""
        return (not self.precedes(x, y)) and (not self.precedes(y, x))

    def successors(self, x):
        """Return set of direct successors of node x."""
        rets = set()
        for y in self._vertices:
            if y != x and (x, y) in self._edges:
                found_one = True
                for z in self._vertices:
                    if z != x and z != y and (x, z) in self._edges and (z, y) in self._edges:
                        found_one = False
                        break
                if found_one:
                    rets.add(y)
        return rets

    def is_maximal(self, x):
        """Return True iff x is maximal."""
        for y in self._vertices:
            if x != y and (x, y) in self._edges:
                return False
        return True

    def is_minimal(self, y):
        """Return True iff x is maximal."""
        for x in self._vertices:
            if x != y and (x, y) in self._edges:
                return False
        return True

    def minimals(self):
        """Return set of minimal elements."""
        return set(x for x in self._vertices if self.is_minimal(x))

    def maximals(self):
        """Return set of maximal elements."""
        return set(x for x in self._vertices if self.is_maximal(x))

    def is_serial(self, witness=False):
        """Return True iff self is a serial product, plus witness."""
        if not self.is_connected():
            return False, None if witness else False
        p1, p2 = set(), set()
        # easy nodes first: minimals to p1; maximals to p2
        for x in self._vertices:
            if self.is_minimal(x):
                p1.add(x)
            elif self.is_maximal(x):
                p2.add(x)
        rest = self._vertices - p1 - p2
        print(p1, p2, rest)
        for parti in mypartitions2(rest):
            #print(p1, p2, parti)
            # parti[0] goes left
            p1c, p2c = p1.copy(), p2.copy()
            for x in parti[0]:
                p1c.add(x)
            for x in parti[1]:
                p2c.add(x)
            print(p1c, p2c)
            e1, e2 = set(), set()
            for x, y in self._edges:
                if x in p1c and y in p1c:
                    e1.add((x, y))
                elif x in p2c and y in p2c:
                    e2.add((x, y))
            P1 = Poset(p1c, e1)
            P2 = Poset(p2c, e2)
            print(P1, P2, P1.serial(P2))
            if self == P1.serial(P2):
                return True, (P1, P2) if witness else True
            #parti[0] goes right
            p1c, p2c = p1.copy(), p2.copy()
            for x in parti[0]:
                p2c.add(x)
            for x in parti[1]:
                p1c.add(x)
            e1, e2 = set(), set()
            for x, y in self._edges:
                if x in p1c and y in p1c:
                    e1.add((x, y))
                elif x in p2c and y in p2c:
                    e2.add((x, y))
            P1 = Poset(p1c, e1)
            P2 = Poset(p2c, e2)
            print(P1, P2)
            if self == P1.serial(P2):
                return True, (P1, P2) if witness else True
        return False, None if witness else False

    def is_connected(self, print_counterexample=False):
        """Return True iff self is (weakly) connected."""
        if len(self._vertices) <= 1: # empty and 1-point posets are connected
            return True
        rest = self._vertices.copy()
        first = rest.pop()
        coco = {first} # connected component of self._vertices[0]
        done = False
        while not done:
            done = True
            for x in coco.copy():
                for y in rest.copy():
                    if (x, y) in self._edges or (y, x) in self._edges:
                        coco.add(y)
                        rest.remove(y)
                        done = False
        if rest == set():
            return True
        else:
            if print_counterexample:
                r = rest.pop()
                print("Counterexample:", first, r)
            return False

    def is_stepseq(self, print_counterexample=False):
        """Return True iff self is a step sequence."""
        # check if self.concurrent is transitive
        for x in self._vertices:
            for y in self._vertices:
                if self.concurrent(x, y):
                    for z in self._vertices:
                        if z != x:
                            if self.concurrent(y, z):
                                if not self.concurrent(x, z):
                                    if print_counterexample:
                                        print("Counterexample:", x, y, z)
                                    return False
        return True

    def is_intervalorder(self, print_counterexample=False):
        """Return True iff self is an interval order."""
        # check if self._edges contains 2+2
        for x in self._vertices:
            for y in self._vertices:
                if x != y and self.precedes(x, y):
                    for z in self._vertices:
                        for w in self._vertices:
                            if z != w and self.precedes(z, w):
                                if (not self.precedes(x, w)) and (not self.precedes(z, y)):
                                    if print_counterexample:
                                        print("Counterexample:", x, y, z, w)
                                    return False
        return True

    def order_extensions(self):
        """Return all order extensions of self."""
        retval = {self}
        for x in self._vertices:
            for y in self._vertices:
                if not (x, y) in self._edges:
                    try:
                        selfext = Poset(self._vertices, self._edges | {(x, y)})
                        retval |= selfext.order_extensions()
                    except PartialOrderNotAntisymmetricError:
                        pass
        return retval

    def step_extensions(self):
        """Return all step extensions of self."""
        if self.is_stepseq():
            retval = {self}
        else:
            retval = set()
        for x in self._vertices:
            for y in self._vertices:
                if not (x, y) in self._edges:
                    try:
                        selfext = Poset(self._vertices, self._edges | {(x, y)})
                        retval |= selfext.step_extensions()
                    except PartialOrderNotAntisymmetricError:
                        pass
        return retval

    def interval_extensions(self):
        """Return all interval-order extensions of self."""
        if self.is_intervalorder():
            retval = {self}
        else:
            retval = set()
        for x in self._vertices:
            for y in self._vertices:
                if not (x, y) in self._edges:
                    #print(x, y)
                    try:
                        selfext = Poset(self._vertices, self._edges | {(x, y)})
                        retval |= selfext.interval_extensions()
                    except PartialOrderNotAntisymmetricError:
                        pass
        return retval

    def minrep(self):
        """Return minimal representation (Hasse diagram) of edges."""
        edges = self._edges.copy()
        # remove reflexive edges
        for (x, y) in self._edges:
            if x == y:
                edges.remove((x, y))
        # transitive reduction
        done = False
        while not done:
            done = True
            e2 = edges.copy()
            for (x, y) in e2:
                for (a, b) in e2:
                    if a == x:
                        for (c, d) in e2:
                            if b == c and d == y:
                                done = False
                                edges.discard((x, y))
        return edges

    def parallel(self, other, *, _return_mappings=False):
        """Return parallel product of self and other.
        if _return_mapping, then als return mappings from old to new vertices.
        """
        # create vertices; keep track where they come from
        vertices, count, oldnew1, oldnew2 = set(), 0, dict(), dict()
        for v in self._vertices:
            count += 1
            vertices.add(count)
            oldnew1[v] = count
        for v in other._vertices:
            count += 1
            vertices.add(count)
            oldnew2[v] = count
        # create edges
        edges = set()
        for (x, y) in self._edges:
            edges.add((oldnew1[x], oldnew1[y]))
        for (x, y) in other._edges:
            edges.add((oldnew2[x], oldnew2[y]))
        if _return_mappings:
            return Poset(vertices, edges), oldnew1, oldnew2
        else:
            return Poset(vertices, edges)

    def serial(self, other): #TODO make new version which doesn't rename vertices
        """Return serial product of self before other."""
        # create vertices; keep track where they come from
        vertices, count, oldnew1, oldnew2 = set(), 0, dict(), dict()
        for v in self._vertices:
            count += 1
            vertices.add(count)
            oldnew1[v] = count
        for v in other._vertices:
            count += 1
            vertices.add(count)
            oldnew2[v] = count
        # create edges
        edges = set()
        for (x, y) in self._edges:
            edges.add((oldnew1[x], oldnew1[y]))
        for (x, y) in other._edges:
            edges.add((oldnew2[x], oldnew2[y]))
        for x in self._vertices:
            for y in other._vertices:
                edges.add((oldnew1[x], oldnew2[y]))
        return Poset(vertices, edges)

    def is_iso(self, other):
        """Return True iff self is isomorphic to other."""
        # obvious things first
        if len(self._vertices) != len(other._vertices):
            return False
        if len(self._edges) != len(other._edges):
            return False
        selfv = tuple(self._vertices)
        for otherv in itertools.permutations(other._vertices):
            permgood = True
            # check edges
            for i in range(len(selfv)):
                for j in range(len(selfv)):
                    if i != j:
                        if (selfv[i], selfv[j]) in self._edges and (otherv[i], otherv[j]) not in other._edges:
                            permgood = False
                            break
                        if (selfv[i], selfv[j]) not in self._edges and (otherv[i], otherv[j]) in other._edges:
                            permgood = False
                            break
                if not permgood:
                    break
            if permgood: return True
        return False

def mypartitions2(s):
    """yield all partitions of s into two parts, but including empty parts."""
    if len(s) == 0:
        yield [[], []]
    else:
        yield [list(s), []]
        yield [[], list(s)]
        yield from mitertools.set_partitions(s,2)

def write_to_file(posets, fname):
    with open(fname, 'w') as f:
        for p in posets:
            f.write(str(p) + '\n')

def read_posets_from_file(fname):
    ps = set()
    with open(fname, 'r') as f:
        for l in f.readlines():
            ps.add(Poset.from_str(l))
    return ps

def read_posets_from_file_olavi(fname):
    ps = []
    with open(fname) as f:
        for l in f:
            if l[:5] == 'Graph':
                p = Poset({1,2,3,4,5,6,7},set())
                rem_point, forbidden = None, None
            elif l[:3] == 'rem':
                parts = l.strip().split()
                rem_point = int(parts[1])
                forbidden = parts[3]
            elif l == 'end\n':
                ps.append((p, rem_point, forbidden))
            elif l == '\n':
                continue
            else: # l is tuple 'p, q'
                x, y = l.strip().split(', ')
                x, y = int(x), int(y)
                p.add_edge((x, y))
    return ps

def write_posets_to_tex_olavi(ps, fname):
    """Write tex file fname showing forbidden substructures.
    - ps: list of tuples (poset, rem_point, forbidden)
    """
    ps.sort(key = lambda x: x[2]) # sort by forbidden substructure
    with open(fname, 'w') as f:
        f.write("\
\\documentclass[a4paper]{article}\n\
\\usepackage{tikz}\n\
\\usetikzlibrary{arrows}\n\
\\tikzset{->, auto, >=latex'}\n\
\\newcommand*\\intpt{$\\circ$}\n\
\\begin{document}\n\n")
        f.write(f'\
{len(ps)} posets on 7 points which are not gp.\n\
Forbidden substructure in black (removing red part).\n\n')
        count = 0
        for p, rem_point, forbidden in ps:
            count += 1
            f.write('\\bigskip\n\\begin{tikzpicture}\n')
            f.write(f'\\node at (-1,0) {{{count:}}};\n')
            f.write(p.to_tikz_olavi(rem_point, forbidden))
            f.write('\\end{tikzpicture}\n\n')
        f.write('\\end{document}\n')


def testposet():
    t = Poset({1,2,3,4}, {(1,2), (3,4)})
    print('t:', t)
    print('t.step_extensions():', t.step_extensions())
    print('t.interval_extensions():', t.interval_extensions())
    print('t.order_extensions()', t.order_extensions())
    print('t.parallel(t):', t.parallel(t))
    print('t.serial(t):', t.serial(t))
    print('t.is_connected():', t.is_connected(print_counterexample=True))
    print('t.minimals():', t.minimals())
    print('t.maximals():', t.maximals())
    print(t.successors(1))
    t = Poset({1, 2, 3, 4, 5, 6}, {(5, 4), (1, 3), (5, 6), (1, 6), (2, 3), (2, 4)})
    print(t.hasse_diagram())
    return t.to_tikz()
