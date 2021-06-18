#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 12:41:52 2019

@author: uli
"""

from poset import *
from iposet import *
import itertools

def all_posets(n):
    """Generate all posets with n points, up to isomorphism."""
    # handle special cases n=0 and n=1
    if n == 0:
        return {Poset(set(), set())}
    if n == 1:
        return {Poset({1}, set())}
    # now for n>=2
    vertices = set(range(1, n+1))
    alledges = set()
    for x in vertices:
        for y in vertices:
            if x != y:
                alledges.add((x,y))
    #print(alledges)
    allposets = set()
    for numedges in range(len(alledges) + 1):
        for edges_t in itertools.combinations(alledges, numedges):
            #print(edges)
            edges = set(edges_t)
            try:
                p = Poset(vertices, edges)
            except PartialOrderNotAntisymmetricError:
                continue
            seen = False
            for q in allposets:
                if q.is_iso(p):
                    seen = True
                    break
            if not seen:
                allposets.add(p)
    return allposets

def all_connected_posets(n):
    """Generate all connected posets with n points, up to isomorphism."""
    # handle special cases n=0 and n=1
    if n == 0:
        return {Poset(set(), set())}
    if n == 1:
        return {Poset({1}, set())}
    # now for n>=2
    vertices = set(range(1, n+1))
    alledges = set()
    for x in vertices:
        for y in vertices:
            if x != y:
                alledges.add((x,y))
    #print(alledges)
    allposets = set()
    for numedges in range(len(alledges) + 1):
        for edges_t in itertools.combinations(alledges, numedges):
            #print(edges)
            edges = set(edges_t)
            try:
                p = Poset(vertices, edges)
            except PartialOrderNotAntisymmetricError:
                continue
            if not p.is_connected():
                continue
            seen = False
            for q in allposets:
                if p.is_iso(q):
                    seen = True
                    break
            if not seen:
                allposets.add(p)
    return allposets

def sp_posets(n):
    """Generate all sp-posets with n points, up to isomorphism.
    Uses recursion without memoization: slow!
    """
    # base cases
    if n == 0:
        return {Poset(set(), set())}
    if n == 1:
        return {Poset({1}, set())}
    # inductive case
    allposets = set()
    for n1 in range(1, n):
        n2 = n - n1
        ps1 = sp_posets(n1)
        ps2 = sp_posets(n2)
        for p1 in ps1:
            for p2 in ps2:
                p = p1.serial(p2)
                seen = False
                for q in allposets:
                    if p.is_iso(q):
                        seen = True
                        break
                if not seen:
                    allposets.add(p)
                p = p1.parallel(p2)
                seen = False
                for q in allposets:
                    if p.is_iso(q):
                        seen = True
                        break
                if not seen:
                    allposets.add(p)
    return allposets

def sp_posets_memo(n):
    """Generate all sp-posets with n points, up to isomorphism.
    Uses memoization.
    """
    allposets = [set() for _ in range(n+1)]
    allposets[0] = {Poset(set(), set())}
    if n == 0:
        return allposets[n]
    allposets[1] = {Poset({1}, set())}
    if n == 1:
        return allposets[n]
    for numpoints in range(2, n+1):
        for n1 in range(1, numpoints):
            n2 = numpoints - n1
            #print(numpoints, n1, n2)
            ps1 = allposets[n1]
            ps2 = allposets[n2]
            for p1 in ps1:
                for p2 in ps2:
                    p = p1.serial(p2)
                    seen = False
                    for q in allposets[numpoints]:
                        if p.is_iso(q):
                            seen = True
                            break
                    if not seen:
                        allposets[numpoints].add(p)
                    p = p1.parallel(p2)
                    seen = False
                    for q in allposets[numpoints]:
                        if p.is_iso(q):
                            seen = True
                            break
                    if not seen:
                        allposets[numpoints].add(p)
    return allposets[n]

def all_siposets(n):
    """Generate all iposets with n points, with only starting interfaces,
    up to isomorphism.
    """
    # handle special cases n=0 and n=1
    if n == 0:
        return {Iposet(Poset(set(), set()), (), ())}
    if n == 1:
        alliposets = set()
        for s in [(), (1,)]:
            alliposets.add(Iposet(Poset({1}, set()), s, ()))
        return alliposets
    # now for n>=2
    vertices = set(range(1, n+1))
    alledges = set()
    for x in vertices:
        for y in vertices:
            if x != y:
                alledges.add((x,y))
    #print(alledges)
    alliposets = set()
    for numedges in range(len(alledges) + 1):
        for edges_t in itertools.combinations(alledges, numedges):
            #print(edges)
            edges = set(edges_t)
            try:
                p = Poset(vertices, edges)
            except PartialOrderNotAntisymmetricError:
                continue
            for nums in range(n+1):
                for s1 in itertools.combinations(tuple(range(1, n+1)), nums):
                    for s in itertools.permutations(s1):
                        try:
                            ip = Iposet(p, s, ())
                        except StartingInterfaceNotMinimalError:
                            continue
                        seen = False
                        for iq in alliposets:
                            if ip.is_iso(iq):
                                seen = True
                                break
                        if not seen:
                            alliposets.add(ip)
    return alliposets

def all_tiposets(n):
    """Generate all iposets with n points, with only terminating interfaces,
    up to isomorphism.
    """
    # handle special cases n=0 and n=1
    if n == 0:
        return {Iposet(Poset(set(), set()), (), ())}
    if n == 1:
        alliposets = set()
        for t in [(), (1,)]:
            alliposets.add(Iposet(Poset({1}, set()), (), t))
        return alliposets
    # now for n>=2
    vertices = set(range(1, n+1))
    alledges = set()
    for x in vertices:
        for y in vertices:
            if x != y:
                alledges.add((x,y))
    #print(alledges)
    alliposets = set()
    for numedges in range(len(alledges) + 1):
        for edges_t in itertools.combinations(alledges, numedges):
            #print(edges)
            edges = set(edges_t)
            try:
                p = Poset(vertices, edges)
            except PartialOrderNotAntisymmetricError:
                continue
            for numt in range(n+1):
                for t1 in itertools.combinations(tuple(range(1, n+1)), numt):
                    for t in itertools.permutations(t1):
                        try:
                            ip = Iposet(p, (), t)
                        except TerminatingInterfaceNotMaximalError:
                            continue
                        seen = False
                        for iq in alliposets:
                            if ip.is_iso(iq):
                                seen = True
                                break
                        if not seen:
                            alliposets.add(ip)
    return alliposets

def all_iposets(n):
    """Generate all iposets with n points up to isomorphism."""
    # handle special cases n=0 and n=1
    if n == 0:
        return {Iposet(Poset(set(), set()), (), ())}
    if n == 1:
        alliposets = set()
        for s in [(), (1,)]:
            for t in [(), (1,)]:
                alliposets.add(Iposet(Poset({1}, set()), s, t))
        return alliposets
    # now for n>=2
    vertices = set(range(1, n+1))
    alledges = set()
    for x in vertices:
        for y in vertices:
            if x != y:
                alledges.add((x,y))
    #print(alledges)
    alliposets = set()
    for numedges in range(len(alledges) + 1):
        for edges_t in itertools.combinations(alledges, numedges):
            #print(edges)
            edges = set(edges_t)
            try:
                p = Poset(vertices, edges)
            except PartialOrderNotAntisymmetricError:
                continue
            for nums in range(n+1):
                for s1 in itertools.combinations(tuple(range(1, n+1)), nums):
                    for s in itertools.permutations(s1):
                        for numt in range(n+1):
                            for t1 in itertools.combinations(tuple(range(1, n+1)), numt):
                                for t in itertools.permutations(t1):
                                    try:
                                        ip = Iposet(p, s, t)
                                    except (StartingInterfaceNotMinimalError, TerminatingInterfaceNotMaximalError):
                                        continue
                                    seen = False
                                    for iq in alliposets:
                                        if ip.is_iso(iq):
                                            seen = True
                                            break
                                    if not seen:
                                        alliposets.add(ip)
    return alliposets

def gp_posets(n):
    """Return all gp-iposets without interfaces, with n points, up to
    isomorphism.
    """
    alliposets = gp_iposets_memo(n)
    allposets = set()
    for ip in alliposets:
        #print(ip)
        if ip.is_poset():
            #print(ip)
            p = ip._poset
            seen = False
            for q in allposets:
                if q.is_iso(p):
                    seen = True
                    break
            if not seen:
                allposets.add(p)
    return allposets

def gp_iposets(n):
    """Generate all gp-iposets with n points, up to isomorphism.
    Uses recursion without memoization: slow!
    """
    # base cases
    if n == 0:
        return {Iposet(Poset(set(), set()), (), ())}
    if n == 1:
        alliposets = set()
        for s in [(), (1,)]:
            for t in [(), (1,)]:
                alliposets.add(Iposet(Poset({1}, set()), s, t))
        return alliposets
    # inductive case
    alliposets = set()
    for n1 in range(1, n):
        for n2 in range(1, n):
            ips1 = gp_iposets(n1)
            ips2 = gp_iposets(n2)
            for ip1 in ips1:
                for ip2 in ips2:
                    try:
                        ip = ip1.gluing(ip2)
                    except GluingInterfacesDontMatchError:
                        continue
                    if len(ip._poset._vertices) == n:
                        seen = False
                        for iq in alliposets:
                            if ip.is_iso(iq):
                                seen = True
                                break
                        if not seen:
                            alliposets.add(ip)
    for n1 in range(1, n):
        n2 = n - n1
        ips1 = gp_iposets(n1)
        ips2 = gp_iposets(n2)
        for ip1 in ips1:
            for ip2 in ips2:
                ip = ip1.parallel(ip2)
                seen = False
                for iq in alliposets:
                    if ip.is_iso(iq):
                        seen = True
                        break
                if not seen:
                    alliposets.add(ip)
    return alliposets

def gp_iposets_memo(n):
    """Generate all gp-iposets with n points, up to isomorphism.
    Uses memoization
    """
    alliposets = [set() for _ in range(n+1)]
    alliposets[0] = {Iposet(Poset(set(), set()), (), ())}
    if n == 0:
        return alliposets[0]
    for s in [(), (1,)]:
        for t in [(), (1,)]:
            alliposets[1].add(Iposet(Poset({1}, set()), s, t))
    if n ==1:
        return alliposets[1]
    for numpoints in range(2, n+1):
        for n1 in range(1, numpoints):
            for n2 in range(1, numpoints):
                #print(numpoints, n1, n2)
                ips1 = alliposets[n1]
                ips2 = alliposets[n2]
                #print(ips1, ips2)
                for ip1 in ips1:
                    for ip2 in ips2:
                        try:
                            ip = ip1.gluing(ip2)
                        except GluingInterfacesDontMatchError:
                            continue
                        ipsize = len(ip._poset._vertices)
                        if ipsize <= n:
                            seen = False
                            for iq in alliposets[ipsize]:
                                if iq.is_iso(ip):
                                    seen = True
                                    break
                            if not seen:
                                alliposets[ipsize].add(ip)
        for n1 in range(1, numpoints):
            n2 = numpoints - n1
            #print(numpoints, n1, n2)
            ips1 = alliposets[n1]
            ips2 = alliposets[n2]
            for ip1 in ips1:
                for ip2 in ips2:
                    ip = ip1.parallel(ip2)
                    seen = False
                    for iq in alliposets[numpoints]:
                        if iq.is_iso(ip):
                            seen = True
                            break
                    if not seen:
                        alliposets[numpoints].add(ip)
    return alliposets[n]

def gp_iposets_by_st(n):
    """Generate all gp-iposets with n points, up to isomorphism.
    Differentiate by number of interfaces
    """
    alliposets = [dict() for _ in range(n+1)]
    alliposets[0] = {(0, 0): {Iposet(Poset(set(), set()), (), ())}}
    if n == 0:
        return alliposets[0]
    p = Poset({1}, set())
    alliposets[1] = {(0, 0): {Iposet(p, (), ())},
                     (1, 0): {Iposet(p, (1,), ())},
                     (0, 1): {Iposet(p, (), (1,))},
                     (1, 1): {Iposet(p, (1,), (1,))}}
    if n == 1:
        return alliposets[1]
    for numpoints in range(2, n+1):
        for n1 in range(1, numpoints):
            for n2 in range(1, numpoints):
                for t1 in range(max(n1+n2-n,0), min(n1,n2)+1):
                    s2 = t1
                    ipsize = n1 + n2 - t1
                    for s1 in range(n1+1):
                        for t2 in range(n2+1):
                            #print(numpoints, n1, n2)
                            ips1 = alliposets[n1][(s1, t1)]
                            ips2 = alliposets[n2][(s2, t2)]
                            #print(ips1, ips2)
                            for ip1 in ips1:
                                for ip2 in ips2:
                                    ip = ip1.gluing(ip2)
                                    if (s1, t2) not in alliposets[ipsize]:
                                        alliposets[ipsize][(s1, t2)] = {ip}
                                    else:
                                        seen = False
                                        for iq in alliposets[ipsize][(s1, t2)]:
                                            if iq.is_iso(ip):
                                                seen = True
                                                break
                                        if not seen:
                                            alliposets[ipsize][(s1, t2)].add(ip)
        for n1 in range(1, numpoints):
            n2 = numpoints - n1
            for s1 in range(n1+1):
                for t1 in range(n1+1):
                    for s2 in range(n2+1):
                        for t2 in range(n2+1):
                            #print(numpoints, n1, n2)
                            ips1 = alliposets[n1][(s1, t1)]
                            ips2 = alliposets[n2][(s2, t2)]
                            s = s1 + s2
                            t = t1 + t2
                            for ip1 in ips1:
                                for ip2 in ips2:
                                    ip = ip1.parallel(ip2)
                                    if (s, t) not in alliposets[numpoints]:
                                        alliposets[numpoints][(s, t)] = {ip}
                                    else:
                                        seen = False
                                        for iq in alliposets[numpoints][(s, t)]:
                                            if iq.is_iso(ip):
                                                seen = True
                                                break
                                        if not seen:
                                            alliposets[numpoints][(s, t)].add(ip)
    return alliposets[n]

def gpc_iposets_by_st(n):
    """Generate all very weakly connected gp-iposets with n points, up to isomorphism.
    Differentiate by number of interfaces
    """
    alliposets = [dict() for _ in range(n+1)]
    alliposets[0] = {(0, 0): {Iposet(Poset(set(), set()), (), ())}}
    if n == 0:
        return alliposets[0]
    p = Poset({1}, set())
    alliposets[1] = {(0, 0): {Iposet(p, (), ())},
                     (1, 0): {Iposet(p, (1,), ())},
                     (0, 1): {Iposet(p, (), (1,))},
                     (1, 1): {Iposet(p, (1,), (1,))}}
    if n == 1:
        return alliposets[1]
    for numpoints in range(2, n):
        for n1 in range(1, numpoints):
            for n2 in range(1, numpoints):
                for t1 in range(max(n1+n2-n,0), min(n1,n2)+1):
                    s2 = t1
                    ipsize = n1 + n2 - t1
                    for s1 in range(n1+1):
                        for t2 in range(n2+1):
                            #print(numpoints, n1, n2)
                            ips1 = alliposets[n1][(s1, t1)]
                            ips2 = alliposets[n2][(s2, t2)]
                            #print(ips1, ips2)
                            for ip1 in ips1:
                                for ip2 in ips2:
                                    ip = ip1.gluing(ip2)
                                    if (s1, t2) not in alliposets[ipsize]:
                                        alliposets[ipsize][(s1, t2)] = {ip}
                                    else:
                                        seen = False
                                        for iq in alliposets[ipsize][(s1, t2)]:
                                            if iq.is_iso(ip):
                                                seen = True
                                                break
                                        if not seen:
                                            alliposets[ipsize][(s1, t2)].add(ip)
        for n1 in range(1, numpoints):
            n2 = numpoints - n1
            for s1 in range(n1+1):
                for t1 in range(n1+1):
                    for s2 in range(n2+1):
                        for t2 in range(n2+1):
                            #print(numpoints, n1, n2)
                            ips1 = alliposets[n1][(s1, t1)]
                            ips2 = alliposets[n2][(s2, t2)]
                            s = s1 + s2
                            t = t1 + t2
                            for ip1 in ips1:
                                for ip2 in ips2:
                                    ip = ip1.parallel(ip2)
                                    if (s, t) not in alliposets[numpoints]:
                                        alliposets[numpoints][(s, t)] = {ip}
                                    else:
                                        seen = False
                                        for iq in alliposets[numpoints][(s, t)]:
                                            if iq.is_iso(ip):
                                                seen = True
                                                break
                                        if not seen:
                                            alliposets[numpoints][(s, t)].add(ip)
    for numpoints in range(n, n+1):
        """last step only uses gluing; no parallel"""
        for n1 in range(1, numpoints):
            for n2 in range(1, numpoints):
                for t1 in range(max(n1+n2-n,0), min(n1,n2)+1):
                    s2 = t1
                    ipsize = n1 + n2 - t1
                    for s1 in range(n1+1):
                        for t2 in range(n2+1):
                            #print(numpoints, n1, n2)
                            ips1 = alliposets[n1][(s1, t1)]
                            ips2 = alliposets[n2][(s2, t2)]
                            #print(ips1, ips2)
                            for ip1 in ips1:
                                for ip2 in ips2:
                                    ip = ip1.gluing(ip2)
                                    if (s1, t2) not in alliposets[ipsize]:
                                        alliposets[ipsize][(s1, t2)] = {ip}
                                    else:
                                        seen = False
                                        for iq in alliposets[ipsize][(s1, t2)]:
                                            if iq.is_iso(ip):
                                                seen = True
                                                break
                                        if not seen:
                                            alliposets[ipsize][(s1, t2)].add(ip)
    return alliposets[n]

def gpi_decompositions(n):
    """Return decompositions for all gp-iposets with n points.
    Find decompositions with fewest interfaces.
    """
    alliposets = [{} for _ in range(n+1)]
    alliposets[0] = {Iposet(Poset(set(), set()), (), ()): (None, None, None)}
    if n == 0:
        return alliposets[0]
    for s in [(), (1,)]:
        for t in [(), (1,)]:
            alliposets[1][Iposet(Poset({1}, set()), s, t)] = (None, None, None)
    if n == 1:
        return alliposets[1]
    for numpoints in range(2, n+1):
        for n1 in range(1, numpoints):
            for n2 in range(1, numpoints):
                #print(numpoints, n1, n2)
                ips1 = alliposets[n1]
                ips2 = alliposets[n2]
                #print(ips1, ips2)
                for ip1 in ips1:
                    for ip2 in ips2:
                        try:
                            ip = ip1.gluing(ip2)
                        except GluingInterfacesDontMatchError:
                            continue
                        ipsize = len(ip._poset._vertices)
                        if ipsize <= n:
                            seen = False
                            for iq in alliposets[ipsize]:
                                if iq.is_iso(ip):
                                    seen = True
                                    # do we have fewer interfaces in the comp?
                                    if alliposets[ipsize][iq][0] is not None:
                                        if len(ip1._t) < len(alliposets[ipsize][iq][0]._t):
                                            del alliposets[ipsize][iq]
                                            alliposets[ipsize][ip] = (ip1, 'glue', ip2)
                                    break
                            if not seen:
                                alliposets[ipsize][ip] = (ip1, 'glue', ip2)
        for n1 in range(1, numpoints):
            n2 = numpoints - n1
            #print(numpoints, n1, n2)
            ips1 = alliposets[n1]
            ips2 = alliposets[n2]
            for ip1 in ips1:
                for ip2 in ips2:
                    ip = ip1.parallel(ip2)
                    seen = False
                    for iq in alliposets[numpoints]:
                        if iq.is_iso(ip):
                            seen = True
                            break
                    if not seen:
                        alliposets[numpoints][ip] = (ip1, 'par', ip2)
    return alliposets[n]

def gpc_decompositions(n):
    """Return decompositions for all connected gp-posets with n points.
    Timing: n=4 : .5s
            n=5 : 50s
            n=6 : 280min
    """
    allipos = gpi_decompositions(n)
    retd = {}
    for ip in allipos:
        if ip.is_poset() and ip.is_connected():
            #print(ip, ip.to_poset())
            retd[ip] = (allipos[ip][0], allipos[ip][2])
    return retd

def gpc_decompositions_to_file(n, fname):
    """Write decompositions for all connected gp-posets with n points to file.
    """
    dcs = gpc_decompositions(n)
    with open(fname, 'w') as f:
        for p in dcs:
            f.write(str(p) + ': ' + str(dcs[p]) + '\n')

def convert_decomp_to_tikz(fprefix):
    """Convert decompositions in `fprefix` to tikz code.
    Reads `fprefix`.txt, writes `fprefix`.tex."""
    decomps = dict()
    with open(fprefix + '.txt') as f:
        for line in f:
            if line != '\n':
                ip, decomp = line.split(':')
                decomps[eval(ip)] = eval(decomp)
    with open(fprefix + '.tex', 'w') as f:
        f.write('% generated by ramics.py\n')
        count = 0
        for ip in decomps:
            count += 1
            f.write('{} &\n'.format(count))
            f.write('\\begin{tikzpicture}\n')
            f.write(ip.to_tikz())
            f.write('\\end{tikzpicture}\n')
            for i in range(2):
                f.write('&\n')
                f.write('\\begin{tikzpicture}\n')
                f.write(decomps[ip][i].to_tikz())
                f.write('\\end{tikzpicture}\n')
            f.write('\\\\\\hline\n')


def gpc_decompositions_to_tikz(n, fname):
    """Write decompositions for all connected gp-posets with n points to file,
    as a tabular of tikzpictures."""
    assert n >= 2
    dcs = gpc_decompositions(n)
    with open(fname, 'w') as f:
        f.write('% generated by ramics.py\n')
        count = 0
        for ip in dcs:
            count += 1
            f.write('{} &\n'.format(count))
            f.write('\\begin{tikzpicture}\n')
            f.write(ip.to_tikz())
            f.write('\\end{tikzpicture}\n')
            for i in range(2):
                f.write('&\n')
                f.write('\\begin{tikzpicture}\n')
                f.write(dcs[ip][i].to_tikz())
                f.write('\\end{tikzpicture}\n')
            f.write('\\\\\\hline\n')

def all_posets_to_file():
    """timing:
    range(5) : .14s
    range(6) : 61s
    range(7) : > 3h
    """
    for n in range(5):
        write_to_file(all_posets(n), 'posets'+str(n)+'.txt')

def gp_posets_to_file():
    """timing:
    range(5) : .5s
    range(6) : 47s
    range(7) : 5h
    """
    for n in range(5):
        write_to_file(gp_posets(n), 'gpposets'+str(n)+'.txt')

def connected_gp_posets(n):
    """Return all connected gp-posets with n points.
    Needs file of gp-posets created by gp_posets_to_file().
    """
    assert n <= 6
    gps = read_posets_from_file('gpposets{}.txt'.format(n))
    rets = set()
    for p in gps:
        if p.is_connected():
            rets.add(p)
    return rets

def missing_posets6():
    """Find the 5 posets on 6 points which are not gp.
    Timing: ca. 2h
    """
    gps = read_posets_from_file('gpposets6.txt')
    vertices = set(range(1, 7))
    alledges = set()
    for x in vertices:
        for y in vertices:
            if x != y:
                alledges.add((x,y))
    #print(alledges)
    missing = set()
    for numedges in range(len(alledges) + 1):
        for edges_t in itertools.combinations(alledges, numedges):
            edges = set(edges_t)
            #print(edges)
            try:
                p = Poset(vertices, edges)
            except PartialOrderNotAntisymmetricError:
                continue
            seen = False
            for q in gps:
                if q.is_iso(p):
                    seen = True
                    break
            if not seen: # it's new! but is it iso to some other new?
                for q in missing:
                    if q.is_iso(p):
                        seen = True
                        break
                if not seen:
                    print(p)
                    missing.add(p)
    write_to_file(missing, 'posets6nogp.txt')

def parse_argv(args):
    """Return int(argv[1]) iff this was the only arg and the int() call
    succeeds; otherwise return None."""
    if len(args) == 2:
        try:
            return int(args[1])
        except ValueError:
            pass
    return None



