using LightGraphs
using Combinatorics
using StaticArrays
import LightGraphs.transitiveclosure
import LightGraphs.intersect
import Base.string
import Base.Threads.@spawn
import Base.==
import Base.reverse
import LightGraphs.vertices
import LightGraphs.nv
import LightGraphs.edges
import LightGraphs.ne
include("LightPosets.jl")

"""
Library for operations on iposets: posets with interfaces.
Iposets are structs consisting of two tuples s and t and a poset P.
s: starting interface
t: terminating interface
Nodes in s are minimal in P, nodes in t are maximal in P.
!! Enforcement of these properties is up to the caller !!
"""

const MAXNODES = 16

struct GluingInterfacesDontMatchError <: Exception end

struct Iposet
    s::Tuple{Vararg{Int}}
    t::Tuple{Vararg{Int}}
    poset::Poset
#    Iposet(s, t, p) = ismin(p, s) && ismax(p, t) ? new(s, t, p) : error("Interfaces not minimal/maximal")
#    Iposet(s, t, p) = new(s, t, p)
end

"""Type for vertex profiles."""
#Vprof = Array{Tuple{Int, Int}}
Vprof = MVector{MAXNODES, Tuple{Int, Int}}
@inline function initvprof()
    #return Vprof(undef, n)
    return Vprof((0, 0) for _ in 1:MAXNODES)
end

"""Identity on n points."""
function Idiposet(n)
    st = Tuple(x for x in 1:n)
    return Iposet(st, st, Poset(n))
end

"""The four singletons."""
function singletons()
    return (Iposet("1 0"),
            Iposet("1 0 0"),
            Iposet("1 0 - 0"),
            Iposet("1 0 0 0"))
end

"""Compute inhash of vertex v in digraph g:
number of incoming edges + 16 * their inhashes"""
function inhash(g::Poset, v::Int)
    n = inneighbors(g, v)
    nhashes = 0
    for neighbour in n
        nhashes += inhash(g, neighbour)
    end
    return length(n) + nhashes << 4
end

"""Compute outhash of vertex v in digraph g:
number of outgoing edges + 16 * their outhashes"""
function outhash(g::Poset, v::Int)
    n = outneighbors(g, v)
    nhashes = 0
    for neighbour in n
        nhashes += outhash(g, neighbour)
    end
    return length(n) + nhashes << 4
end

"""Check if tuple is bijection."""
function bijective(a::Tuple{Vararg{Int}})
    n = length(a)
    seen = falses(n)
    @inbounds for i in 1:n
        if a[i] > n
            return false
        end
        seen[a[i]] = true
    end
    return !(false in seen)
end

"""Check if tuple a is an isomorphism from g to s.
Assumes that g and s have same number of vertices and same number
of edges, and that a is a bijection
"""
function isomorphic(g::Poset, s::Poset, a::Tuple{Vararg{Int}})
    for e in edges(g)
        if !(has_edge(s, a[src(e)], a[dst(e)]))
            return false
        end
    end
    return true
end

"""Check if posets g and s are isomorphic."""
function isomorphic(g::Poset, s::Poset)
    #Start with the easy stuff
    if g == s
        return true
    end
    n = nv(g)
    if n != nv(s) || ne(g) != ne(s)
        return false
    end
    #Make arrays of node invariants
    gvprof = initvprof()
    svprof = initvprof()
    for v in 1:n
        gvprof[v] = (inhash(g, v), outhash(g, v))
        svprof[v] = (inhash(s, v), outhash(s, v))
    end
    #Now call the iso checker with the profiles
    return isomorphic(g, gvprof, s, svprof)
end

"""Check if g and s are isomorphic, using precomputed vertex profiles."""
function isomorphic(g::Poset, gvprof::Vprof, 
                    s::Poset, svprof::Vprof)
    #Start with the easy stuff
    if g == s && gvprof == svprof
        return true
    end
    n = nv(g)
    if n != nv(s) || ne(g) != ne(s)
        return false
    end
    #Check that there exists a vertex profile bijection
    used = zeros(Bool, n)
    @inbounds for i in 1:n
        found = false
        for j in 1:n
            if gvprof[i] == svprof[j] && !used[j]
                used[j] = true
                found = true
                break
            end
        end
        if !found
            return false
        end
    end
    # #Break for testing
    # println("Found hash iso:\n", string(g), '\n', string(s))
    # return false
    # #END Break
    #Construct an array node -> list of possible images in iso
    targets = Array{Array{Int}}(undef, n)
    @inbounds for v in 1:n
        targets[v] = []
    end
    @inbounds for v in 1:n
        for i in 1:n
            if gvprof[v] == svprof[i]
                push!(targets[v], i)
            end
        end
    end
    #Check all target permutations if they are isos
    for pos_isom in Iterators.product(targets...)
        if bijective(pos_isom) && isomorphic(g, s, pos_isom)
            return true
        end
    end
    return false
end

"""Check if s is isomorphic to an induced subgraph of g.
Return a tuple (Bool, Array{Int}) to signify success and the
subset of nodes that give the induced subgraph"""
function subgraph(g::Poset, s::Poset)
    n = nv(s)
    svprof = initvprof()
    for v in 1:n
        svprof[v] = (inhash(s, v), outhash(s, v))
    end
    for nodes in combinations(1:nv(g), n)
        subgraph = induced_subgraph(g, nodes)[1]
        subvprof = initvprof()
        for v in 1:n
            subvprof[v] = (inhash(subgraph, v), outhash(subgraph, v))
        end
        if isomorphic(subgraph, subvprof, s, svprof)
            return (true, nodes)
        end
    end
    return (false, [])
end

"""Iposet with empty interfaces."""
function Iposet(p::Poset)
    return Iposet((), (), p)
end

"""Iposets are equal iff their components are."""
function ==(p::Iposet, q::Iposet)
    return p.s == q.s && p.t == q.t && p.poset == q.poset
end

"""Vertices."""
function vertices(p::Iposet)
    return vertices(p.poset)
end
function nv(p::Iposet)
    return nv(p.poset)
end

"""Edges."""
function edges(p::Iposet)
    return edges(p.poset)
end
function ne(p::Iposet)
    return ne(p.poset)
end

"""Check if iposet is weakly connected."""
function isconnected(g::Iposet)
    return is_connected(g.poset)
end

"""
Properties of iposets:
leftwinkowski
rightwinkowski
winkowski
discrete
starter
terminator
almostconnected
intfconsistent
"""

"""Check if all min points are in the starting interface."""
function leftwinkowski(g::Iposet)
    for v in min(g.poset)
        if !(v in g.s)
            return false
        end
    end
    return true
end

"""Check if all max points are in the terminating interface."""
function rightwinkowski(g::Iposet)
    for v in max(g.poset)
        if !(v in g.t)
            return false
        end
    end
    return true
end

"""Check if g is left and right winkowski."""
function winkowski(g::Iposet)
    return leftwinkowski(g) && rightwinkowski(g)
end

"""Check if p is discrete."""
function discrete(p::Union{Poset, Iposet})
    return ne(p) == 0
end
# function discrete(g::Iposet)
#     return ne(g) == 0
# end

"""Check if g is starter."""
function starter(g::Iposet)
    return discrete(g) && rightwinkowski(g)
end

"""Check if g is terminator."""
function terminator(g::Iposet)
    return discrete(g) && leftwinkowski(g)
end

"""Check if iposet is "almost connected": one connected
component plus possibly some isolated points."""
function isalmostconnected(p::Iposet)
    big_one = false
    for comp in weakly_connected_components(p.poset)
        if length(comp) > 1 && big_one == false
            big_one = true
        elseif length(comp) <= 1
            continue
        else
            return false
        end
    end
    return true
end

"""Check if iposet is interface consistent."""
function intfconsistent(p::Iposet)
    st = Array{Union{Int,Nothing}}(undef, length(p.s))
    # find index of each source in targets
    for i in eachindex(st)
        st[i] = findfirst(x -> x==p.s[i], p.t)
    end
    #println(st)
    return issorted((x for x in st if !isnothing(x)))
end

"""
Operations on single iposets:
transitiveclosure
hassediagram
reverse
"""

"""Transitive closure of iposet. No reflexive edges added"""
function transitiveclosure(g::Iposet)
    return Iposet(g.s, g.t, transitiveclosure(g.poset, false))
end
function transitiveclosure!(g::Iposet)
    transitiveclosure!(g.poset, false)
end

"""Iposet version of hassediagram."""
function hassediagram(g::Iposet)
    return Iposet(g.s, g.t, transitivereduction(g.poset))
end
function hassediagram!(g::Iposet)
    transitivereduction!(g.poset)
end

"""Return reverse of iposet (all edges reversed, sources and targets
inverted).
"""
function reverse(g::Iposet)
    return Iposet(g.t, g.s, reverse(g.poset))
end
function reverse!(g::Iposet)
    reverse!(g.poset)
    g.s, g.t = g.t, g.s
end
function reverse(tip::Tuple{Iposet, Vprof})
    vproflen = length(tip[2])
    nvprof = initvprof()
    @inbounds for j in 1:vproflen
        nvprof[j] = (tip[2][j][2], tip[2][j][1])
    end
    return reverse(tip[1]), nvprof
end

"""
Operations on multiple (i)posets:
glue
parallel
wglue
riglue
liglue
"""

"""Glue iposets g and p along their interfaces."""
function glue(p::Iposet, q::Iposet)
    il = length(p.t)
    il != length(q.s) && throw(GluingInterfacesDontMatchError)
    np = nv(p)
    nq = nv(q)
    #resulting poset, empty for now
    rposet = Poset(np + nq - il)
    #add edges from p
    for e in edges(p)
        add_edge!(rposet, e)
    end
    #add edges from q
    qmap = zeros(Int, nq) #qmap[vertex in q.poset] = vertex in rposet
    @inbounds for i in 1:il
        qmap[q.s[i]] = p.t[i]
    end
    counter = 1
    @inbounds for v in 1:nq
        if qmap[v] == 0
            qmap[v] = np + counter
            counter += 1
        end
    end
    for e in edges(q)
        add_edge!(rposet, qmap[src(e)], qmap[dst(e)])
    end
    #add edges from p to q
    @inbounds for v1 in vertices(p)
        if !(v1 in p.t)
            for v2 in (np + 1):(np + nq - il)
                add_edge!(rposet, v1, v2)
            end
        end
    end
    #construct target interface
    tl = length(q.t)
    rt = Array{Int}(undef, tl)
    @inbounds for i in 1:tl
        rt[i] = qmap[q.t[i]]
    end
    return Iposet(p.s, Tuple(rt), rposet)
end
function glue(ips::Vector{Iposet})
    length(ips) == 0 && return Idiposet(0)
    length(ips) == 1 && return ips[1]
    return glue(ips[1], glue(ips[2:end]))
end

"""Parallel composition of iposets g and s."""
function parallel(p::Iposet, q::Iposet)
    #TODO use LightGraphs.blockdiag ?
    psize = nv(p.poset)
    #resulting poset, empty for now
    rposet = Poset(psize + nv(q.poset))
    #add edges from p
    for e in edges(p.poset)
        add_edge!(rposet, e)
    end
    #add edges from q
    for e in edges(q.poset)
        add_edge!(rposet, src(e) + psize, dst(e) + psize)
    end
    #construct new source interface
    ns = Array{Int}(undef, length(p.s) + length(q.s))
    for i in 1:length(p.s)
        ns[i] = p.s[i]
    end
    for i in 1:length(q.s)
        ns[i+length(p.s)] = q.s[i] + psize
    end
    nt = Array{Int}(undef, length(p.t) + length(q.t))
    for i in 1:length(p.t)
        nt[i] = p.t[i]
    end
    for i in 1:length(q.t)
        nt[i+length(p.t)] = q.t[i] + psize
    end
    return Iposet(Tuple(ns), Tuple(nt), rposet)
end
function parallel(ips::Array{Iposet})
    length(ips) == 0 && return Idiposet(0)
    length(ips) == 1 && return ips[1]
    return parallel(ips[1], parallel(ips[2:end]))
end

"""Return all non-isomorphic Winkoswki gluings of two posets."""
function wglue(p::Poset, q::Poset)
    maxp = collect(max(p))
    minq = collect(min(q))
    il = length(maxp)
    il != length(minq) && throw(GluingInterfacesDontMatchError)
    discrete(p) && return (q,)
    discrete(q) && return (p,)
    np = nv(p)
    nq = nv(q)
    res = Array{Tuple{Poset, Vprof}}(undef, 0)
    for maxpp in permutations(maxp)
        #resulting poset, empty for now
        r = Poset(np + nq - il)
        #add edges from p
        for e in edges(p)
            add_edge!(r, e)
        end
        #add edges from q
        qmap = zeros(Int, nq) #qmap[vertex in q] = vertex in r
        @inbounds for i in 1:il
            qmap[minq[i]] = maxpp[i]
        end
        counter = 1
        @inbounds for v in 1:nq
            if qmap[v] == 0
                qmap[v] = np + counter
                counter += 1
            end
        end
        for e in edges(q)
            add_edge!(r, qmap[src(e)], qmap[dst(e)])
        end
        #add edges from p to q
        @inbounds for v1 in vertices(p)
            if !(v1 in maxp)
                for v2 in (np + 1):(np + nq - il)
                    add_edge!(r, v1, v2)
                end
            end
        end
        pushuptoiso!(res, r)
    end
    return (x[1] for x in res)
end

"""Right-interior gluing of two posets."""
function riglue(p::Poset, q::Poset)
    np = nv(p)
    nq = nv(q)
    #resulting poset, empty for now
    r = Poset(np + nq)
    #add edges from p
    for e in edges(p)
        add_edge!(r, e)
    end
    #add edges from q
    for e in edges(q)
        add_edge!(r, src(e)+np, dst(e)+np)
    end
    #add edges from p to q
    for v2 in vertices(q)
        if !isminimal(q, v2)
            for v1 in vertices(p)
                add_edge!(r, v1, v2+np)
            end
        end
    end
    return r
end

"""Left-interior gluing of two posets."""
function liglue(p::Poset, q::Poset)
    np = nv(p)
    nq = nv(q)
    #resulting poset, empty for now
    r = Poset(np + nq)
    #add edges from p
    for e in edges(p)
        add_edge!(r, e)
    end
    #add edges from q
    for e in edges(q)
        add_edge!(r, src(e)+np, dst(e)+np)
    end
    #add edges from p to q
    for v1 in vertices(p)
        if !ismaximal(p, v1)
            for v2 in vertices(q)
                add_edge!(r, v1, v2+np)
            end
        end
    end
    return r
end

"""
Functions for loading and saving iposets from files:
string(Iposet)
Iposet(string)
saveiposets
loadiposets
"""

"""Convert iposet to string representation.
Format: n k e1 .. ek [s [t]]
n: number of points, numbered from 0 to n-1 in hexadecimal
k: number of edges in Hasse diagram
e1 .. ek: edges: format ab: a = source point, b = target point
s: sources: format - or a+:
   a+: string containing all points in source interface
   -: empty source interface; only if t exists
t: targets: format a+ like for s (may not be -)
"""
function string(p::Iposet)
    g = hassediagram(p)
    resl = Array{String}(undef, 0)
    push!(resl, string(nv(g), base=16))
    push!(resl, string(ne(g), base=16))
    for e in edges(g)
        push!(resl, string(src(e)-1, base=16) * string(dst(e)-1, base=16))
    end
    if length(g.s) == 0 && length(g.t) > 0
        push!(resl, "-")
    end
    if length(g.s) > 0
        sl = Array{String}(undef, 0)
        for s in g.s
            push!(sl, string(s-1, base=16))
        end
        push!(resl, join(sl))
    end
    if length(g.t) > 0
        tl = Array{String}(undef, 0)
        for t in g.t
            push!(tl, string(t-1, base=16))
        end
        push!(resl, join(tl))
    end
    return join(resl, " ")
end
function string(p::Poset)
    return string(Iposet(p))
end

"""Convert string to iposet. See above for format.
Old file format has node numbers starting at 1 instead of 0
Use hex numbers if hex; otherwise decimal
If comma, then numbers in edges are separated by ","
"""
function Iposet(line::String; hex=true, comma=false)
    oldformat = false
    addone::Int = (oldformat ? 0 : 1)
    hd::Int = (hex ? 16 : 10)
    nums = split(line)
    np = parse(Int, nums[1], base=hd)
    ne = parse(Int, nums[2], base=hd)
    dg = Poset(np)
    for edges_num in nums[3:ne+2]
        if comma
            edges_num = split(edges_num, ",")
        end
        add_edge!(dg, parse(Int, edges_num[1], base=hd)+addone, parse(Int, edges_num[2], base=hd)+addone)
    end
    s::Tuple{Vararg{Int}} = ()
    t::Tuple{Vararg{Int}} = ()
    if length(nums) == ne+2 # no sources neither targets
        nothing
    elseif length(nums) == ne+3 # sources, but no targets
        s = Tuple([parse(Int, x, base=hd)+addone for x in nums[ne+3]])
    elseif length(nums) == ne+4 && nums[ne+3] == "-" # targets, but no sources
        t = Tuple([parse(Int, x, base=hd)+addone for x in nums[ne+4]])
    else # sources and targets
        s = Tuple([parse(Int, x, base=hd)+addone for x in nums[ne+3]])
        t = Tuple([parse(Int, x, base=hd)+addone for x in nums[ne+4]])
    end
    return transitiveclosure(Iposet(s, t, dg))
end

"""Save array of iposets to file."""
function saveiposets(a::Array{T}, filename::String) where T <: Union{Poset, Iposet}
    open(filename, "w") do file
        for ips in a
            write(file, string(ips)*"\n")
        end
    end
end

"""Save array of iposets plus vprofs to file, dropping vprofs."""
function saveiposets(a::Array{Tuple{Iposet, Vprof}}, filename::String)
    open(filename, "w") do file
        for tips in a
            write(file, string(tips[1])*"\n")
        end
    end
end

"""Iterator which returns each iposet in file, ignoring comments."""
function eachiposet(filename::String; hex=true, comma=false)
    return (Iposet(x, hex=hex, comma=comma)
            for x in eachline(filename) if !occursin(r"^\s*(?:#|$)", x))
end

"""Load array of iposets from file."""
function loadiposets(filename::String; hex=true, comma=false)
    res = Array{Iposet}(undef, 0)
    for ip in eachiposet(filename, hex=hex, comma=comma)
        push!(res, ip)
    end
    return res
end
function loadiposets!(res::Array{Iposet}, filename::String; hex=true, comma=false)
    for ip in eachiposet(filename, hex=hex, comma=comma)
        push!(res, ip)
    end
end

"""
Isomorphism checking for iposets:
inhash
outhash
isomorphic(g, s)
isomorphic(g, g_vprof, s, s_vprof)
"""

"""Inhash and outhash for iposets: like for posets, but interfaces get
unique negative numbers.
"""
function inhash(g::Iposet, v::Int)
    if v in g.t && !(v in g.s)
        return 0
    end
    @inbounds for i in 1:length(g.s)
        if g.s[i] == v
            return -i
        end
    end
    return inhash(g.poset, v)
end

function outhash(g::Iposet, v::Int)
    if v in g.s && !(v in g.t)
        return 0
    end
    @inbounds for i in 1:length(g.t)
        if g.t[i] == v
            return -i
        end
    end
    return outhash(g.poset, v)
end

"""Check if iposets g and s are isomorphic."""
function isomorphic(g::Iposet, s::Iposet)
    #Start with the easy stuff
    if g == s
        return true
    end
    n = nv(g)
    if n != nv(s) || ne(g) != ne(s)
        return false
    end
    if length(g.s) != length(s.s) || length(g.t) != length(s.t)
        return false
    end
    #Make arrays of node invariants
    gvprofs = initvprof()
    svprofs = initvprof()
    for v in 1:n
        gvprofs[v] = (inhash(g, v), outhash(g, v))
        svprofs[v] = (inhash(s, v), outhash(s, v))
    end
    #Now call the iso checker with the profiles
    return isomorphic(g, gvprofs, s, svprofs)
end

"""Check if g and s are isomorphic, using precomputed vertex profiles."""
function isomorphic(g::Iposet, gvprofs::Vprof, 
                    s::Iposet, svprofs::Vprof)
    #Start with the easy stuff
    if g == s && gvprofs == svprofs
        return true
    end
    if length(g.s) != length(s.s) || length(g.t) != length(s.t)
        return false
    end
    #Interface info is in the profiles,
    #so we just call the corresponding poset iso
    return isomorphic(g.poset, gvprofs, s.poset, svprofs)
end

function containsuptoiso(iqs, ip)
    n = nv(ip)
    vprof = initvprof()
    @inbounds for v in 1:n
        vprof[v] = (inhash(ip, v), outhash(ip, v))
    end
    for tiq in iqs
        if isomorphic(tiq[1], tiq[2], ip, vprof)
            return true
        end
    end
    return false
end

"""push (i)poset (and computed vprof) into array of (i)posets and
vprofs, unless array already contains an isomorph.  Return false if
array contained an isomorph, otherwise true
"""
#function pushuptoiso!(iqs::Array{Tuple{Iposet, Vprof}}, ip)
function pushuptoiso!(iqs, ip)
    n = nv(ip)
    vprof = initvprof()
    @inbounds for v in 1:n
        vprof[v] = (inhash(ip, v), outhash(ip, v))
    end
    pushuptoiso!(iqs, ip, vprof)
end
function pushuptoiso!(iqs, ip, vprof::Vprof)
    for tiq in iqs
        if isomorphic(tiq[1], tiq[2], ip, vprof)
            return false
        end
    end
    push!(iqs, (ip, vprof))
    return true
end
#now with a lock passed as arg, to tighten locking
function pushuptoiso!(iqs, ip, mylock::ReentrantLock)
    n = nv(ip)
    vprof = initvprof()
    @inbounds for v in 1:n
        vprof[v] = (inhash(ip, v), outhash(ip, v))
    end
    lock(mylock)
    pushuptoiso!(iqs, ip, vprof)
    unlock(mylock)
    return true
end

"""Reduce an array of iposets up to iso."""
function reduceiso(ips::Array{Iposet})
    vprofs = Array{Vprof}(undef, length(ips))
    for i in eachindex(ips)
        np = nv(ips[i])
        vprofs[i] = initvprof()
        for v in 1:np
            vprofs[i][v] = (inhash(ips[i], v), outhash(ips[i], v))
        end
    end
    res = Array{Tuple{Iposet, Vprof}}(undef, 0)
    for i in eachindex(ips)
        pushuptoiso!(res, ips[i], vprofs[i])
    end
    return (x[1] for x in res)
end
function reduceiso(ps::T) where T <: Union{Array{Poset}, Array{Poset{Int}}}
    return (y.poset for y in reduceiso([Iposet(x) for x in ps]))
end

"""Return n-th power of an iterator.
TODO Is there a better way to do this?
"""
# function power_old(itr, n::Int)
#     n < 1 && error("Only positive powers defined.")
#     n == 1 && return itr
#     # return (flatten(x) for x in product(power(itr, n-1), itr))
#     return product(power(itr, n-1), itr)
# end
function power(itr, n::Int)
    return product(repeated(itr, n)...)
end

"""Return order dimension of poset."""
function dimension(p::Poset; witness=false)
    linearorder(p) && return 1
    discrete(p) && return 2
    lex = linearextensions(p)
    #@show lex
    #println(length(lex))
    for n in 2:length(lex)
        #println(n)
        for tlex in power(lex, n)
            #@show tlex
            #println(collect(tlex))
            q = intersect(collect(tlex))
            if p == q
                witness && println("Witness: ", [string(x) for x in tlex])
                return n
            end
        end
    end
end

