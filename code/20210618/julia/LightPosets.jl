using LightGraphs
import LightGraphs.transitiveclosure!
using Combinatorics

"""

Library for operations on posets.

Posets are `SimpleDiGraph` objects which are acyclic and transitively
closed.

Posets are *irreflexive*: there are no self-loops.

Enforcement of these properties is up to the caller; but they may be
checked using `isposet()`.

"""

"""Posets are simple digraphs."""
Poset = SimpleDiGraph

"""Create linear order from vector."""
function LinearOrder(lorder::Vector{Int})
    p = Poset(length(lorder))
    for i in 1:length(lorder)-1
        add_edge!(p, lorder[i], lorder[i+1])
    end
    transitiveclosure!(p)
    return p
end

"""Check if digraph is transitive."""
function istransitive(g::SimpleDiGraph)
    return g == transitiveclosure(g)
end

"""Check if digraph g is antisymmetric.
Assumes that g is transitive and irreflexive."""
function isantisymmetric(g::SimpleDiGraph)
    for e in edges(g)
        if has_edge(g, reverse(e))
            return false
        end
    end
    return true
end

"""Check that p is a poset."""
function isposet(p::SimpleDiGraph)
    return ! has_self_loops(p) && istransitive(p) && isantisymmetric(p)
end

"""Check that v is minimal in p."""
function isminimal(p::Poset, v::Int)
    return length(inneighbors(p, v)) == 0
end

"""Check that v is maximal in p."""
function ismaximal(p::Poset, v::Int)
    return length(outneighbors(p, v)) == 0
end

"""Return all minimal nodes in p."""
function min(p::Poset)
    return [v for v in vertices(p) if isminimal(p, v)]
end

"""Return all maximal nodes in p."""
function max(p::Poset)
    return [v for v in vertices(p) if ismaximal(p, v)]
end

"""Check if x and y are incomparable in p."""
function concurrent(p::Poset, x::Int, y::Int)
    return (!has_edge(p, x, y)) && (!has_edge(p, y, x))
end

"""Check whether x is above some y in ys."""
function isabove(p::Poset, x::Int, ys::Array{Int})
    for y in ys
        if has_edge(p, (x,y))
            return true
        end
    end
    return false
end

"""Check whether x is below some y in ys."""
function isbelow(p::Poset, x::Int, ys::Array{Int})
    for y in ys
        if has_edge(p, (y,x))
            return true
        end
    end
    return false
end

"""Compute minimum filtration of poset.
Return: array which maps nodes to filtration level"""
function minfiltration(p::Poset)
    rest = [v for v in vertices(p)]
    filter = Array{Int}(undef, nv(p))
    level = 1
    while rest != []
        removeids = Array{Int}(undef, 0)
        for i in eachindex(rest)
            v = rest[i]
            minimal = true
            for x in inneighbors(p, v)
                if x in rest
                    #found incoming edge in rest
                    minimal = false
                    break
                end
            end
            if minimal
                filter[v] = level
                push!(removeids, i)
            end
        end
        #removeids is sorted. use searchsortedfirst?
        rest = [rest[i] for i in eachindex(rest) if ! (i in removeids)]
        level += 1
    end
    return filter
end

"""Compute maximum filtration of poset.
Return: array which maps nodes to filtration level"""
function maxfiltration(p::Poset)
    rest = [v for v in vertices(p)]
    filter = Array{Int}(undef, nv(p))
    level = 1
    while rest != []
        removeids = Array{Int}(undef, 0)
        for i in eachindex(rest)
            v = rest[i]
            maximal = true
            for x in outneighbors(p, v)
                if x in rest
                    #found outgoing edge in rest
                    maximal = false
                    break
                end
            end
            if maximal
                filter[v] = level
                push!(removeids, i)
            end
        end
        #removeids is sorted. use searchsortedfirst?
        rest = [rest[i] for i in eachindex(rest) if ! (i in removeids)]
        level += 1
    end
    return filter
end

"""Compute transitive reduction (Hasse diagram) of poset.
Result is not transitive, hence not a poset
"""
function transitivereduction!(p::Poset)
    redundant = Array{LightGraphs.SimpleGraphs.SimpleEdge}(undef, 0)
    for e in edges(p)
        for pos_intermediary in outneighbors(p, src(e))
            if pos_intermediary == dst(e)
                continue
            end
            if dst(e) in outneighbors(p, pos_intermediary)
                push!(redundant, e)
                break
            end
        end
    end
    for e in redundant
        rem_edge!(p, e)
    end
end
function transitivereduction(p::Poset)::SimpleDiGraph
    q = copy(p)
    transitivereduction!(q)
    return q
end
function hassediagram(p::Poset)::SimpleDiGraph
    return transitivereduction(p)
end

"""Check whether array of nodes is an antichain in p."""
function isantichain(p::Poset, nodes::Array{Int})
    return ne(induced_subgraph(p, nodes)[1]) == 0
end

"""Check whether array of nodes is a maximal antichain in p."""
function ismaximalantichain(p::Poset, nodes::Array{Int})
    if ! isantichain(p, nodes)
        return false
    end
    for v in vertices(p)
        if ! (v in nodes)
            if isantichain(p, vcat(nodes,[v]))
                return false
            end
        end
    end
    return true
end

"""Return array of all maximal antichains.
Attention may be slow: iterates over all subsets of nodes
TODO: Implement Bron-Kerbosch
"""
function maximalantichains(p::Poset)
    res = Array{Array{Int}}(undef, 0)
    for subs in powerset(vertices(p))
        if ismaximalantichain(p, subs)
            push!(res, subs)
        end
    end
    return res
end

"""Check whether array of nodes is a simple path in g."""
function issimplepath(g::AbstractGraph, nodes::Array{Int})
    for o in permutations(nodes)
        good = true
        for i in 1:length(nodes)-1
            if !has_edge(g, o[i], o[i+1])
                good = false
                break
            end
        end
        if good
            return true #return (true, o)
        end
    end
    return false
end

"""Compute height of g: the length of a longest simple path."""
function height(g::AbstractGraph)
    nv(g) <= 1 && return 0 #not enough nodes...
    n = 2
    while n <= nv(g)
        foundone = false #found simple path of length n-1
        for c in combinations(vertices(g), n)
            if issimplepath(g, c)
                foundone = true
                break
            end
        end
        !foundone && break
        n += 1
    end
    return n-2
end

"""Compute "weak height" of poset: the height of its symmetrization."""
function weakheight(p::Poset)
    return height(SimpleGraph(p))
end

"""Check if p is step sequence."""
function stepsequence(p::Poset; witness=false)
    #check if concurrency relation is transitive
    for x in vertices(p)
        for y in vertices(p)
            if concurrent(p, x, y)
                for z in vertices(p)
                    if z != x
                        if concurrent(p, y, z)
                            if ! concurrent(p, x, z)
                                witness && println("Witness: $x, $y, $z")
                                return false
                            end
                        end
                    end
                end
            end
        end
    end
    return true
end

"""Check if p is interval order."""
function intervalorder(p::Poset; witness=false)
    #check if p contains a 2+2
    for x in vertices(p)
        for y in vertices(p)
            if has_edge(p, x, y)
                for z in vertices(p)
                    if z != x
                        for w in vertices(p)
                            if (w != y) && has_edge(p, z, w)
                                if (!has_edge(p, x, w)) && (!has_edge(p, z, y))
                                    witness && println("Witness: $x, $y, $z, $w")
                                    return false
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return true
end

"""Check if p is linear order."""
function linearorder(p::Poset; witness=false)
    #check if concurrency relation is empty
    for x in vertices(p)
        for y in vertices(p)
            if (x != y) && concurrent(p, x, y)
                witness && println("Witness: $x, $y")
                return false
            end
        end
    end
    return true
end

"""Return all extensions of p which satisfy property.
Result is *not* reduced up to iso
"""
#TODO does not work for nv(p)>5 ??
function extensions(p::Poset, property::Function)::Vector{Poset}
    if property(p)
        res = [p]
    else
        res = []
    end
    for x in vertices(p)
        for y in vertices(p)
            if (x != y) && (!has_edge(p, x, y))
                q = copy(p)
                add_edge!(q, x, y)
                transitiveclosure!(q)
                if isantisymmetric(q)
                    #println(string(q))
                    append!(res, extensions(q, property))
                end
            end
        end
    end
    return res
end
# function orderextensions(p)
#     return extensions(p, x->true)
# end
# function stepextensions(p)
#     return extensions(p, x->stepsequence(x))
# end
# function linearextensions(p)
#     return extensions(p, x->linearorder(x))
# end
# function intervalextensions(p)
#     return extensions(p, x->intervalorder(x))
# end

#https://networkx.org/documentation/stable/_modules/networkx/algorithms/dag.html#all_topological_sorts
function linearextensions(G::Poset)
    res = Vector{Vector{Int}}(undef, 0)
    count = copy(indegree(G))
    D = copy(min(G))
    bases = []
    current_sort = []
    while true
        # assert all([count[v] == 0 for v in D])
        if length(current_sort) == nv(G)
            #Push current sort and clean up
            push!(res, current_sort)
            while length(current_sort) > 0
                #assert len(bases) == len(current_sort)
                q = pop!(current_sort)
                for j in outneighbors(G, q)
                    count[j] += 1
                end
                while (length(D) > 0) && (count[D[end]] > 0)
                    pop!(D)
                end
                pushfirst!(D, q)
                if D[end] == bases[end]
                    pop!(bases)
                else
                    break
                end
            end
        else
            length(D) == 0 && error("Input contains cycle!?")
            #Choose and remove a minimal node
            q = pop!(D)
            for j in outneighbors(G, q)
                count[j] -= 1
                #assert count[j] >= 0
                count[j] == 0 && push!(D, j)
            end
            push!(current_sort, q)
            length(bases) < length(current_sort) && push!(bases, q)
        end
        length(bases) == 0 && break
    end
    return (LinearOrder(x) for x in res)
end

function intersect(gs::Vector{T}) where T <: AbstractGraph
    length(gs) == 0 && return T(0)
    length(gs) == 1 && return gs[1]
    return intersect(gs[1], intersect(gs[2:end]))
end

