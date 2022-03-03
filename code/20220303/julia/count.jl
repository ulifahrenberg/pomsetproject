#using IterTools
#using Base.Iterators
include("analysis.jl")
include("generate.jl")
#include("iposets.jl")
using InteractiveUtils
#using Combinatorics
using DataStructures

const n = parse(Int, ARGS[1])

"""Count number of bijections which are checked for isomorphism."""
function c_isomorphic(g::Poset, gvprof::Vprof, 
                    s::Poset, svprof::Vprof)
    #Start with the easy stuff
    if g == s && gvprof == svprof
        return 0
    end
    n = nv(g)
    if n != nv(s) || ne(g) != ne(s)
        return 0
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
            return 0
        end
    end
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
    count = 0
    for pos_isom in Iterators.product(targets...)
        if bijective(pos_isom)
            count += 1
        end
    end
    return count
end

"""Count number of bijections which are checked for isomorphism."""
function countperms(filename; hex=true)
    ps = Array{Tuple{Poset, Vprof}}(undef, 0)
    for ip in eachiposet(filename, hex=hex)
        p = ip.poset
        vprof = initvprof()
        for v in 1:nv(p)
            vprof[v] = (inhash(p, v), outhash(p, v))
        end
        push!(ps, (p, vprof))
    end
    counts = DefaultDict(0)
    for (p1, vprof1) in ps
        for (p2, vprof2) in ps
            nperms = c_isomorphic(p1, vprof1, p2, vprof2)
            counts[nperms] += 1
        end
    end
    return counts
end

#println(countperms(string("data/posets",n,".txt"), hex=false))
