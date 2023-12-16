using Base.Iterators
include("iposets.jl")

"""
Functions for generating posets and iposets of various types
"""

"""
Functions for loading caches:
loadiposetsstrat!
initiposets
"""

"""Load array of iposets from file, stratified by number of edges,
numbers of interfaces, and number of nodes.  Also generate and store
vprofs and keep track of `filled` array
"""
function loadiposetsstrat!(res::Array{Array{Tuple{Iposet, Vprof}}}, filled::Array{Bool}, filename::String; hex=true)
    for ip in eachiposet(filename, hex=hex)
        np = nv(ip.poset)
        ipe = ne(ip.poset)
        kp = length(ip.s)
        lp = length(ip.t)
        vprof = initvprof()
        @inbounds for v in 1:np
            vprof[v] = (inhash(ip, v), outhash(ip, v))
        end
        push!(res[ipe+1, kp+1, lp+1, np], (ip, vprof))
        #If we see one, we see all
        filled[kp + 1, lp + 1, np] = true
    end
end

"""Initialization for closure computations.  Returns alliposets,
filled, and locks arrays for use in gpclosure(n, k, l, alliposets,
filled, locks). Uses files in data/ if fromcache
"""
function initiposets(n::Int, seed::String, fromcache=true;
                     cacheprefix="data/", cachesuffix=".txt",
                     hex=true, quiet=false)
    alliposets = Array{Array{Tuple{Iposet, Vprof}}}(undef, (n*(n-1))÷2 + 1, n+1, n+1, n)
    for i in eachindex(alliposets)
        alliposets[i] = []
    end
    filled = zeros(Bool, n+1, n+1, n)
    cacheprefix *= seed
    #Always load seed1
    fname = string(cacheprefix,"1",cachesuffix)
    if !isfile(fname)
        error("Did not find $fname")
        return
    end
    !quiet && println("Found $fname")
    loadiposetsstrat!(alliposets, filled, fname, hex=hex)
    if fromcache
        !quiet && println("Loading iposets from cache...")
        for k in 2:n
            fname = string(cacheprefix,k,cachesuffix)
            if !isfile(fname)
                !quiet && println("Did not find $fname. Continuing")
                continue
            end
            !quiet && println("Found $fname")
            loadiposetsstrat!(alliposets, filled, fname, hex=hex)
        end
        !quiet && println("Done loading.")
    else
        !quiet && println("Not using cache.")
    end
    locks = Array{ReentrantLock}(undef, (n*(n-1))÷2 + 2, n+1, n+1, n)
    for i in eachindex(locks)
        locks[i] = ReentrantLock()
    end
    return alliposets, filled, locks
end

"""
Functions for closure computations:
gluingclosure
paraclosure
gpclosure
gluings
paras
gpwinkowskis

The arrays alliposets and filled below are for memoisation.
alliposets: num_edg+1 x num_sourc+1 x num_targ+1 x num_points
locks: same
filled: num_sourc+1 x num_targ+1 x num_points
"""

"""Recursively compute gluing closure."""
function gluingclosure(n::Int, k::Int, l::Int, alliposets, filled, locks)
    #@printf("Fill ratio: %f\n", count(x -> x, filled)/length(filled))
    #Return memoized if exists
    if filled[k+1, l+1, n]
        return (x[1] for x in vcat(alliposets[:, k+1, l+1, n]...))
    end
    #If opposites exist, return opposites
    if filled[l+1, k+1, n]
        # lock(locks[end, k+1, l+1, n])
        # for ne in 1:size(alliposets, 1)
        #     for tip in alliposets[ne, l+1, k+1, n]
        #         push!(alliposets[ne, k+1, l+1, n], reverse(tip))
        #     end
        # end
        # filled[k+1, l+1, n] = true
        # unlock(locks[end, k+1, l+1, n])
        # return (x[1] for x in vcat(alliposets[:, k+1, l+1, n]...))
        #Do NOT store opposites
        return (reverse(x[1]) for x in vcat(alliposets[:, l+1, k+1, n]...))
    end
    #Otherwise, generate recursively
    #println("Generating gpi-", n, "-", k, "-", l, " recursively")
    lock(locks[end, k+1, l+1, n])
    Threads.@threads for n1 in 1:(n-1)
        Threads.@threads for n2 in 1:(n-1)
            u = n1 + n2 - n
            if u < 0 || u > n1 || u > n2 || k > n1 || l > n2
                continue
            end
            for ip1 in gluingclosure(n1, k, u, alliposets, filled, locks)
                for ip2 in gluingclosure(n2, u, l, alliposets, filled, locks)
                    ip = glue(ip1, ip2)
                    ipe = ne(ip)
                    pushuptoiso!(alliposets[ipe+1, k+1, l+1, n], ip, locks[ipe+1, k+1, l+1, n])
                end
            end
        end
    end
    filled[k+1, l+1, n] = true
    unlock(locks[end, k+1, l+1, n])
    return (x[1] for x in vcat(alliposets[:, k+1, l+1, n]...))
end

"""Recursively compute parallel closure."""
function paraclosure(n::Int, k::Int, l::Int, alliposets, filled, locks)
    #Return memoized if exists
    if filled[k+1, l+1, n]
        return (x[1] for x in vcat(alliposets[:, k+1, l+1, n]...))
    end
    #This is not symmetric
    # if filled[l+1, k+1, n]
    #     return (reverse(x[1]) for x in vcat(alliposets[:, l+1, k+1, n]...))
    # end
    #Otherwise, generate recursively
    lock(locks[end, k+1, l+1, n])
    Threads.@threads for n1 in 1:(n-1)
        n2 = n - n1
        Threads.@threads for ip1s in 0:k
            ip2s = k - ip1s
            Threads.@threads for ip1t in 0:l
                ip2t = l - ip1t
                if n1 < ip1s || n1 < ip1t || n2 < ip2s || n2 < ip2t
                    continue
                end
                for ip1 in paraclosure(n1, ip1s, ip1t, alliposets, filled, locks)
                    for ip2 in paraclosure(n2, ip2s, ip2t, alliposets, filled, locks)
                        ip = parallel(ip1, ip2)
                        ipe = ne(ip)
                        pushuptoiso!(alliposets[ipe+1, k+1, l+1, n], ip, locks[ipe+1, k+1, l+1, n])
                    end
                end
            end
        end
    end
    filled[k+1, l+1, n] = true
    unlock(locks[end, k+1, l+1, n])
    return (x[1] for x in vcat(alliposets[:, k+1, l+1, n]...))
end

"""Recursively compute gluing-parallel closure."""
function gpclosure(n::Int, k::Int, l::Int, alliposets, filled, locks)
    #Return memoized if exists
    if filled[k+1, l+1, n]
        return (x[1] for x in vcat(alliposets[:, k+1, l+1, n]...))
    end
    #If opposites exist, return opposites
    if filled[l+1, k+1, n]
        return (reverse(x[1]) for x in vcat(alliposets[:, l+1, k+1, n]...))
    end
    #Otherwise, generate recursively
    #println("Generating gpi-", n, "-", k, "-", l, " recursively")
    lock(locks[end, k+1, l+1, n])
    #First, the gluings
    Threads.@threads for n1 in 1:(n-1)
        Threads.@threads for n2 in 1:(n-1)
            u = n1 + n2 - n
            if u < 0 || u > n1 || u > n2 || k > n1 || l > n2
                continue
            end
            for ip1 in gpclosure(n1, k, u, alliposets, filled, locks)
                for ip2 in gpclosure(n2, u, l, alliposets, filled, locks)
                    ip = glue(ip1, ip2)
                    ipe = ne(ip)
                    pushuptoiso!(alliposets[ipe+1, k+1, l+1, n], ip, locks[ipe+1, k+1, l+1, n])
                end
            end
        end
    end
    #Now, the parallel compositions
    Threads.@threads for n1 in 1:(n-1)
        n2 = n - n1
        Threads.@threads for ip1s in 0:k
            ip2s = k - ip1s
            Threads.@threads for ip1t in 0:l
                ip2t = l - ip1t
                if n1 < ip1s || n1 < ip1t || n2 < ip2s || n2 < ip2t
                    continue
                end
                for ip1 in gpclosure(n1, ip1s, ip1t, alliposets, filled, locks)
                    for ip2 in gpclosure(n2, ip2s, ip2t, alliposets, filled, locks)
                        ip = parallel(ip1, ip2)
                        ipe = ne(ip)
                        pushuptoiso!(alliposets[ipe+1, k+1, l+1, n], ip, locks[ipe+1, k+1, l+1, n])
                    end
                end
            end
        end
    end
    filled[k+1, l+1, n] = true
    unlock(locks[end, k+1, l+1, n])
    return (x[1] for x in vcat(alliposets[:, k+1, l+1, n]...))
end

"""Generate all gluings of left and right and put into ips, up to iso."""
function gluings(n::Int, k::Int, l::Int, left, right, ips)
    Threads.@threads for n1 in 1:n
        Threads.@threads for n2 in 1:n
            u = n1 + n2 - n
            if u < 0 || u > n1 || u > n2 || k > n1 || l > n2
                continue
            end
            for ip1 in (tip[1] for tip in vcat(left[:, k+1, u+1, n1]...))
                for ip2 in (tip[1] for tip in vcat(right[:, u+1, l+1, n2]...))
                    ip = glue(ip1, ip2)
                    ipe = ne(ip)
                    pushuptoiso!(ips[ipe+1, k+1, l+1, n], ip)
                end
            end
        end
    end
    return (x[1] for x in vcat(ips[:, k+1, l+1, n]...))
end

"""Generate all parallel products of left and right and put into ips,
up to iso.
"""
function paras(n::Int, k::Int, l::Int, left, right, ips)
    Threads.@threads for n1 in 1:(n-1)
        n2 = n - n1
        Threads.@threads for ip1s in 0:k
            ip2s = k - ip1s
            Threads.@threads for ip1t in 0:l
                ip2t = l - ip1t
                if n1 < ip1s || n1 < ip1t || n2 < ip2s || n2 < ip2t
                    continue
                end
                for ip1 in (tip[1] for tip in vcat(left[:, ip1s+1, ip1t+1, n1]...))
                    for ip2 in (tip[1] for tip in vcat(right[:, ip2s+1, ip2t+1, n2]...))
                        ip = parallel(ip1, ip2)
                        ipe = ne(ip)
                        pushuptoiso!(ips[ipe+1, k+1, l+1, n], ip)
                    end
                end
            end
        end
    end
    return (x[1] for x in vcat(ips[:, k+1, l+1, n]...))
end

"""Generate gp-winkowskis with profile n, k, l.
Needs starters and terminators up to n nodes
"""
function gpwinkowskis(n::Int, k::Int, l::Int, gpstar, gpterm, gpwink, filled, locks)
    #Return memoized if exists
    if filled[k+1, l+1, n]
        return (x[1] for x in vcat(gpwink[:, k+1, l+1, n]...))
    end
    #If opposites exist, return opposites
    if filled[l+1, k+1, n]
        return (reverse(x[1]) for x in vcat(gpwink[:, l+1, k+1, n]...))
    end
    #Otherwise, generate recursively
    lock(locks[end, k+1, l+1, n])
    #First, the gluings
    Threads.@threads for n1 in 1:(n-1)
        Threads.@threads for n2 in 1:(n-1)
            u = n1 + n2 - n
            if u < 0 || u > n1 || u > n2 || k > n1 || l > n2
                continue
            end
            for ip1 in flatten((
                gpwinkowskis(n1, k, u, gpstar, gpterm, gpwink, filled, locks),
                (tip[1] for tip in vcat(gpterm[:, k+1, u+1, n1]...))))
                for ip2 in flatten((
                    gpwinkowskis(n2, u, l, gpstar, gpterm, gpwink, filled, locks),
                    (tip[1] for tip in vcat(gpstar[:, u+1, l+1, n2]...))))
                    ip = glue(ip1, ip2)
                    ipe = ne(ip)
                    #lock(locks[ipe+1, k+1, l+1, n])
                    pushuptoiso!(gpwink[ipe+1, k+1, l+1, n], ip, locks[ipe+1, k+1, l+1, n])
                    #unlock(locks[ipe+1, k+1, l+1, n])
                end
            end
        end
    end
    #Now, the parallel compositions
    Threads.@threads for n1 in 1:(n-1)
        n2 = n - n1
        Threads.@threads for ip1s in 0:k
            ip2s = k - ip1s
            Threads.@threads for ip1t in 0:l
                ip2t = l - ip1t
                if n1 < ip1s || n1 < ip1t || n2 < ip2s || n2 < ip2t
                    continue
                end
                for ip1 in gpwinkowskis(n1, ip1s, ip1t, gpstar, gpterm, gpwink, filled, locks)
                    for ip2 in gpwinkowskis(n2, ip2s, ip2t, gpstar, gpterm, gpwink, filled, locks)
                        ip = parallel(ip1, ip2)
                        ipe = ne(ip)
                        pushuptoiso!(gpwink[ipe+1, k+1, l+1, n], ip, locks[ipe+1, k+1, l+1, n])
                    end
                end
            end
        end
    end
    filled[k+1, l+1, n] = true
    unlock(locks[end, k+1, l+1, n])
    return (x[1] for x in vcat(gpwink[:, k+1, l+1, n]...))
end

"""
Functions for generating posets and iposets:
posets
iposets
gpposets
gpiposets
gpintorders
spintorders
spposets
gpstarters
gpterminators
gpdiscrete
gpwinkowskis
gpleftwinkowskis
gpifromlw
"""

"""Generate all posets of size n. Return as array
timing (without cache):
6: .3s
7: 32s
8: 220min (?)
"""
#TODO THREADING??
function posets(n; fromcache=true)
    if fromcache
        return (x.poset for x in loadiposets(string("data/posets",n,".txt"), hex=false))
    end
    n == 0 && return [Poset(0)]
    n == 1 && return [Poset(1)]
    max_edges = (n*(n-1))÷2
    all_edges = Array{Tuple{Int, Int}}(undef, max_edges)
    en = 1
    #Generate all possible edges. Guarantee irreflexive and antisymmetric
    @inbounds for i in 1:(n-1)
        for j in (i+1):n
            all_edges[en] = (i, j)
            en += 1
        end
    end
    #println(all_edges)
    #res is graded by number of edges
    res = Array{Array{Tuple{Poset, Vprof}}}(undef, max_edges + 1)
    @inbounds for i in eachindex(res)
        res[i] = []
    end
    @inbounds for numedges in 0:max_edges
        #Get all the ways to distribute the edges
        for distribution in combinations(all_edges, numedges)
            g = Poset(n)
            for e in distribution
                add_edge!(g, e[1], e[2])
            end
            transitiveclosure!(g, false)
            #Now g is transitive, irreflexive and antisymmetric
            pne = ne(g)
            pushuptoiso!(res[pne+1], g)
        end
    end
    return (x[1] for x in vcat(res...))
end

"""Generate all iposets with sources k, targets l, and underlying poset p."""
function iposetsone(k::Int, l::Int, p::Poset)
    mins = min(p)
    maxes = max(p)
    if discrete(p) #If p is discrete, no need to permute sources
        pos_sources = [[mins[i] for i in 1:k]] #assert length(mins) < k
    else
        pos_sources = permutations(mins, k)
    end        
    pos_targets = permutations(maxes, l)
    ips = Array{Tuple{Iposet, Vprof}}(undef, 0)
    for combination in Iterators.product(pos_sources, pos_targets)
        nip = Iposet(Tuple(combination[1]), Tuple(combination[2]), p)
        pushuptoiso!(ips, nip)
    end
    return ips
end

function iposets(k::Int, l::Int, posets)
    iposets = Array{Array{Tuple{Iposet, Vprof}}}(undef, 0)
    for p in posets
        push!(iposets, iposetsone(k, l, p))
    end
    return (x[1] for x in vcat(iposets...))
end

function iposets(n::Int, posets)
    elems = Array{Array{Iposet}}(undef, n+1, n+1)
    for i in eachindex(elems)
        elems[i] = []
    end
    #First create for all k<=l
    Threads.@threads for k in 0:n
        Threads.@threads for l in k:n
            elems[k+1, l+1] = collect(iposets(k, l, posets))
        end
    end
    #Now the rest, using opposition
    Threads.@threads for k in 1:n
        Threads.@threads for l in 0:(k-1)
            elems[k+1, l+1] = [reverse(x) for x in elems[l+1, k+1]]
        end
    end
    return vcat(elems...)
end

"""Return all iposets with n nodes, k sources and l targets."""
function iposets(n::Int, k::Int, l::Int; fromcache=true)
    allposets = collect(posets(n, fromcache=fromcache))
    return iposets(k, l, allposets)
end

"""Return all iposets with n nodes.
timings (with posets cache):
4:    .001s
5:    .02s
6:   2.1s
7:  495s = 8min15s
"""
function iposets(n::Int; fromcache=true)
    allposets = collect(posets(n, fromcache=fromcache))
    return iposets(n, allposets)
end

"""Return all gp-posets with n nodes.
timings (without cache):
7:  .7s
8: 43s
"""
function gpposets(n, fromcache=true)
    return gpiposets(n, 0, 0, fromcache)
end

"""Return all gp-iposets with n nodes.
timings (without cache):
5:    .06s
6:    2.1s
7: 3min50s
"""
function gpiposets(n, fromcache=true)
    alliposets, filled, locks = initiposets(n, "gpi", fromcache)
    elems = Array{Array{Iposet}}(undef, n+1, n+1)
    for i in eachindex(elems)
        elems[i] = []
    end
    #First create for all k<=l
    Threads.@threads for k in 1:(n+1)
        Threads.@threads for l in k:(n+1)
            elems[k, l] = collect(gpclosure(n, k-1, l-1, alliposets, filled, locks))
        end
    end
    #Now the rest, using opposition
    Threads.@threads for k in 2:(n+1)
        Threads.@threads for l in 1:(k-1)
            elems[k, l] = [reverse(x) for x in elems[l, k]]
        end
    end
    return vcat(elems...)
end

"""Return all gp-iposets with n nodes, k sources, and l targets."""
function gpiposets(n, k, l, fromcache=true)
    alliposets, filled, locks = initiposets(n, "gpi", fromcache)
    return gpclosure(n, k, l, alliposets, filled, locks)
end

"""Return all gpi-interval orders with n nodes, k sources, and l targets."""
function gpintorders(n, k, l)
    alliposets,_, locks = initiposets(n, "gpdisc", quiet=true)
    filled = zeros(Bool, n+1, n+1, n)
    return gluingclosure(n, k, l, alliposets, filled, locks)
end

"""Return all interval orders with n nodes. Uses gpdiscrete cache.
timings:
5:   .005s
6:   .03s
7:   .7s
8:   18s
9:  13min40s = 820s
10: 20h ?
"""
function intorders(n)
    return gpintorders(n,0,0)
end

"""Return all sp-interval orders with n nodes.
timings (using gpintorders cache):
5:   .002s
6:   .006s
7:   .04s
8:   .36s
9:   5.4s
"""
function spintorders(n)
    alliposets,_, locks = initiposets(n, "gpint", cachesuffix="-0-0.txt")
    filled = zeros(Bool, n+1, n+1, n)
    return gpclosure(n, 0, 0, alliposets, filled, locks)
end

"""Return all sp-posets with n nodes.
timings (without cache):
 8:  .15s
 9:  1.3s
10: 22s
11:  ??
"""
function spposets(n, fromcache=true)
    alliposets, filled, locks = initiposets(n, "sp", fromcache)
    return gpclosure(n, 0, 0, alliposets, filled, locks)
end

"""Generate gp-starters on n points."""
function gpstarters(n, fromcache=true)
    alliposets, filled, locks = initiposets(n, "gpstar", fromcache)
    elems = Array{Array{Iposet}}(undef, n+1, n+1)
    for i in eachindex(elems)
        elems[i] = []
    end
    Threads.@threads for k in 1:(n+1)
        Threads.@threads for l in 1:(n+1)
            elems[k, l] = collect(paraclosure(n, k-1, l-1, alliposets, filled, locks))
        end
    end
    return vcat(elems...)
end
    
"""Generate gp-terminators on n points."""
function gpterminators(n, fromcache=true)
    alliposets, filled, locks = initiposets(n, "gpterm", fromcache)
    elems = Array{Array{Iposet}}(undef, n+1, n+1)
    for i in eachindex(elems)
        elems[i] = []
    end
    Threads.@threads for k in 1:(n+1)
        Threads.@threads for l in 1:(n+1)
            elems[k, l] = collect(paraclosure(n, k-1, l-1, alliposets, filled, locks))
        end
    end
    return vcat(elems...)
end

"""Generate all discrete gp-iposets on n points."""
function gpdiscrete(n)
    seed = singletons()
    res = Array{Tuple{Iposet, Vprof}}(undef, 0)
    for t in product((seed for _ in 1:n)...)
        pushuptoiso!(res, parallel(collect(t)))
    end
    return (x[1] for x in res)
end

"""Generate gp-winkowskis with n points. Use data/ cache.
Needs gpstar and gpterm cache up to n-1 (these are not generated if missing!)
timing: 5   .004s
        6   .1s
        7  5s
        8  12min (?)
"""
function gpwinkowskis(n, fromcache=true)
    gpstar,_,_ = initiposets(n, "gpstar", quiet=true)
    gpterm,_,_ = initiposets(n, "gpterm", quiet=true)
    gpwink, filled, locks = initiposets(n, "gpwink", fromcache)
    elems = Array{Array{Iposet}}(undef, n+1, n+1)
    for i in eachindex(elems)
        elems[i] = []
    end
    #First create for all k<=l
    Threads.@threads for k in 1:(n+1)
        Threads.@threads for l in k:(n+1)
            elems[k, l] = collect(gpwinkowskis(n, k-1, l-1, gpstar, gpterm, gpwink, filled, locks))
            #saveiposets(elems[k, l], string("gpwink",n,"-",k-1,"-",l-1,".txt"))
        end
    end
    #Now the rest, using opposition
    Threads.@threads for k in 2:(n+1)
        Threads.@threads for l in 1:(k-1)
            elems[k, l] = [reverse(x) for x in elems[l, k]]
        end
    end
    return vcat(elems...)
end

"""Generate gp-leftwinkowskis by gluing winkowskis with terminators.
Needs winkowskis with n nodes in gpwink and gplw, and terminators up
to n nodes in gpterm.
timings:
5     .04s
6     .3s
7     22s
8   113min
"""
function gpleftwinkowskis(n)
    gpterm,_,_ = initiposets(n, "gpterm", quiet=true)
    gpwink,_,_ = initiposets(n, "gpwink", quiet=true)
    gplw = copy(gpwink)
    elems = Array{Array{Iposet}}(undef, n+1, n+1)
    for i in eachindex(elems)
        elems[i] = []
    end
    Threads.@threads for k in 1:(n+1)
        Threads.@threads for l in 1:(n+1)
            elems[k, l] = collect(gluings(n, k-1, l-1, gpwink, gpterm, gplw))
            #saveiposets(elems[k, l], string("gplw",n,"-",k-1,"-",l-1,".txt"))
        end
    end
    return vcat(elems...)
end

"""Generate gpi-posets by gluing starters with gp-leftwinkowskis.
Needs leftwinkowskis with n nodes in gplwink and gpi, and starters up
to n nodes in gpstar
timings:
5   .05s
6   1.5s
7   4min45s
"""
function gpifromlw(n)
    gpstar,_,_ = initiposets(n, "gpstar", quiet=true)
    gplwink,_,_ = initiposets(n, "gplw", quiet=true)
    gpi = copy(gplwink)
    elems = Array{Array{Iposet}}(undef, n+1, n+1)
    for i in eachindex(elems)
        elems[i] = []
    end
    Threads.@threads for k in 1:(n+1)
        Threads.@threads for l in k:(n+1)
            elems[k, l] = collect(gluings(n, k-1, l-1, gpstar, gplwink, gpi))
            #saveiposets(elems[k, l], string("gpifromlw",n,"-",k-1,"-",l-1,".txt"))
        end
    end
    #now the opposites
    for k in 2:(n+1)
        for l in 1:(k-1)
            elems[k, l] = [reverse(x) for x in elems[l, k]]
        end
    end
    return vcat(elems...)
end

"""
Other stuff, old and untested.
"""

"""Generate all connected gluings, without interfaces, of leftK and
rightK for K<n.  Split output into files gluingsn-X-Y.txt which are
subsets of gpn obtained by gluing leftX with rightY.  Each
gluingsn-X-Y.txt is reduced up to iso, but elements of different files
may be isomorphic
"""
#UNTESTED
function connectedgluings(n; leftfileprefix="left", rightfileprefix="right")
    k, l = 0, 0
    lefts,_,_ = initiposets(n-1, leftfileprefix)
    rights,_,_ = initiposets(n-1, rightfileprefix)
    for n1 in 1:(n-1)
        for n2 in n1:(n-1)
            u = n1 + n2 - n
            if u < 0 || u > n1 || u > n2 || k > n1 || l > n2
                continue
            end
            res = Array{Tuple{Iposet, Vprof}}(undef, 0)
            fname = string("gluings",n,"-",n1,"-",n2,".txt")
            println("Generating $fname...")
            for ip1 in (tip[1] for tip in vcat(lefts[:, k+1, u+1, n1]...))
                for ip2 in (tip[1] for tip in vcat(rights[:, u+1, l+1, n2]...))
                    ip = glue(ip1, ip2)
                    if !connected(ip)
                        continue
                    end
                    pushuptoiso!(res, ip)
                end
            end
            println("Saving $fname...")
            saveiposets(res, fname)
        end
    end
end

"""Generate gluingsn-l-k from gluingsn-k-l."""
#UNTESTED
function mirror(n, k, l; fileprefix="gluings")
    infile = string("data/",fileprefix,n,"-",k,"-",l,".txt")
    outfile = string(fileprefix,n,"-",l,"-",k,".txt")
    ips = loadiposets(infile)
    for ip in ips
        reverse!(ip)
    end
    saveiposets(ips, outfile)
end


