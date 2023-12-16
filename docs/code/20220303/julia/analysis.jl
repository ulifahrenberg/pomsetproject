include("iposets.jl")

"""Compute difference a-b up to isomorphism: return all iposets in a
which do not have an isomorph in b.
Assumes that elements of a are pairwise non-isomorphic
"""
function diffuptoiso(a::Array{Iposet}, b::Array{Iposet})
    #Prepare vprofs for b
    bvprofs = Array{Vprof}(undef, length(b))
    @inbounds for i in 1:length(b)
        bvprofs[i] = initvprof()
        for v in vertices(b[i])
            bvprofs[i][v] = (inhash(b[i], v), outhash(b[i], v))
        end
    end
    #Now check every iposet in a
    res = Array{Iposet}(undef, 0)
    for p in a
        #println("a: ", string(p))
        vprof = initvprof()
        for v in vertices(p)
            vprof[v] = (inhash(p, v), outhash(p, v))
        end
        foundone = false
        @inbounds for i in 1:length(b)
            #println("b: ", string(b[i]))
            if isomorphic(p, vprof, b[i], bvprofs[i])
                #println("Good")
                foundone = true
                #Remove iposets from b which have been found in a
                #Depending on b, this might be faster
                #deleteat!(b, i)
                #deleteat!(bvprofs, i)
                break
            end
        end
        if !foundone
            #println("Not found in b")
            push!(res, p)
        end
    end
    return res
end
function diffuptoiso(a::T, b::T) where T <: Union{Array{Poset}, Array{Poset{Int}}}
    return [x.poset for x in diffuptoiso([Iposet(x) for x in a], [Iposet(x) for x in b])]
end

"""Compute symmetric difference (a-b)+(b-a) up to isomorphism: return
reduceuptoiso(a cup b), but faster.  Assumes that elements of a and of
b are pairwise non-isomorphic
"""
function symdiffuptoiso(a::Array{Iposet}, b::Array{Iposet})
    #Prepare vprofs for b
    bvprofs = Array{Vprof}(undef, length(b))
    @inbounds for i in 1:length(b)
        bvprofs[i] = initvprof()
        for v in vertices(b[i])
            bvprofs[i][v] = (inhash(b[i], v), outhash(b[i], v))
        end
    end
    #Now filter b for every iposet in a
    for p in a
        #println("a: ", string(p))
        vprof = initvprof()
        for v in vertices(p)
            vprof[v] = (inhash(p, v), outhash(p, v))
        end
        @inbounds for i in 1:length(b)
            #println("b: ", string(b[i]))
            if isomorphic(p, vprof, b[i], bvprofs[i])
                #println("Good")
                deleteat!(b, i)
                deleteat!(bvprofs, i)
                break
            end
        end
    end
    return vcat(a, b)
end

"""From an array of iposets, return all associated posets, up to iso."""
function posets(a::Array{Iposet})
    return reduceiso([x.poset for x in a])
end

"""Find all posets on n points which are not gp.
timing: 7  1.2s
        8  29s
        9  30min ?
"""
function posetsnotgp(n; fileprefix="data/gpi", filesuffix="-0-0.txt")
    ps = collect(posets(n))
    filename = string(fileprefix, n, filesuffix)
    gps = [ip.poset for ip in eachiposet(filename)]
    return diffuptoiso(ps, gps)
end

"""Return all p in posets which have no induced subgraph in subs."""
function nosubs(posets, subs)
    res = Array{Poset}(undef, 0)
    for p in posets
        hasnosub = true
        for s in subs
            sg, _ = subgraph(p, s)
            if sg
                hasnosub = false
                break
            end
        end
        hasnosub && push!(res, p)
    end
    return res
end

"""Find forbidden substructures for gp-posets on up to nmax points.
timing: 6  .01s
        7  .3s
        8  24s
        9  4300s
"""
function findforbiddensubs(nmax; fileprefix="data/gpi", filesuffix="-0-0.txt")
    n = 5
    fsubs = Array{Poset}(undef, 0)
    while true
        n += 1
        n > nmax && break
        pngs = posetsnotgp(n, fileprefix=fileprefix, filesuffix=filesuffix)
        newfsubs = nosubs(pngs, fsubs)
        if !isempty(newfsubs)
            println("Found new forbidden substructure(s) on $n points:")
            for s in newfsubs
                println(string(s))
            end
            append!(fsubs, newfsubs)
        end
    end
    return fsubs
end

"""Split count of iposet file up into k-l source-target pairs."""
function splitcount(n, fileprefix="gpi", printname="\\GPI"; dirprefix="data/")
    infilename = string(dirprefix*fileprefix, n, ".txt")
    res = zeros(Int64, n+1, n+1)
    for ip in eachiposet(infilename)
        k = length(ip.s)
        l = length(ip.t)
        res[k+1, l+1] += 1
    end
    print("\\begin{tabular}[t]{r|")
    for l in 0:n
        print("r")
    end
    println("}")
    print("\$",printname,"(",n,")\$\\!\\!")
    for l in 0:n
        print(" & ",l)
    end
    println(" \\\\\\hline")
    for k in 0:n
        print(k)
        for l in 0:(k-1)
            print(" &")
        end
        for l in k:n
            print(" & ",res[k+1, l+1])
        end
        println(" \\\\")
    end
    println("\\end{tabular}")
end

"""Sift iposet file for property."""
function siftfile(property::Function, infilename, outfilename)
    saveiposets([ip for ip in eachiposet(infilename)
                 if property(ip)], outfilename)
end

"""Sift iposet file for winkowskis."""
function siftwink(n, infileprefix="data/gpi", outfileprefix="siftwink")
    infilename = string(infileprefix, n, ".txt")
    outfilename = string(outfileprefix, n, ".txt")
    siftfile(winkowski, infilename, outfilename)
end

"""Sift iposet file for discrete non-starters and right winkowskis
with empty left interface.
"""
function siftleft(n, infileprefix="data/gpi", outfileprefix="siftleft")
    infilename = string(infileprefix, n, ".txt")
    outfilename = string(outfileprefix, n, ".txt")
    property(ip) = length(ip.s)==0 && (discrete(ip) || rightwinkowski(ip)) && !starter(ip)
    siftfile(property, infilename, outfilename)
end

"""Sift iposet file for discrete non-terminators and left winkowskis
with empty right interface.
"""
function siftright(n, infileprefix="data/gpi", outfileprefix="siftright")
    infilename = string(infileprefix, n, ".txt")
    outfilename = string(outfileprefix, n, ".txt")
    property(ip) = length(ip.t)==0 && (discrete(ip) || leftwinkowski(ip)) && !terminator(ip)
    siftfile(property, infilename, outfilename)
end

"""Sift iposet file for (weakly) connected and empty interfaces."""
function siftconn(n, infileprefix="data/gpi", outfileprefix="siftconn")
    infilename = string(infileprefix, n, ".txt")
    outfilename = string(outfileprefix, n, ".txt")
    property(ip) = length(ip.s)==0 && length(ip.t)==0 && connected(ip)
    siftfile(property, infilename, outfilename)
end

"""Sift iposet file for interface consistency."""
function siftintfcons(n, infileprefix="data/iposets", outfileprefix="ici")
    infilename = string(infileprefix, n, ".txt")
    outfilename = string(outfileprefix, n, ".txt")
    siftfile(intfconsistent, infilename, outfilename)
end

"""Sift iposet file for interface consistent Winkowskis."""
function siftintfcwink(n, infileprefix="data/iposets", outfileprefix="icw")
    infilename = string(infileprefix, n, ".txt")
    outfilename = string(outfileprefix, n, ".txt")
    property(ip) = intfconsistent(ip) && winkowski(ip)
    siftfile(property, infilename, outfilename)
end

"""Reduce iposets in file up to iso.
timing (depending on #iposets in file and #iposets in output):
 25k   3k     .8s
 14k  14k     7s
 25k  25k    10s
500k  26k   2min
"""
function reduceuptoiso(filename; hex=true, comma=false)
    seen = Array{Tuple{Iposet, Vprof}}(undef, 0)
    return (x for x in eachiposet(filename, hex=hex, comma=comma)
            if pushuptoiso!(seen, x))
end

"""Try to deglue iposet along a maximal antichain."""
#Only works for iposets without interfaces
function deglue(ip::Iposet)
    if ip.s != () || ip.t != ()
        error("Input has interfaces")
        return
    end
    println("Trying to deglue ", string(ip))
    res = Array{Tuple{Iposet, Iposet}}(undef, 0)
    p = ip.poset
    for a in maximalantichains(p)
        if a == min(p) || a == max(p)
            continue
        end
        #println(a)
        #Find left and right iposets for decomp
        vleft = [v for v in vertices(p) if isabove(p, v, a)]
        vright = [v for v in vertices(p) if isbelow(p, v, a)]
        left, lmap = induced_subgraph(p, vcat(vleft, a))
        right, rmap = induced_subgraph(p, vcat(vright, a))
        #println(vleft," ",vright," ",lmap," ",rmap)
        #Find embedding of a into left and right
        aleft = Array{Int}(undef, length(a))
        for i in eachindex(a)
            for j in eachindex(lmap)
                if lmap[j] == a[i]
                    aleft[i] = j
                end
            end
        end
        aright = Array{Int}(undef, length(a))
        for i in eachindex(a)
            for j in eachindex(rmap)
                if rmap[j] == a[i]
                    aright[i] = j
                end
            end
        end
        #println(aleft," ",aright)
        #Now we have all data
        ip1 = Iposet((), Tuple(aleft), left)
        ip2 = Iposet(Tuple(aright), (), right)
        print("Checking (", string(ip1),") glue (", string(ip2),"): ")
        iq = glue(ip1, ip2)
        if isomorphic(ip, iq)
            println("ok!")
            push!(res, (ip1, ip2))
        else
            println("no")
        end
    end
    length(res) == 0 && println("None found!")
    return res
end
function deglue(p::Poset)
    return deglue(Iposet(p))
end

"""Try to find serial decomposition."""
#Only works for iposets without interfaces
function deserial(ip::Iposet)
    if ip.s != () || ip.t != ()
        error("Input has interfaces")
        return
    end
    println("Trying to de-serialize ", string(ip))
    res = Array{Tuple{Iposet, Iposet}}(undef, 0)
    p = ip.poset
    for (vleft, vright) in partitions(vertices(p), 2)
        #println(left," ",right)
        left, _ = induced_subgraph(p, vleft)
        right, _ = induced_subgraph(p, vright)
        ip1 = Iposet(left)
        ip2 = Iposet(right)
        print("Checking (", string(ip1),") serial (", string(ip2),"): ")
        iq = glue(ip1, ip2)
        if isomorphic(ip, iq)
            println("ok!")
            push!(res, (ip1, ip2))
        else
            println("no")
        end
    end
    length(res) == 0 && println("None found!")
    return res
end
function deserial(p::Poset)
    return deserial(Iposet(p))
end

