include("iposets.jl")

"""
Various closure computations for sets of (i)posets:

gluingclosure
paraclosure
gpclosure
gpsclosure
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

"""Recursively compute gluing-parallel-symmetric closure.
Like gpclosure, but permutes interfaces arbitrarily
"""
function gpsclosure(n::Int, k::Int, l::Int, alliposets, filled, locks)
    #Return memoized if exists
    if filled[k+1, l+1, n]
        return (x[1] for x in vcat(alliposets[:, k+1, l+1, n]...))
    end
    #If opposites exist, return opposites
    if filled[l+1, k+1, n]
        return (reverse(x[1]) for x in vcat(alliposets[:, l+1, k+1, n]...))
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
            for ip1 in gpsclosure(n1, k, u, alliposets, filled, locks)
                for ip2 in gpsclosure(n2, u, l, alliposets, filled, locks)
                    ip = glue(ip1, ip2)
                    ipe = ne(ip)
                    #Now permute interfaces
                    seen = Array{Tuple{Iposet, Vprof}}(undef, 0)
                    for pip in (Iposet(Tuple(ps), Tuple(pt), ip.poset)
                                for ps in permutations(ip.s)
                                for pt in permutations(ip.t)
                                if pushuptoiso!(seen, Iposet(Tuple(ps), Tuple(pt), ip.poset)))
                        pushuptoiso!(alliposets[ipe+1, k+1, l+1, n], pip, locks[ipe+1, k+1, l+1, n])
                    end
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
                for ip1 in gpsclosure(n1, ip1s, ip1t, alliposets, filled, locks)
                    for ip2 in gpsclosure(n2, ip2s, ip2t, alliposets, filled, locks)
                        ip = parallel(ip1, ip2)
                        ipe = ne(ip)
                        #Now permute interfaces
                        seen = Array{Tuple{Iposet, Vprof}}(undef, 0)
                        for pip in (Iposet(Tuple(ps), Tuple(pt), ip.poset)
                                    for ps in permutations(ip.s)
                                    for pt in permutations(ip.t)
                                    if pushuptoiso!(seen, Iposet(Tuple(ps), Tuple(pt), ip.poset)))
                            pushuptoiso!(alliposets[ipe+1, k+1, l+1, n], pip, locks[ipe+1, k+1, l+1, n])
                        end
                    end
                end
            end
        end
    end
    filled[k+1, l+1, n] = true
    unlock(locks[end, k+1, l+1, n])
    return (x[1] for x in vcat(alliposets[:, k+1, l+1, n]...))
end

"""Recursively compute gluing-parallel-symmetric closure.
Use poset-only generation algo
BUG puts some 2-point posets into wrong place
Iposet("2 0") gets put into both 2-2-2 (good) and 2-1-1 (wrong)
"""
function gpsclosure2(n, k, l, posets, filled, locks)
    #Return memoized if exists
    if filled[k+1, l+1, n]
        return (x[1] for x in vcat(posets[:, k+1, l+1, n]...))
    end
    #If opposites exist, return opposites
    if filled[l+1, k+1, n]
        return (reverse(x[1]) for x in vcat(posets[:, l+1, k+1, n]...))
    end
    #Otherwise, generate recursively
    lock(locks[end, k+1, l+1, n])
    #First, the parallel compositions
    #Threads.@threads for n1 in 1:(n-1)
    for n1 in 1:(n-1)
        n2 = n - n1
        #Threads.@threads for ip1s in 1:k
        for ip1s in 1:k
            ip2s = k - ip1s
            #Threads.@threads for ip1t in 1:l
            for ip1t in 1:l
                ip2t = l - ip1t
                if n1 < ip1s || n1 < ip1t || n2 < ip2s || n2 < ip2t
                    continue
                end
                for p1 in gpsclosure2(n1, ip1s, ip1t, posets, filled, locks)
                    for p2 in gpsclosure2(n2, ip2s, ip2t, posets, filled, locks)
                        p = parallel(p1, p2)
                        #println(string(p))
                        pe = ne(p)
                        pushuptoiso!(posets[pe+1, k+1, l+1, n], p, locks[pe+1, k+1, l+1, n])
                    end
                end
            end
        end
    end
    #Second, the Winkowski gluings
    #Threads.@threads for n1 in 1:(n-1)
    for n1 in 1:(n-1)
        #Threads.@threads for n2 in 1:(n-1)
        for n2 in 1:(n-1)
            u = n1 + n2 - n
            if u < 1 || u > n1 || u > n2 || k > n1 || l > n2
                continue
            end
            for p1 in gpsclosure2(n1, k, u, posets, filled, locks)
                for p2 in gpsclosure2(n2, u, l, posets, filled, locks)
                    for p in wglue(p1, p2)
                        #println(string(p))
                        pe = ne(p)
                        pushuptoiso!(posets[pe+1, k+1, l+1, n], p, locks[pe+1, k+1, l+1, n])
                    end
                end
            end
        end
    end
    #Third, the right interior gluings
    #Threads.@threads for n1 in 1:(n-1)
    for n1 in 1:(n-1)
        n2 = n - n1
        # n1 + u = k
        u = k - n1
        #println("$n1, $n2, $k, $u, $l")
        if u < 1 || u > n2 || l > n2
            continue
        end
        p1 = Poset(n1)
        for p2 in gpsclosure2(n2, u, l, posets, filled, locks)
            discrete(p2) && continue
            p = riglue(p1, p2)
            #println(string(p))
            pe = ne(p)
            pushuptoiso!(posets[pe+1, k+1, l+1, n], p, locks[pe+1, k+1, l+1, n])
        end
    end
    #Fourth, the left interior gluings
    #Threads.@threads for n1 in 1:(n-1)
    for n1 in 1:(n-1)
        n2 = n - n1
        # n2 + u = l
        u = l - n2
        #println("$n1, $n2, $k, $u, $l")
        if u < 1 || k > n1 || u > n1
            continue
        end
        p2 = Poset(n2)
        for p1 in gpsclosure2(n1, k, u, posets, filled, locks)
            discrete(p1) && continue
            p = liglue(p1, p2)
            #println(string(p))
            pe = ne(p)
            pushuptoiso!(posets[pe+1, k+1, l+1, n], p, locks[pe+1, k+1, l+1, n])
        end
    end
    #Finally, some serial products
    if k > 0 && l > 0 && k + l == n
        p1 = Poset(k)
        p2 = Poset(l)
        p = glue(Iposet(p1), Iposet(p2)).poset
        pe = ne(p)
        pushuptoiso!(posets[pe+1, k+1, l+1, n], p, locks[pe+1, k+1, l+1, n])
    end
    filled[k+1, l+1, n] = true
    unlock(locks[end, k+1, l+1, n])
    return (x[1] for x in vcat(posets[:, k+1, l+1, n]...))
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
                    pushuptoiso!(gpwink[ipe+1, k+1, l+1, n], ip, locks[ipe+1, k+1, l+1, n])
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

