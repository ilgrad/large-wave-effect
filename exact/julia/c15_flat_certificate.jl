#!/usr/bin/env julia
# No flat point on the orbit closure H_15, hence B_15 < sqrt(15) (paper, Theorem bn-saturation (4)).
# Finds and checks, in exact integer arithmetic, a Nullstellensatz certificate
#     D * y^4 = sum_j N_j * m_j * f_{l_j}                                                        (*)
# with integers D != 0 and N_j, monomials m_j, and f_1..f_7 the cleared autocorrelations below.
# Julia Base only.  Run:
#     julia --startup-file=no exact/julia/c15_flat_certificate.jl [certificate.gp]
# The optional file receives D and the triples (l_j, m_j, N_j) for the independent PARI/GP check
# exact/pari/c15_flat_check.gp.
#
# mu_n = 2 cos(2 pi n/15) = 2 - lambda_n.  With t = mu_1 (minimal polynomial t^4 - t^3 - 4t^2 + 4t + 1),
# mu_n = D_n(t), D_0 = 2, D_1 = t, D_{n+1} = t D_n - D_{n-1} (2 cos na = D_n(2 cos a)); reduced and written
# in the Q-basis (1, mu_1, mu_2 = t^2 - 2, mu_3 = t^3 - 3t), mu_n has integer coordinates c(n).  So
# e^{it mu_n} = X^{c(n)} with X = (x, y, z, u) = (e^{it}, e^{it mu_1}, e^{it mu_2}, e^{it mu_3}), and X(t) is
# dense in T^4: every point of H_15 is (X^{c(n)})_n times a unimodular factor, X in T^4.  It is flat iff the
# autocorrelations R_l = sum_k X^{c(k+l) - c(k)}, l = 1..7, vanish (on T^4 the conjugate of X^a is X^{-a});
# f_l is R_l times the least monomial that makes it a polynomial.  By (*), y = 0 at every common zero of
# the f_l, and y != 0 on T^4.
#
# Method: (*) is a linear system in the coefficients of the multipliers of degree <= 18 - deg f_l (y^4 is
# not in the span below degree 18).  A nonsingular block is found modulo a prime P, solved over Q by
# P-adic (Dixon) lifting and rational reconstruction, and the resulting identity is checked over Z.

const P = 1048573      # prime < 2^20: products < 2^40, so lazy sums of <= 4711 of them stay below 2^53
const PF = Float64(P)  # and Float64 arithmetic on these integers is exact
const DEG = 18

say(s) = (println(s); flush(stdout))

function mu_coords()
    function red(a)  # reduce modulo t^4 = t^3 + 4t^2 - 4t - 1
        a = copy(a)
        while length(a) > 4
            c = pop!(a)
            k = length(a)
            a[k] += c
            a[k-1] += 4c
            a[k-2] -= 4c
            a[k-3] -= c
        end
        vcat(a, zeros(Int, 4 - length(a)))
    end
    D = Vector{Vector{Int}}(undef, 15)
    D[1] = [2, 0, 0, 0]
    D[2] = [0, 1, 0, 0]
    for n = 2:14
        D[n+1] = red(vcat(0, D[n]) .- vcat(D[n-1], 0))
    end
    [(a[1] + 2a[3], a[2] + 3a[4], a[3], a[4]) for a in D]
end

function generators()
    C = mu_coords()
    F = Vector{Vector{NTuple{4,Int}}}()
    for l = 1:7
        E = [C[mod(k+l, 15)+1] .- C[k+1] for k = 0:14]
        lo = ntuple(j -> minimum(e[j] for e in E), 4)
        push!(F, [e .- lo for e in E])
    end
    F
end

monos(d) = [(a, b, c, d - a - b - c) for a = 0:d for b = 0:(d-a) for c = 0:(d-a-b)]
invmod_p(a) = invmod(mod(Int(a), P), P)
modp(x::Float64) = mod(x, PF)

function reduce_row!(
    row::Vector{Float64},
    pivrows::Vector{Vector{Float64}},
    pivcol::Vector{Int},
)
    ne = length(row)
    for k in eachindex(pivrows)
        f = modp(row[pivcol[k]])
        f == 0 && continue
        pr = pivrows[k]
        @inbounds @simd for i = 1:ne
            row[i] = muladd(-f, pr[i], row[i])
        end
    end
    @inbounds for i = 1:ne
        row[i] = modp(row[i])
    end
    row
end

# row-reduce A^T mod P (rows = unknowns): pivot unknowns, pivot equations, and whether y^4 is in the span
function eliminate(cols::Vector{Vector{Int}}, ne::Int, target::Int)
    pivrows = Vector{Vector{Float64}}()
    pivcol = Int[]
    pivunk = Int[]
    row = zeros(Float64, ne)
    for j in eachindex(cols)
        fill!(row, 0.0)
        for i in cols[j]
            row[i] += 1
        end
        reduce_row!(row, pivrows, pivcol)
        c = findfirst(!iszero, row)
        c === nothing && continue
        iv = Float64(invmod_p(row[c]))
        push!(pivrows, [modp(v * iv) for v in row])
        push!(pivcol, c)
        push!(pivunk, j)
    end
    fill!(row, 0.0)
    row[target] = 1
    pivcol, pivunk, all(iszero, reduce_row!(row, pivrows, pivcol))
end

# LU with row pivoting of the dense block mod P, rows stored as vectors, L below the diagonal
function lu_modp!(S::Vector{Vector{Float64}})
    r = length(S)
    perm = collect(1:r)
    for k = 1:r
        @inbounds for j = k:r
            S[k][j] = modp(S[k][j])
        end
        if S[k][k] == 0
            q = findfirst(i -> modp(S[i][k]) != 0, (k+1):r)
            q === nothing && error("block singular mod P")
            q += k
            S[k], S[q] = S[q], S[k]
            perm[k], perm[q] = perm[q], perm[k]
            @inbounds for j = k:r
                S[k][j] = modp(S[k][j])
            end
        end
        iv = Float64(invmod_p(S[k][k]))
        prow = S[k]
        for i = (k+1):r
            si = S[i]
            f = modp(si[k])
            if f == 0
                si[k] = 0.0
                continue
            end
            lik = modp(f * iv)
            si[k] = lik
            @inbounds @simd for j = (k+1):r
                si[j] = muladd(-lik, prow[j], si[j])
            end
        end
    end
    perm, [Float64(invmod_p(S[i][i])) for i = 1:r]
end

function lusolve(S::Vector{Vector{Float64}}, dinv::Vector{Float64}, b::Vector{Int64})
    r = length(S)
    y = Float64.(mod.(b, P))
    for i = 1:r
        s = y[i]
        si = S[i]
        @inbounds @simd for k = 1:(i-1)
            s = muladd(-si[k], y[k], s)
        end
        y[i] = modp(s)
    end
    for i = r:-1:1
        s = y[i]
        si = S[i]
        @inbounds @simd for k = (i+1):r
            s = muladd(-si[k], y[k], s)
        end
        y[i] = modp(modp(s) * dinv[i])
    end
    Int64.(y)
end

function ratrec(a::BigInt, M::BigInt)
    bound = isqrt(M ÷ 2)
    r0, r1, s0, s1 = M, mod(a, M), BigInt(0), BigInt(1)
    while r1 > bound
        q = r0 ÷ r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    end
    (s1 == 0 || abs(s1) > bound) && return nothing
    return s1 > 0 ? (r1, s1) : (-r1, -s1)
end

# all entries over one growing common denominator D, so that N = D * c is integral
function common_reconstruction(acc::Vector{BigInt}, M::BigInt)
    half = M ÷ 2
    bound = isqrt(half)
    D = BigInt(1)
    N = Vector{BigInt}(undef, length(acc))
    for i in eachindex(acc)
        a = mod(D * acc[i], M)
        a > half && (a -= M)
        if abs(a) > bound
            rr = ratrec(mod(D * acc[i], M), M)
            rr === nothing && return nothing, N
            num, den = rr
            D *= den
            for k = 1:(i-1)
                N[k] *= den
            end
            a = num
        end
        N[i] = a
        D > bound && return nothing, N
    end
    D, N
end

function main(outfile)
    t0 = time()
    F = generators()
    @assert all(length(unique(f)) == 15 for f in F)
    degs = [maximum(sum, f) for f in F]
    eqs = vcat([monos(d) for d = 0:DEG]...)
    eqidx = Dict(m => i for (i, m) in enumerate(eqs))
    unk = [(l, m) for l = 1:7 for d = 0:(DEG-degs[l]) for m in monos(d)]
    target = eqidx[(0, 4, 0, 0)]
    say(
        "deg f_l = $degs; $(length(unk)) multiplier coefficients, $(length(eqs)) monomials of degree <= $DEG",
    )
    cols = [[eqidx[m .+ e] for e in F[l]] for (l, m) in unk]   # sparse columns, all entries 1

    pivcol, pivunk, inspan = eliminate(cols, length(eqs), target)
    r = length(pivcol)
    say(
        "rank mod P = $r, y^4 in the span mod P: $inspan ($(round(time() - t0, digits=1)) s)",
    )
    inspan || return false

    ridx = Dict(e => i for (i, e) in enumerate(pivcol))
    Scols = [[ridx[i] for i in cols[j] if haskey(ridx, i)] for j in pivunk]
    S = [zeros(Float64, r) for _ = 1:r]
    for (jj, sc) in enumerate(Scols), i in sc
        S[i][jj] += 1
    end
    perm, dinv = lu_modp!(S)

    b = zeros(Int64, r)
    haskey(ridx, target) && (b[ridx[target]] = 1)
    acc = [BigInt(0) for _ = 1:r]
    modulus = BigInt(1)
    step = 0
    nextcheck = 64
    while step < 8192
        x = lusolve(S, dinv, b[perm])
        for i = 1:r
            acc[i] += modulus * x[i]
        end
        for (jj, sc) in enumerate(Scols)    # b <- (b - S x) / P, exact
            xj = x[jj]
            xj == 0 && continue
            for i in sc
                b[i] -= xj
            end
        end
        @assert all(v -> mod(v, P) == 0, b)
        b .÷= P
        modulus *= P
        step += 1
        step < nextcheck && continue
        nextcheck *= 2
        D, N = common_reconstruction(acc, modulus)
        D === nothing && continue
        g = gcd(D, gcd(N))
        D ÷= g
        N .÷= g
        # the identity (*) over Z, on every monomial of degree <= DEG
        tot = Dict{Int,BigInt}()
        for (jj, j) in enumerate(pivunk)
            N[jj] == 0 && continue
            for i in cols[j]
                tot[i] = get(tot, i, BigInt(0)) + N[jj]
            end
        end
        get(tot, target, BigInt(0)) == D || continue
        all(p -> p.first == target || p.second == 0, tot) || continue
        nz = count(!iszero, N)
        say(
            "identity (*) holds over Z after $step lifting steps: $nz terms, D of $(ndigits(D)) digits, " *
            "|N_j| up to $(maximum(ndigits, N)) digits ($(round(time() - t0, digits=1)) s)",
        )
        if outfile !== nothing
            open(outfile, "w") do io
                println(io, D)
                for (jj, j) in enumerate(pivunk)
                    N[jj] == 0 && continue
                    l, m = unk[j]
                    println(io, "[", l, ",", join(m, ","), ",", N[jj], "]")
                end
            end
        end
        return true
    end
    false
end

ok = main(isempty(ARGS) ? nothing : ARGS[1])
say(
    ok ?
    "RESULT: PASS -- D*y^4 lies in the C_15 flatness ideal over Z, so H_15 has no flat point and B_15 < sqrt(15)" :
    "RESULT: FAIL -- no certificate found",
)
