using Printf

abstract type Compound end

struct Gas <: Compound
    name::String
end

struct Liquid <: Compound
    name::String
end

struct Solid <: Compound
    name::String
end

struct Solute <: Compound
    name::String
end

name(c::Compound) = c.name
statesymbol(::Gas) = "(g)"
statesymbol(::Liquid) = "(l)"
statesymbol(::Solid) = "(s)"
statesymbol(::Solute) = "(aq)"

function Base.show(io::IO, c::Compound)
    print(io, name(c), " ", statesymbol(c))
end

struct Reaction
    reactants::Vector{Pair{Compound,Int}}
    products::Vector{Pair{Compound,Int}}

    kf::Float64
    kr::Float64

    orders_f::Vector{Int}
    orders_r::Vector{Int}
end

reg_compound::Regex = r"(?:(\d+) )?([\S]+) \((s|l|g|aq)\)"

function parse_compound(s)
    m = match(reg_compound, s)

    if isnothing(m)
        throw("Illegal compound string: $(s)")
    end

    n = if !isnothing(m[1])
        parse(Int, m[1])
    else
        1
    end

    if m[3] == "g"
        Gas(m[2]) => n
    elseif m[3] == "l"
        Liquid(m[2]) => n
    elseif m[3] == "s"
        Solid(m[2]) => n
    elseif m[3] == "aq"
        Solute(m[2]) => n
    else
        throw("Illegal phase: $(m[3])")
    end
end

reg_reaction::Regex = r"(?:(?:\d+ )?[\S]+ \((?:s|l|g|aq)\))|(?: \+ )|(?: <?->? )"

# CO₃²⁻ (aq) + 2 H⁺ (aq) <-> CO₂ (g) + H₂O (l)
function parse_reaction(s)
    reactants = Pair{Compound,Int}[]
    products = Pair{Compound,Int}[]

    cur_list = reactants

    expect_term = true

    arrows = (" -> ", " <- ", " <-> ")

    for m in eachmatch(reg_reaction, s)
        if !expect_term && m.match ∈ arrows
            cur_list = products
        elseif !expect_term
            if m.match != " + "
                throw("Illegal separator: \"$(m.match)\"")
            end
        elseif expect_term
            push!(cur_list, parse_compound(m.match))
        else
            throw("Illegal term: '\"$(m.match)\"")
        end

        expect_term = !expect_term
    end

    reactants, products
end

# K = kf / kr

function Reaction(s;
    K=nothing, kf=nothing, kr=nothing,
    orders_f=nothing, orders_r=nothing)

    if !isnothing(K) && !isnothing(kf) && !isnothing(kr)
        throw("Can't specify both rates and eq. constant at the same time")
    elseif !isnothing(K) && isnothing(kf) && isnothing(kr)
        kf = sqrt(K)
        kr = 1 / kf
    elseif !isnothing(K) && !isnothing(kf) && isnothing(kr)
        kr = kf / K
    elseif !isnothing(K) && isnothing(kf) && !isnothing(kr)
        kf = K / kr
    elseif isnothing(K) && !isnothing(kf) && !isnothing(kr)
        # K = kf / kr
    else
        throw("Not enough data specified to determine rates")
    end

    rx, px = parse_reaction(s)

    if isnothing(orders_f)
        orders_f = [n for (_, n) in rx]
    end

    if isnothing(orders_r)
        orders_r = [n for (_, n) in px]
    end

    Reaction(rx, px, kf, kr, orders_f, orders_r)
end

function show_reactionside(io::IO, side)
    isfirst = true
    for (c, n) in side
        if !isfirst
            print(io, " + ")
        end
        if n != 1
            print(io, n, ' ')
        end
        print(io, name(c), ' ', statesymbol(c))
        isfirst = false
    end
end

function Base.show(io::IO, r::Reaction)
    show_reactionside(io, r.reactants)

    if iszero(r.kr)
        print(io, " → ")
    elseif iszero(r.kf)
        print(io, " ← ")
    elseif !iszero(r.kf) && !iszero(r.kr)
        print(io, " ⇋  ")
    end

    show_reactionside(io, r.products)

    @printf io "  kf = %.2e, kr = %.2e" r.kf r.kr
end

abstract type GasContainer end

struct ConstantVolume <: GasContainer
    volume::Float64
    temperature::Float64
end

struct ConstantPressure <: GasContainer
    pressure::Float64
    temperature::Float64
end

function Base.show(io::IO, c::ConstantVolume)
    @printf io "Gas at constant volume: V = %.2f L, T = %.2f K" c.volume c.temperature
end

function Base.show(io::IO, c::ConstantPressure)
    @printf io "Gas at constant pressure: P = %.2f atm, T = %.2f" c.volume c.temperature
end

# Assuming unchanging molar fraction of liquids and solids for now
struct ReactionNetwork
    reactions::Vector{Reaction}
    compounds::Vector{Compound}

    solution_volume::Float64
    gascontainer::GasContainer
end

function Base.show(io::IO, network::ReactionNetwork)
    println(io, "Reactions:")
    for (i, rx) in enumerate(network.reactions)
        @printf io "(%2d) " i
        println(io, rx)
    end

    println(io, "\nCompounds:")
    for (i, c) in enumerate(network.compounds)
        @printf io "(%2d) " i
        println(io, c)
    end

    @printf io "Solution volume: %.2f L\n" network.solution_volume
    println(io, network.gascontainer)
end
