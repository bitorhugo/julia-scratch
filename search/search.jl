using DataStructures

const Graph = Dict{Int64, Vector{Tuple{Int64, Float64}}}

function reconstruct_path(paths, source, target)

    source == target && return [source]

    !haskey(paths, target) && return

    trace = [target]
    current = target

    while paths[current] != source && isfinite(paths[current])
        current = paths[current]
        push!(trace, current)
    end

    push!(trace, source)

    return reverse(trace)
end

function h(node, target, coords)
    xn, yn = coords[node]
    xt, yt = coords[target]
    return sqrt((xt - xn)^2 + (yt - yn)^2)
end

function astar(g::Dict{T, Vector{Tuple{T, W}}},
                      source::T,
                      target::T,
                      coords::Dict{Int64, Tuple{W, W}}) where {T, W<:Real}

    q::PriorityQueue{Int, Float64} = PriorityQueue(source => h(source, target, coords))
    distances::Dict{Int, Float64} = Dict(k => Inf for k in keys(g))
    visited::BitSet = BitSet()
    hops::Dict{Int, Int} = Dict()

    distances[source] = 0

    while !isempty(q)
        current_node, current_cost = dequeue_pair!(q)

        if current_node == target
            break
        end

        for (neighbor, edge_cost) in g[current_node]
            additive_cost = distances[current_node] + edge_cost

            if !(neighbor in visited)
                if distances[neighbor] > additive_cost
                    distances[neighbor] = additive_cost
                    hops[neighbor] = current_node
                    q[neighbor] = additive_cost + h(neighbor, target, coords) # heuristically speaking
                end
            end
        end

        push!(visited, current_node)
    end

    return (; distances, hops, visited)
end

function dijkstra(g::Dict{T, Vector{Tuple{T, W}}}, source::T) where {T, W<:Real}
    q::PriorityQueue{Int, Float64} = PriorityQueue(source => 0)
    distances::Dict{Int, Float64} = Dict(k => Inf for k in keys(g))
    visited::BitSet = BitSet()
    hops::Dict{Int, Int} = Dict()

    while !isempty(q)
        current_node, current_cost = dequeue_pair!(q)

        distances[current_node] = current_cost

        for (neighbor, edge_cost) in g[current_node]
            additive_cost = current_cost + edge_cost

            if !(neighbor in visited)
                if distances[neighbor] > additive_cost
                    distances[neighbor] = additive_cost
                    hops[neighbor] = current_node
                    q[neighbor] = additive_cost
                end
            end
        end

        push!(visited, current_node)
    end

    return (; distances, hops, visited)
end

# linear
function shortpath(g, source)
    dist::Dict{Int, Float64} = Dict([k => k == source ? 0 : Inf for (k, v) in g])
    visited::BitSet = BitSet()
    prev::Dict{Int, Int} = Dict()

    function relax(u, v, w)
        cost = dist[u] + w
        if dist[v] > cost
            dist[v] = cost
            prev[v] = u
        end
    end

    while length(visited) < length(dist)
        u = argmin(k -> dist[k], [k for k in keys(dist) if !(k in visited)])

        if !isfinite(dist[u])
            break
        end

        push!(visited, u)

        for (v, w) in g[u]
            if !(v in visited)
                relax(u, v, w)
            end
        end
    end

    return (dist, prev)
end

# programs

g = Graph(
    1 => [(2, 1.0), (4, 1.0)],
    2 => [(1, 1.0), (3, 1.0)],
    3 => [(2, 1.0), (6, 1.0)],
    4 => [(1, 1.0), (7, 1.0)],
    5 => [],
    6 => [(3, 1.0), (9, 1.0)],
    7 => [(4, 1.0), (8, 1.0)],
    8 => [(7, 1.0), (9, 1.0)],
    9 => [(6, 1.0), (8, 1.0)],
)

## to use in heurisitic search
coords = Dict(
    1 => (0.0, 2.0), 2 => (1.0, 2.0), 3 => (2.0, 2.0),
    4 => (0.0, 1.0), 5 => (1.0, 1.0), 6 => (2.0, 1.0),
    7 => (0.0, 0.0), 8 => (1.0, 0.0), 9 => (2.0, 0.0),
)
