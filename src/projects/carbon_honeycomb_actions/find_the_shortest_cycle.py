import numpy as np
from collections import deque


def find_all_shortest_cycles(common_points_matrix: np.ndarray) -> list[list[int]]:
    n: int = common_points_matrix.shape[0]

    # Helper to find shortest cycle length
    def find_shortest_cycle_length() -> int | float:
        shortest = float('inf')
        for start in range(n):
            dist = [-1] * n
            parent = [-1] * n
            dist[start] = 0
            queue = deque([start])

            while queue:
                u = queue.popleft()
                for v in range(n):
                    if not common_points_matrix[u, v]:
                        continue
                    if dist[v] == -1:
                        # Tree edge
                        dist[v] = dist[u] + 1
                        parent[v] = u
                        queue.append(v)
                    elif parent[u] != v and parent[v] != u:
                        # Back-edge found indicating a cycle
                        # Cycle length = dist[u] + dist[v] + 1
                        cycle_length = dist[u] + dist[v] + 1
                        if cycle_length < shortest:
                            shortest = cycle_length
        return shortest if shortest != float('inf') else -1

    shortest_cycle_length = find_shortest_cycle_length()
    if shortest_cycle_length == -1:
        # No cycles found
        return []

    # Helper to reconstruct cycle from a back-edge
    def reconstruct_cycle(u, v, parent_u, parent_v):
        # u and v are endpoints where back-edge was found
        # We know dist[u] + dist[v] + 1 = shortest_cycle_length
        cycle_path_u = []
        cycle_path_v = []

        # Trace back u to start
        cur = u
        while cur != -1:
            cycle_path_u.append(cur)
            cur = parent_u[cur]

        # Trace back v to start
        cur = v
        while cur != -1:
            cycle_path_v.append(cur)
            cur = parent_v[cur]

        # Reverse paths to get from root to u and root to v
        cycle_path_u.reverse()
        cycle_path_v.reverse()

        # Find the first common node in cycle_path_u and cycle_path_v
        # Actually, we know they share no direct intersection except at the cycle closing
        # A simpler approach: since the cycle is formed by joining these two paths plus the edge (u,v)
        # We'll find from start forward how they diverge
        i = 0
        while i < len(cycle_path_u) and i < len(cycle_path_v) and cycle_path_u[i] == cycle_path_v[i]:
            i += 1

        # Now cycle is formed by the unique part of cycle_path_u + reverse of unique part of cycle_path_v plus the connecting edge
        # Actually, we must be careful to ensure the cycle is correct:
        # The cycle runs as u->...->root->...->v and then directly v->u forms the cycle.
        # We'll take the portion from i-1 onward because i-1 is the last common node
        # Actually, let's do a direct approach:
        # The shortest cycle is: the path from u up to common ancestor + path from v up to common ancestor + edge (u,v)
        # But we recognized a back-edge in BFS from a single source; to simplify, we do BFS from each node separately.

        # Let's do a simpler approach: Just run BFS from each node, store parent arrays,
        # and whenever we find a back-edge (u,v), reconstruct from u->start and v->start, find intersection:
        ancestor_set = set(cycle_path_u)
        # Find first intersection in cycle_path_v from start:
        # The earliest intersection point is the LCA
        intersect_node = None
        for node in cycle_path_v:
            if node in ancestor_set:
                intersect_node = node
                break

        # Now form cycle: from intersect_node to u, then (u,v), then from v to intersect_node
        # Extract subpath from intersect_node to u:
        iu = cycle_path_u.index(intersect_node)
        iv = cycle_path_v.index(intersect_node)

        subpath_u = cycle_path_u[iu:]  # path from intersect_node to u
        subpath_v = cycle_path_v[iv:]  # path from intersect_node to v
        # The cycle goes: intersect_node -> ... -> u -> v -> ... -> intersect_node
        # But currently subpath_u is from intersect_node to u, we need it to go from u to intersect_node reversed.
        # Actually, we know u was the endpoint of the back-edge. The cycle was detected at edge (u,v), so:
        # cycle = subpath_u (from intersect_node to u) + v + reverse(subpath_v except intersect_node)
        # Wait, we must be careful:
        # subpath_u: intersect_node ... u
        # subpath_v: intersect_node ... v
        # The cycle formed: u->...->intersect_node->...->v->u
        # So the cycle is subpath_u going forward from intersect_node to u,
        # then add v, then add subpath_v from v back to intersect_node in reverse order.

        # Actually, we want a consistent direction. Let's say we start cycle at intersect_node:
        # cycle: intersect_node->...->u->v->...->intersect_node
        # subpath_u: [intersect_node, ..., u]
        # subpath_v: [intersect_node, ..., v]
        # If we go forward along subpath_u (intersect_node to u), then edge (u,v),
        # then we want to go from v back to intersect_node. But subpath_v goes forward from intersect_node to v.
        # To go from v back to intersect_node, we reverse subpath_v and skip intersect_node:
        cycle = subpath_u[:]  # intersect_node ... u
        cycle.append(v)
        reversed_subpath_v = subpath_v[::-1]
        # reversed_subpath_v: v ... intersect_node
        # skip the intersect_node at the end since it's already in the start
        # add v->...->intersect_node but we've added v already, so skip first element?
        cycle.extend(reversed_subpath_v[1:])

        # Now cycle might have duplicates if we're not careful. Let's ensure unique by checking length:
        # The BFS ensures no repeated node in a path. Since we carefully constructed from parent arrays, no repetition should occur.

        return cycle

    # Next, we find all shortest cycles
    shortest_cycles = set()

    def record_shortest_cycles_from_start(start):
        dist = [-1] * n
        parent = [-1] * n
        dist[start] = 0
        queue = deque([start])

        while queue:
            u = queue.popleft()
            for v in range(n):
                if not common_points_matrix[u, v]:
                    continue
                if dist[v] == -1:
                    dist[v] = dist[u] + 1
                    parent[v] = u
                    queue.append(v)
                elif parent[u] != v and parent[v] != u:
                    # Found a cycle
                    cycle_length = dist[u] + dist[v] + 1
                    if cycle_length == shortest_cycle_length:
                        # Reconstruct cycle
                        # We need parent arrays from both ends. Let's do BFS from 'v' too to get its parent array
                        # Or we can do a trick: run BFS from each node separately and only reconstruct cycle from that BFS's perspective.

                        # We will run BFS from 'v' as well to get its parents:
                        dist2 = [-1]*n
                        parent2 = [-1]*n
                        dist2[v] = 0
                        q2 = deque([v])
                        # Run BFS from v just to get a parent chain
                        while q2:
                            x = q2.popleft()
                            for y in range(n):
                                if common_points_matrix[x, y] and dist2[y] == -1:
                                    dist2[y] = dist2[x] + 1
                                    parent2[y] = x
                                    q2.append(y)
                                    if y == u:
                                        # Once we reached u from v, we have both parent chains
                                        break

                        cycle = reconstruct_cycle(u, v, parent, parent2)
                        # Canonicalize cycle to avoid duplicates: sort and store as tuple
                        cycle = tuple(sorted(cycle))
                        shortest_cycles.add(cycle)

    for start in range(n):
        record_shortest_cycles_from_start(start)

    # Convert each tuple back to a list for the final result
    return [list(cycle) for cycle in shortest_cycles]
