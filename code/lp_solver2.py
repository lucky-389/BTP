from ortools.init.python import init
from ortools.linear_solver import pywraplp

class LPSolver:
    def __init__(self):
        self.n = 0
        self.m = 0
        self.V = set()
        self.E = set()  # (u,v) with u < v
        self.W = {}
        self.adj = {}
        self.ce = {}

    def add_node(self, v):
        self.V.add(v)
        self.adj[v] = set()
        self.n += 1
        
    def add_edge(self, u, v, cost=0):
        if u > v:
            u, v = v, u
        self.E.add((u, v))
        self.adj[u].add(v)
        self.adj[v].add(u)
        self.ce[(u, v)] = cost
        self.m += 1

    def remove_node(self, u):
        if u not in self.V:
            return
        self.V.remove(u)
        if u in self.W:
            self.W.pop(u)
        for v in self.adj[u]:
            self.adj[v].remove(u)
            if (u, v) in self.E:
                self.E.remove((u, v))
            else:
                self.E.remove((v, u))
        self.adj.pop(u)
        self.n -= 1

    def remove_edge(self, u, v):
        if u > v:
            u, v = v, u
        if (u, v) not in self.E:
            return
        self.E.remove((u, v))
        self.adj[u].remove(v)
        self.adj[v].remove(u)
        self.m -= 1
        
    def set_bound(self, v, bound):
        self.W[v] = bound

    def unset_bound(self, v):
        if v in self.W:
            self.W.pop(v)

    def get_bound(self, v):
        if v not in self.W:
            return None
        return self.W[v]
    
    def get_degree(self, v):
        return len(self.adj[v])

    # considering all subset of vertices for cut constraints
    def solve(self):
        n = self.n
        m = self.m
        V = list(self.V)
        index_of = {V[i]: i for i in range(n)}
        W = [self.W[V[i]] if (V[i] in self.W) else m for i in range(n)]
        E = []
        rE = []
        adj = [[] for _ in range(n)]
        ce = []
        edge_index = {}
        for e in self.E:
            u = index_of[e[0]]
            v = index_of[e[1]]
            E.append([u, v])
            rE.append([v, u])
            edge_index[(u,v)] = len(E)-1
            adj[u].append(v)
            adj[v].append(u)
            ce.append(self.ce[e])

        solver = pywraplp.Solver.CreateSolver('GLOP')

        # Variables
        x = [solver.NumVar(0, 1, f"x{i}") for i in range(m)]
        f = [[solver.NumVar(0,1, f"f{E[i]}_{j}") for j in range(n)] for i in range(m)]
        rf = [[solver.NumVar(0,1, f"f{rE[i]}_{j}") for j in range(n)] for i in range(m)]

        # bfs tree, towards root (includes all edges)
        root = 0
        # incoming = [[] for _ in range(n)]
        # outgoing = [[] for _ in range(n)]
        # queue = [root]
        # visited = [False]*n
        # visited[root] = True
        # while queue:
        #     u = queue.pop(0)
        #     for v in adj[u]:
        #         ei = edge_index[(u,v)] if (u,v) in edge_index else edge_index[(v,u)]
        #         if not visited[v]:
        #             visited[v] = True
        #             incoming[u].append(ei)
        #             outgoing[v].append(ei)
        #             queue.append(v)
        #         else:
        #             if ei not in outgoing[v] and ei not in incoming[v]:
        #                 outgoing[v].append(ei)
        #                 incoming[u].append(ei)

        # printing
        print("vertices:", V)
        print("edges:", E)

        # Constraints
        # degree bound to each node
        inf = solver.infinity()
        for u,bu in self.W.items():
            u = index_of[u]
            constraint = solver.Constraint(0, bu, "cb"+str(u))
            for v in adj[u]:
                if [u,v] in E:
                    idx = E.index([u,v])
                    constraint.SetCoefficient(x[idx], 1)
                else:
                    idx = E.index([v,u])
                    constraint.SetCoefficient(x[idx], 1)
        
        # flow constraints
        # fe_v <= xe
        for i in range(m):
            for j in range(n):
                constraint = solver.Constraint(-inf, 0, "cf1_"+str(i)+"_"+str(j))
                constraint.SetCoefficient(f[i][j], 1)
                constraint.SetCoefficient(x[i], -1)

        for i in range(m):
            for j in range(n):
                constraint = solver.Constraint(-inf, 0, "cf2_"+str(i)+"_"+str(j))
                constraint.SetCoefficient(rf[i][j], 1)
                constraint.SetCoefficient(x[i], -1)
        
        # sum of incoming flow = sum of outgoing flow
        for i in range(1,n):
            for j in range(n):
                if j == root:
                    constraint = solver.Constraint(1, 1, "cf3_"+str(i)+"_"+str(j))
                    # for ei in outgoing[j]:
                    #     constraint.SetCoefficient(f[ei][i], 1)
                    for ei in range(m):
                        if (E[ei][1] == root):
                            constraint.SetCoefficient(f[ei][i], 1)
                    for ei in range(m):
                        if (rE[ei][1] == root):
                            constraint.SetCoefficient(rf[ei][i], 1)
                    for ei in range(m):
                        if (E[ei][0] == root):
                            constraint.SetCoefficient(f[ei][i],-1)
                    for ei in range(m):
                        if (rE[ei][0] == root):
                            constraint.SetCoefficient(rf[ei][i],-1)
                elif j == i:
                    constraint = solver.Constraint(1, 1, "cf3_"+str(i)+"_"+str(j))
                    for ei in range(m):
                        if (E[ei][0] == j):
                            constraint.SetCoefficient(f[ei][i], 1)
                    for ei in range(m):
                        if (rE[ei][0] == j):
                            constraint.SetCoefficient(rf[ei][i], 1)
                    for ei in range(m):
                        if (E[ei][1] == j):
                            constraint.SetCoefficient(f[ei][i],-1)
                    for ei in range(m):
                        if (rE[ei][1] == j):
                            constraint.SetCoefficient(rf[ei][i],-1)
                else:
                    constraint = solver.Constraint(0, 0, "cf3_"+str(i)+"_"+str(j))
                    for ei in range(m):
                        if (E[ei][0] == j):
                            constraint.SetCoefficient(f[ei][i], 1)
                    for ei in range(m):
                        if (rE[ei][0] == j):
                            constraint.SetCoefficient(rf[ei][i], 1)
                    for ei in range(m):
                        if (E[ei][1] == j):
                            constraint.SetCoefficient(f[ei][i], -1)
                    for ei in range(m):
                        if (rE[ei][1] == j):
                            constraint.SetCoefficient(rf[ei][i], -1)

        
        # sum of x equals n-1
        constraint = solver.Constraint(n-1, n-1, "ct")
        for i in range(m):
            constraint.SetCoefficient(x[i], 1)

        # Objective
        objective = solver.Objective()
        for i in range(m):
            objective.SetCoefficient(x[i], ce[i])
        objective.SetMinimization()

        status = solver.Solve()

        non_zero_edges = []
        # print("Objective value =", objective.Value())
        for i in range(m):
            if x[i].solution_value() > 0.0:
                non_zero_edges.append([V[E[i][0]], V[E[i][1]]])
            print(f"x{i} = {x[i].solution_value()}")

        # print flow values
        for i in range(m):
            for j in range(n):
                if f[i][j].solution_value() != 0.0:
                    print(f"f{E[i]}_{j} = {f[i][j].solution_value()}")
        for i in range(m):
            for j in range(n):
                if rf[i][j].solution_value() != 0.0:
                    print(f"f{rE[i]}_{j} = {rf[i][j].solution_value()}")
        return objective.Value(), non_zero_edges
    

# if __name__ == "__main__":
#     lp = LPSolver()
#     n = 5
#     for i in range(1, n+1):
#         lp.add_node(i)
#     edges = [(1,2),(1,3),(2,3),(2,4),(3,4),(4,5)]
#     costs = [1,2,2,1,3,1]
#     for i in range(len(edges)):
#         lp.add_edge(edges[i][0], edges[i][1], costs[i])
#     lp.set_bound(1, 2)
#     lp.set_bound(2, 2)
#     lp.set_bound(3, 2)
#     lp.set_bound(4, 2)
#     lp.set_bound(5, 1)
#     print(lp.solve())
    

if __name__ == "__main__":
    lp = LPSolver()
    n = 4
    for i in range(1, n+1):
        lp.add_node(i)
    edges = [(1,2),(1,3),(2,3),(2,4),(3,4)]
    costs = [1,2,5,3,4]
    for i in range(len(edges)):
        lp.add_edge(edges[i][0], edges[i][1], costs[i])
    lp.set_bound(1, 1)
    lp.set_bound(2, 2)
    lp.set_bound(3, 1)
    lp.set_bound(4, 2)
    print(lp.solve())
        