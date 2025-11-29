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
        adj = [[] for _ in range(n)]
        ce = []
        for e in self.E:
            u = index_of[e[0]]
            v = index_of[e[1]]
            E.append([u, v])
            adj[u].append(v)
            adj[v].append(u)
            ce.append(self.ce[e])

        solver = pywraplp.Solver.CreateSolver('GLOP')

        # Variables
        x = [solver.NumVar(0, 1, f"x{i}") for i in range(m)]


        # Constraints
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
        
        for bset in range(1<<n):
            if bin(bset).count('1') <= 1:
                continue
            constraint = solver.Constraint(0, bin(bset).count('1') - 1, "cs"+str(bset))
            for i in range(m):
                u = E[i][0]
                v = E[i][1]
                if (bset & (1<<u)) and (bset & (1<<v)):
                    constraint.SetCoefficient(x[i], 1)

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

        return objective.Value(), non_zero_edges
    
        
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