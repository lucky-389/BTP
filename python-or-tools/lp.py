from ortools.init.python import init
from ortools.linear_solver import pywraplp

m = 9 # edges
n = 10 # nodes
costs = [2, 3, 6, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21] # costs
E = [[i,i+1] for i in range(1,n+1)]
W = [[2,i] for i in range(1,n+1)]
adj = [[] for _ in range(n+1)]
for e in E:
    u = e[0]
    v = e[1]
    adj[u].append(v)
    adj[v].append(u)


solver = pywraplp.Solver.CreateSolver('GLOP')
if not solver:
    print('Solver not created.')
    exit(1)

# Variables

for i in range(1,m+1):
    # x1, x2, ..., xE
    vars()['x'+str(i)] = solver.NumVar(0, 1, 'x'+str(i))

# Constraints

inf = solver.infinity()
# sigma xe for each edge e belongs to v <= b(v)
for i in range(len(W)):
    u = W[i][1]
    bu = W[i][0]
    constraint = solver.Constraint(0, bu, "cb"+str(i+1))
    for v in adj[u]:
        if [u,v] in E:    
            idx = E.index([u,v]) + 1
            constraint.SetCoefficient(vars()['x'+str(idx)], 1)
        else:
            idx = E.index([v,u]) + 1
            constraint.SetCoefficient(vars()['x'+str(idx)], 1)

# sigma xe for each edge e belongs to E(S) <= |S| - 1 for each S subset of V with |S| >= 2
for bset in range(1<<n):
    if bin(bset).count('1') <= 1:
        continue
    constraint = solver.Constraint(0, bin(bset).count('1') - 1, "cs"+str(bset))
    for i in range(m):
        u = E[i][0]
        v = E[i][1]
        if (bset & (1<<u)) and (bset & (1<<v)):
            constraint.SetCoefficient(vars()['x'+str(i+1)], 1)

# complete tree constraint
constraint = solver.Constraint(n-1, n-1, "ct")
for i in range(1,m+1):
    constraint.SetCoefficient(vars()['x'+str(i)], 1)

# Objective
objective = solver.Objective()
for i in range(1,m+1):
    objective.SetCoefficient(vars()['x'+str(i)], costs[i-1])
objective.SetMinimization()


# Solve
print(f"Solving with {solver.SolverVersion()}")
result_status = solver.Solve()
print(f"Status: {result_status}")
if result_status != pywraplp.Solver.OPTIMAL:
    print("The problem does not have an optimal solution!")
    if result_status == pywraplp.Solver.FEASIBLE:
        print("A potentially suboptimal solution was found")
    else:
        print("The solver could not solve the problem.")
        exit(1)

# Display the solution
print("Objective value =", objective.Value())
for i in range(1,m+1):
    print(f"x{i} = {vars()['x'+str(i)].solution_value()}")

