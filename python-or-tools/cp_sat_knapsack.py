# cp_sat_knapsack.py
from ortools.sat.python import cp_model

values = [20, 30, 66, 40, 60]
weights = [2, 3, 6, 4, 5]
capacity = 10
n = len(values)

model = cp_model.CpModel()
x = [model.NewBoolVar(f'x{i}') for i in range(n)]

# capacity constraint
model.Add(sum(weights[i] * x[i] for i in range(n)) <= capacity)

# objective: maximize total value
model.Maximize(sum(values[i] * x[i] for i in range(n)))

solver = cp_model.CpSolver()
solver.parameters.max_time_in_seconds = 10  # optional time limit
status = solver.Solve(model)

if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    picked = [i for i in range(n) if solver.Value(x[i]) == 1]
    print("Picked items:", picked)
    print("Total value:", solver.ObjectiveValue())
else:
    print("No solution found; status:", status)
