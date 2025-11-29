from ortools.sat.python import cp_model

n = 8

nodes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

edges = [
    [0, 1, 4],
    [0, 2, 6],
    [0, 3, 5],
    [1, 4, 7],
    [2, 5, 3],
    [2, 6, 8],
    [5, 7, 2],
    [7, 8, 4],
    [8, 9, 9],
    [1, 2, 10],
    [1, 5, 12],
    [3, 4, 11],
    [3, 6, 14],
    [4, 5, 6],
    [6, 7, 15],
    [2, 8, 9],
    [0, 4, 13],
    [0, 5, 5],
    [9, 5, 7],
    [3, 9, 16]
]

degree_bounds = [3, 2, 3, 1, 1, 2, 1, 2, 2, 1]

model = cp_model.CpModel()

x = [model.NewBoolVar(f'x{i}') for i in range(len(edges))]
model