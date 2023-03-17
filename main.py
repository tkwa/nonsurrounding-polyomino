# %%

from ortools.sat.python import cp_model

from pattern_generation import *

# %%

model = cp_model.CpModel()

class PolyominoSATInstance:

    def __init__(self, rows, cols, k):
        self.rows = rows
        self.cols = cols
        self.k = k
        self.model = cp_model.CpModel()
        self.constraint_descriptions = []

        """
        Makes a CPModel containing clauses that define a valid polyomino
        with bounding box of size (rows, cols) and
        maximum BFS depth of k.
        0=empty, 1=BFS root.
        """
        rows, cols, k = self.rows, self.cols, self.k
        model = self.model
        self.vars = [
            [model.NewIntVar(0, k, f"c{row}_{col}") for col in range(cols)]
            for row in range(rows)]
        self.ones = [
            [model.NewBoolVar(f"o{row}_{col}") for col in range(cols)]
            for row in range(rows)]
        self.zeros = [
            [model.NewBoolVar(f"z{row}_{col}") for col in range(cols)]
            for row in range(rows)]

        self.vm1_places = [[['?' for v in range(4)] for col in range(cols)] for row in range(rows)]

        # self.vars, self.ones, self.zeros, self.vm1_places = self.vars, self.ones, self.zeros, self.vm1_places
        # Which square is required to be v-1. 0, 1, 2, 3 = left, right, down, up
        for row in range(rows):
            for col in range(cols):
                if row > 0:
                    self.vm1_places[row][col][0] = model.NewBoolVar(f"vm1{row}_{col}_0")
                if row < rows - 1:
                    self.vm1_places[row][col][1] = model.NewBoolVar(f"vm1{row}_{col}_1")
                if col > 0:
                    self.vm1_places[row][col][2] = model.NewBoolVar(f"vm1{row}_{col}_2")
                if col < cols - 1:
                    self.vm1_places[row][col][3] = model.NewBoolVar(f"vm1{row}_{col}_3")
        
        for row in range(rows):
            for col in range(cols):
                model.Add(self.vars[row][col] == 1).OnlyEnforceIf(self.ones[row][col])
                model.Add(self.vars[row][col] == 0).OnlyEnforceIf(self.zeros[row][col])
                model.Add(self.vars[row][col] != 1).OnlyEnforceIf(self.ones[row][col].Not())
                model.Add(self.vars[row][col] != 0).OnlyEnforceIf(self.zeros[row][col].Not())
        
        # one 1 in the first column
        model.AddExactlyOne([self.ones[row][0] for row in range(rows)])
        # no 1s outside the first column
        for col in range(1, cols):
            for row in range(rows):
                model.AddBoolOr(self.ones[row][col].Not())

        # below a 1 in first column, we must have all 0s
        for i in range(rows): # higher numbers are higher in the grid
            for j in range(i):
                model.AddBoolOr(self.zeros[i][0], self.ones[j][0].Not())

        # any square not equal to 1 or 0 must have at least one neighbor with value v-1
        # this is the only place we use `ones`
        for row in range(rows):
            for col in range(cols):
                vars_to_check = []
                if row > 0:
                    vars_to_check.append(self.vm1_places[row][col][0])
                if row < rows - 1:
                    vars_to_check.append(self.vm1_places[row][col][1])
                if col > 0:
                    vars_to_check.append(self.vm1_places[row][col][2])
                if col < cols - 1:
                    vars_to_check.append(self.vm1_places[row][col][3])
                model.AddBoolOr(vars_to_check).OnlyEnforceIf(self.ones[row][col].Not(), self.zeros[row][col].Not())
                # model.AddBoolAnd([var.Not() for var in vars_to_check]).OnlyEnforceIf(ones[row][col])
                model.AddBoolAnd([var.Not() for var in vars_to_check]).OnlyEnforceIf(self.zeros[row][col])

        # BFS constraints
        # each square with value v can only have v, v+1, v-1, or 0 as neighbors
        # v-1 must be in any neighbor specified by vm1_places
        for row in range(rows):
            for col in range(cols):
                if row > 0:
                    model.AddLinearConstraint(
                        self.vars[row][col] - self.vars[row-1][col], -1, 1
                    ).OnlyEnforceIf(self.zeros[row][col].Not(), self.zeros[row-1][col].Not())
                    model.Add((self.vars[row-1][col] == self.vars[row][col] - 1)).OnlyEnforceIf(self.vm1_places[row][col][0])
                    model.Add(self.vars[row-1][col] != self.vars[row][col] - 1).OnlyEnforceIf(self.vm1_places[row][col][0].Not())
                if row < rows - 1:
                    model.AddLinearConstraint(
                        self.vars[row][col] - self.vars[row+1][col], -1, 1
                    ).OnlyEnforceIf(self.zeros[row][col].Not(), self.zeros[row+1][col].Not())
                    model.Add(self.vars[row+1][col] == self.vars[row][col] - 1).OnlyEnforceIf(self.vm1_places[row][col][1])
                    model.Add(self.vars[row+1][col] != self.vars[row][col] - 1).OnlyEnforceIf(self.vm1_places[row][col][1].Not())
                if col > 0:
                    model.AddLinearConstraint(
                        self.vars[row][col] - self.vars[row][col-1], -1, 1
                    ).OnlyEnforceIf(self.zeros[row][col].Not(), self.zeros[row][col-1].Not())
                    model.Add(self.vars[row][col-1] == self.vars[row][col] - 1).OnlyEnforceIf(self.vm1_places[row][col][2])
                    model.Add(self.vars[row][col-1] != self.vars[row][col] - 1).OnlyEnforceIf(self.vm1_places[row][col][2].Not())
                if col < cols - 1:
                    model.AddLinearConstraint(
                        self.vars[row][col] - self.vars[row][col+1], -1, 1
                    ).OnlyEnforceIf(self.zeros[row][col].Not(), self.zeros[row][col+1].Not())
                    model.Add(self.vars[row][col+1] == self.vars[row][col] - 1).OnlyEnforceIf(self.vm1_places[row][col][3])
                    model.Add(self.vars[row][col+1] != self.vars[row][col] - 1).OnlyEnforceIf(self.vm1_places[row][col][3].Not())
                
        self.constraint_descriptions.append('polyomino')

        # self.vars, self.ones, self.zeros, self.vm1_places = vars, ones, zeros, vm1_places
        # return self

    def add_max_weight_constraint(self, weight:int):
        # There must be `weight` nonzero values
        self.model.Add(
            sum(sum(self.zeros[row][col] for col in range(self.cols)) for row in range(self.rows)) >= self.rows*self.cols - weight)
        self.constraint_descriptions.append(f'max weight {weight}')
        
    def add_corner_constraint(self):
        # Corners must be 0
        rows, cols = self.rows, self.cols
        self.model.AddBoolAnd(
            [self.zeros[0][0], self.zeros[0][cols-1], self.zeros[rows-1][0], self.zeros[rows-1][cols-1]])
        self.constraint_descriptions.append('corners')

    def add_boundary_constraint(self, top_and_right = True):
        # Omino must touch boundary on all sides
        # left (column 0) unnecessary because there's a 1 there
        rows, cols = self.rows, self.cols
        # always enforce bottom
        self.model.AddBoolOr([self.zeros[0][col].Not() for col in range(cols)])
        if top_and_right:
            self.model.AddBoolOr([self.zeros[row][cols-1].Not() for row in range(rows)])
            self.model.AddBoolOr([self.zeros[rows-1][col].Not() for col in range(cols)])
        self.constraint_descriptions.append('boundary')

    def add_surround_constraint(self):
        # Special pattern.
        # There must not be a 0 surrounded by four nonzeros.
        for row in range(1, self.rows-1):
            for col in range(1, self.cols-1):
                self.model.AddBoolOr([
                    self.zeros[row][col].Not(),
                    self.zeros[row-1][col],
                    self.zeros[row+1][col],
                    self.zeros[row][col-1],
                    self.zeros[row][col+1]
                ])
        self.constraint_descriptions.append('surround')

    def add_s2_constraint(self):
        """
        There must not be this pattern:
        000000
        001...
        .1....
        ......
        Only check first/last row and column.
        5x5: 63k -> 81
        """
        # 1s in (row + 1, 0) and (row, 1)
        for row in range(self.rows-1):
            self.model.AddBoolOr(
                [self.zeros[row+1][0], self.zeros[row][1]] +
                [self.zeros[r][0].Not() for r in range(row+1)]
            )
            self.model.AddBoolOr(
                [self.zeros[row+1][-1], self.zeros[row][-2]] +
                [self.zeros[r][-1].Not() for r in range(row+1)]
            )
        
        for row in range(1, self.rows):
            self.model.AddBoolOr(
                [self.zeros[row-1][0], self.zeros[row][1]] +
                [self.zeros[r][0].Not() for r in range(row, self.rows)]
            )
            self.model.AddBoolOr(
                [self.zeros[row-1][-1], self.zeros[row][-2]] +
                [self.zeros[r][-1].Not() for r in range(row, self.rows)]
            )

        for col in range(self.cols-1):
            self.model.AddBoolOr(
                [self.zeros[0][col+1], self.zeros[1][col]] +
                [self.zeros[0][c].Not() for c in range(col+1)]
            )
            self.model.AddBoolOr(
                [self.zeros[-1][col+1], self.zeros[-2][col]] +
                [self.zeros[-1][c].Not() for c in range(col+1)]
            )
        
        for col in range(1, self.cols):
            self.model.AddBoolOr(
                [self.zeros[0][col-1], self.zeros[1][col]] +
                [self.zeros[0][c].Not() for c in range(col, self.cols)]
            )
            self.model.AddBoolOr(
                [self.zeros[-1][col-1], self.zeros[-2][col]] +
                [self.zeros[-1][c].Not() for c in range(col, self.cols)]
            )
    
    def pattern_to_constraint(self, pattern:Pattern):
        self.model.AddBoolOr(
            [self.zeros[i][j].Not() for (i,j) in pattern.bad] +
            [self.zeros[i][j] for (i,j) in pattern.good]
        )
    def add_all_copies_of_pattern(self,pattern:str):
        """
        Add a constraint based on a pattern (multiline string).
        The edges and corners define how the pattern should be extended.
        """
        """
        0000
        001.
        .1..
        ....
        """
        P = Pattern(pattern)
        for variant in P.all_isometric_copies_within_box(self.rows,self.cols):
            self.pattern_to_constraint(variant)




# %%


class VarArraySolutionPrinter(cp_model.CpSolverSolutionCallback):
    """Print intermediate solutions."""

    def __init__(self, pi:PolyominoSATInstance, do_print = True):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.pi = pi
        self.vars, self.ones, self.zeros, self.vm1_places = pi.vars, pi.ones, pi.zeros, pi.vm1_places
        self.__solution_count = 0
        self.print = do_print

    def on_solution_callback(self):
        self.__solution_count += 1
        if not self.print: return
        print('Solution %i' % self.__solution_count)
        for row in range(self.pi.rows)[::-1]: # print from top to bottom
            for col in range(self.pi.cols):
                v = self.Value(self.vars[row][col])
                print('#' if v else ' ', end='')
                o = self.Value(self.ones[row][col])
                z = self.Value(self.zeros[row][col])
                to_print = '?' if o and z else 'o' if o else 'z' if z else ' '
                # print(to_print, end='')
            print()
        print()

        # for row in range(2):
        #     for col in range(2):
        #         for dir in range(4):
        #             vm1 = self.vm1_places[row][col][dir]
        #             print('?' if isinstance(vm1, str) else self.Value(vm1), end='')
        #         print(' ', end='')
        #     print('\n')

    def solution_count(self):
        return self.__solution_count


# %%

omino_instance = PolyominoSATInstance(6,7, 15)
omino_instance.add_boundary_constraint(top_and_right=False)
# omino_instance.add_corner_constraint()
omino_instance.add_max_weight_constraint(24)
omino_instance.add_surround_constraint()
# omino_instance.add_s2_constraint()
for patterns in [
    corner_pattern,
    shift_pattern,
    shifted_U_pattern,
    diamond_contact_pattern,
    offset_corner_pattern,
    offset_corner_pattern_2,
    flipped_U_pattern,
    scoop_pattern,
    deep_scoop_pattern,
]:
    omino_instance.add_all_copies_of_pattern(patterns)

solver = cp_model.CpSolver()
solution_printer = VarArraySolutionPrinter(omino_instance, do_print=True)
# Enumerate all solutions.
solver.parameters.enumerate_all_solutions = True
# Solve.
status = solver.Solve(omino_instance.model, solution_printer)

print('Status = %s' % solver.StatusName(status))
print('Number of solutions found: %i' % solution_printer.solution_count())
print(f'Time: {solver.WallTime():.3f}')
print(f"Info: {solver.ResponseStats()}")

# %%
