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

    def valid_c2_placements(self):
        """Helper function for below. Returns a list of three coordinate pairs (center, side1, side2)."""
        rows, cols = self.rows, self.cols
        for row in range(rows):
            for col in range(cols):
                if row > 0 and col > 0:
                    yield [(row, col), (row-1, col), (row, col-1)]
                if row > 0 and col < cols - 1:
                    yield [(row, col), (row-1, col), (row, col+1)]
                if row < rows - 1 and col > 0:
                    yield [(row, col), (row+1, col), (row, col-1)]
                if row < rows - 1 and col < cols - 1:
                    yield [(row, col), (row+1, col), (row, col+1)]

    def valid_c4_placements(self, all_centers=True):
        """Helper function for below. Returns a list of two coordinate pairs (center, side).
        If all_centers=False, we only check centers off the edge of the grid.
        """
        rows, cols = self.rows, self.cols
        for row in range(rows):
            for col in range(cols):
                if row > 0:
                    yield [(row, col), (row-1, col)]
                if row < rows - 1:
                    yield [(row, col), (row+1, col)]
                if col > 0:
                    yield [(row, col), (row, col-1)]
                if col < cols - 1:
                    yield [(row, col), (row, col+1)]
                    

    def c2_coords_to_check(self, center):
        """Helper function for below. Returns a list of pairs of points."""
        dims =  self.rows, self.cols
        xr = min(center[1], dims[1] - center[1] - 1)
        yr = min(center[0], dims[0] - center[0] - 1)
        for xshift in range(-xr,xr+1):
            for yshift in range(0, yr+1):
                if xshift <= 0 and yshift == 0:
                    continue
                else:
                    yield (center[0] + yshift, center[1] + xshift), (center[0] - yshift, center[1] - xshift)

    def c2_constraints(self):
        """
        Adds constraints corresponding to partitions with C2 symmetry.
        ......
        .01...
        .1....
        ......
        For each placement that fits the above pattern,
        a nonsurrounding omino must contain both cells from at least one pair of cells symmetric about the 0.
        Since the constraint is a disjunction of conjunctions, we need to create new variables.
        """
        for center, side1, side2 in self.valid_c2_placements():
            conjunction_vars = []
            for p1, p2 in self.c2_coords_to_check(center):
                y1, x1 = p1
                y2, x2 = p2
                new_var = self.model.NewBoolVar(f'c2_{x1}_{y1}_{x2}_{y2}')
                self.model.AddBoolOr(self.zeros[y1][x1], self.zeros[y2][x2]).OnlyEnforceIf(new_var.Not())
                self.model.AddBoolAnd(self.zeros[y1][x1].Not(), self.zeros[y2][x2].Not()).OnlyEnforceIf(new_var)
                conjunction_vars.append(new_var)
            yc, xc = center
            ys1, xs1= side1
            ys2, xs2 = side2
            self.model.AddBoolOr(conjunction_vars).OnlyEnforceIf(
                self.zeros[yc][xc], self.zeros[ys1][xs1].Not(), self.zeros[ys2][xs2].Not())
    
    def c4_coords_to_check(self, center):
        """
        Helper function for below. Returns a list of pairs of points.

        """
        dims =  self.rows, self.cols
        xr = min(center[1], dims[1] - center[1] - 1)
        yr = min(center[0], dims[0] - center[0] - 1)
        for xshift in range(-xr,xr+1):
            for yshift in range(0, yr+1):
                if xshift <= 0 and yshift == 0:
                    continue
                else:
                    yield (center[0] + yshift, center[1] + xshift), (center[0] - yshift, center[1] - xshift)

    def c4_constraints(self, all_centers=False):
        """
        Adds constraints corresponding to partitions with C4 symmetry.
        ......
        .0....
        .1....
        ......
        For each placement that fits the above pattern,
        a nonsurrounding omino must contain both cells from at least one pair of cells rotationally symmetric about the 0.
        Since the constraint is a disjunction of conjunctions, we need to create new variables.
        This time, variables will be reused, so we'll track a dictionary of them.

        Also note that 
        """
        for center, side in self.valid_c2_placements():
            conjunction_vars = []
            for p1, p2 in self.c2_coords_to_check(center):
                y1, x1 = p1
                y2, x2 = p2
                new_var = self.model.NewBoolVar(f'c2_{x1}_{y1}_{x2}_{y2}')
                self.model.AddBoolOr(self.zeros[y1][x1], self.zeros[y2][x2]).OnlyEnforceIf(new_var.Not())
                self.model.AddBoolAnd(self.zeros[y1][x1].Not(), self.zeros[y2][x2].Not()).OnlyEnforceIf(new_var)
                conjunction_vars.append(new_var)
            yc, xc = center
            ys1, xs1= side
            self.model.AddBoolOr(conjunction_vars).OnlyEnforceIf(
                self.zeros[yc][xc], self.zeros[ys1][xs1].Not())
    
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
omino_instance.c2_constraints()
for patterns in [
    # shift_pattern, these three already covered by c2_constraints
    # scoop_pattern,
    # diagonal_bump_pattern,
    flipped_U_pattern,
    corner_pattern,
    shifted_U_pattern,
    diamond_contact_pattern,
    offset_corner_pattern,
    offset_corner_pattern_2,
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
# print size of model
print(f"Model size: {omino_instance.model.Proto().ByteSize()}")

# %%
