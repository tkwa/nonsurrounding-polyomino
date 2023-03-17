def str_to_coords(s):
    output = [[] for _ in range(10)]
    l = s.split('\n')[::-1]
    for i in range(len(l)):
        for j in range(len(l[i])):
            if l[i][j] in '0123456789':
                output[int(l[i][j])].append((j,i-1))
    return output

def rotate(point,ccw_turns):
    for _ in range(ccw_turns):
        point = (-point[1],point[0])
    return point
class Pattern:
    def __init__(self,bad,good=None):
        if good is None:
            bad,good = str_to_coords(bad)[:2]
        assert len(good) > 0
        self.bad = sorted(bad)
        self.good = sorted(good)
    def __hash__(self):
        #this lets us track a bunch of patterns in nice datastructures
        return hash((tuple(self.bad),tuple(self.good)))
    def __eq__(self,other):
        return self.bad==other.bad and self.good==other.good
    def __repr__(self):
        return f"Pattern on {len(self.bad)} bad and {len(self.good)} good cells:\n\tbad  = {str(self.bad)[:50]}\n\tgood = {str(self.good)[:50]}"
    def clip(self,xa,xz,ya,yz):
        return Pattern(
            [(u,v) for (u,v) in self.bad  if xa<=u<xz and ya<=v<yz],
            [(u,v) for (u,v) in self.good if xa<=u<xz and ya<=v<yz]
        )
    def translate(self,x,y):
        return Pattern(
            [(u+x,v+y) for (u,v) in self.bad], 
            [(u+x,v+y) for (u,v) in self.good]
        )
    def rotate(self,ccw_turns):
        return Pattern([rotate(p,ccw_turns) for p in self.bad],[rotate(p,ccw_turns) for p in self.good])
    def reflect(self):
        return Pattern([(-u,v) for (u,v) in self.bad],[(-u,v) for (u,v) in self.good])
    def good_bbox(self):
        xs, ys = zip(*self.good)
        return [min(xs),max(xs),min(ys),max(ys)]
    def all_translations_within_box(self,box_x,box_y):
        #box includes points with x = 0, ..., box_x -1 and likewise for y
        xa,xz,ya,yz = self.good_bbox()
        x_width = xz - xa
        y_width = yz - ya
        output = []
        for x_shift in range(-xa, box_x-xz):
            for y_shift in range(-ya, box_y - yz):
                output.append(self.translate(x_shift,y_shift).clip(0,box_x,0,box_y))
        return output
    def all_isometric_copies_within_box(self,box_x,box_y):
        reflected = self.reflect()
        all_orientations = []
        for chiral_thing in [self,reflected]:
            for i in range(4):
                all_orientations.append(chiral_thing.rotate(i))
        all_copies = set()
        for orientation in all_orientations:
            all_copies = all_copies.union(set(orientation.all_translations_within_box(box_x,box_y)))
        return all_copies
    def to_clause():
        #todo: integrate with ortools
        pass
    
#todo: experiment with larger buffers on patterns, see if that helps the solver.
corner_pattern = """
000000000000
01
0
0
0
0
0
0
"""
shift_pattern = """
000000000000
          10000000000000
           1
"""
shifted_U_pattern="""
000000000000
           000000000000
          101
           1
"""
diamond_contact_pattern = """
       000
      00100
     00   00
    00     00
   00       00
  00         00
 00           00
00             00
"""
big_corner_pattern = """
000000000000000000
                1
                1
                111
"""
#code to extract things: l = list(Pattern(shift_pattern).all_isometric_copies_within_box(10,10))