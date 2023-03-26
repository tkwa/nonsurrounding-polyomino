# %%
import math
import random
from PIL import Image,ImageDraw
import IPython.display
import itertools
import time
import book
from fractions import Fraction as frac
import pickle
import pycosat
# %%
def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")
# %%
def dual(l):
    #Dual lattice of a given lattice
    if l=='square':
        return 'square'
    elif l=='triangle':
        return 'hexagon'
    elif l=='hexagon':
        return 'triangle'
def translation_lattice(l):
    #lattice along which points can translate
    if l=='square':
        return 'square'
    else:
        return 'triangle'
def lattice_rotations(l):
    if l=='square':
        return 4
    elif l=='triangle' or l=='hexagon':
        return 6
def lattice_poly(l,c=[0,0,0]):
    if l=='square':
        return [0,0,1,0,1,1,0,1]
    elif l=='triangle':
        return [0.5,-0.288675,0.5,0.288675,0,0.57735,-0.5,0.288675,-0.5,-0.288675,0,-0.57735]
    elif l=='hexagon':
        result = [-0.5,-0.288675,0.5,-0.288675,0,0.57735]
        if c[2]:
            result=[-x for x in result]
        return result
def neighbor_count(l):
    if l=='square':
        return 4
    elif l=='triangle':
        return 6
    elif l=='hexagon':
        return 3
# %%
def translate(P,v):
    #translate P by point/vector 
    if isinstance(P,Polyomino):
        return Polyomino(translate(P.cells,v),P.face_lattice)
    return [e+v for e in P]
def normalize(p):
    #Normalize a list of coordinates (an animal in the plane)
    #so that it is translated to have nonnegative coordinates
    lattice = p[0].lattice
    min_x=min([e.coords[0] for e in p])
    min_y=min([e.coords[1] for e in p])
    return translate(p,Point([-min_x,-min_y],translation_lattice(lattice)))
class Point:
    #Point in a given lattice
    #square: [x,y], normally
    #triangle: [x,y], where x=(1,0) and y=(0.5,sqrt(3)/2)
    #hexagon: [x,y,q], where (x,y) is a point in the triangular lattice and q is 0 or 1 for up/down triangle
    def __init__(self,coords,lattice='square'):
        self.coords=coords
        self.lattice=lattice.lower()
    def __neg__(self):
        if self.lattice=='hexagon':
            return Point([-self.coords[0]-1,-self.coords[1]-1,1-self.coords[2]],self.lattice)
        return Point([-x for x in self.coords],self.lattice)
    def __add__(self,other):
        new_coords = [self.coords[i]+other.coords[i] for i in range(2)]
        if self.lattice=='triangle' and other.lattice=='hexagon':
            return Point(new_coords+[other.coords[2]],'hexagon')
        elif self.lattice=='hexagon' and other.lattice=='triangle':
            return Point(new_coords+[self.coords[2]],'hexagon')
        elif self.lattice=='hexagon' and other.lattice=='hexagon':
            #might fail: hexagon isn't really a lattice!
            #we're aiming for a triangular output here, so want different orientations
            if self.coords[2]!=other.coords[2]:
                return Point([new_coords[i]+1 for i in range(2)],'triangle')
            else:
                return False
        else:
            return Point(new_coords,self.lattice)
    def __sub__(self,other):
        return self+(-other)
    def __str__(self):
        return ("Point(%s): %s" % (self.lattice,str(self.coords))) 
    def __eq__(self,other):
        return self.coords==other.coords and self.lattice==other.lattice
    def __gt__(self,other):
        return self.coords>other.coords
    def __lt__(self,other):
        return self.coords<other.coords
    def __ge__(self,other):
        return self.coords>=other.coords
    def __le__(self,other):
        return self.coords<=other.coords
    def rotate(self):
        #Rotate clockwise at an angle specific to the lattice
        c=self.coords
        if self.lattice=='square':
            return Point([c[1],-c[0]],'square')
        elif self.lattice=='triangle':
            return Point([c[0]+c[1],-c[0]],'triangle')
        elif self.lattice=='hexagon':
            return Point([c[0]+c[1]+c[2],-c[0]-1,1-c[2]],'hexagon')
    def flip(self):
        #flip about the x-axis
        c=self.coords
        if self.lattice=='square':
            return Point([-c[0],c[1]],'square')
        elif self.lattice=='triangle':
            return Point([-c[0]-c[1],c[1]],'triangle')
        elif self.lattice=='hexagon':
            return Point([-c[0]-c[1]-c[2],c[1],c[2]],'hexagon')
    def neighbors(self):
        c=self.coords
        n=[]
        if self.lattice=='square':
            n=[[c[0]+1,c[1]],[c[0]-1,c[1]],[c[0],c[1]+1],[c[0],c[1]-1]]
        elif self.lattice=='triangle':
            n=[[c[0]+1,c[1]],[c[0]-1,c[1]],[c[0],c[1]+1],[c[0],c[1]-1],[c[0]-1,c[1]+1],[c[0]+1,c[1]-1]]
        elif self.lattice=='hexagon':
            n=[[c[0],c[1],1-c[2]],[c[0]+c[2],c[1]-1+c[2],1-c[2]],[c[0]-1+c[2],c[1]+c[2],1-c[2]]]
        return [Point(e,self.lattice) for e in n]
    def dual_neighbors(self):
        #Points in the dual around the point
        c=self.coords
        n=[]
        if self.lattice=='square':
            n=[[c[0],c[1]],[c[0]+1,c[1]],[c[0]+1,c[1]+1],[c[0],c[1]+1]]
        elif self.lattice=='triangle':
            vs=[[0,0,0],[-1,0,0],[-1,0,1],[-1,-1,1],[0,-1,0],[0,-1,1]]
            n=[[c[0]+v[0],c[1]+v[1],v[2]] for v in vs]
        elif self.lattice=='hexagon':
            n=[[c[0],c[1]+1],[c[0]+1,c[1]],[c[0]+c[2],c[1]+c[2]]]
        return [Point(e,dual(self.lattice)) for e in n]
    def x(self):
        c=self.coords
        if self.lattice=='square':
            return c[0]
        elif self.lattice=='triangle':
            return c[0]+c[1]/2
        elif self.lattice=='hexagon':
            return c[0]+c[1]/2+[0.5,1][c[2]]
    def y(self):
        c=self.coords
        if self.lattice=='square':
            return c[1]
        elif self.lattice=='triangle':
            return c[1]*0.866
        elif self.lattice=='hexagon':
            return (c[1]+(1+c[2])/3)*0.866
    def key(self):
        return tuple(self.coords)
    def __hash__(self):
        return hash(self.key())
    def draw(self,draw,translate,scale,fill='red',outline='white'): #Draw the cell corresponding to the point
        poly = lattice_poly(self.lattice,self.coords)
        pos=[self.x(),self.y()]
        true_poly = [scale*(poly[i]+pos[i%2])*[1,-1][i%2] + translate[i%2] for i in range(len(poly))]
        draw.polygon(true_poly,fill=fill,outline=outline)
def pointify(x,lattice):
    #convert coords, or a list of coords, into points
    if isinstance(x,Point):
        return x
    elif type(x)==list and len(x)>0 and type(x[0])==int:
        return Point(x,lattice)
    elif type(x)==list:
        return [pointify(e,lattice) for e in x]
    else:
        return False

# %%
#defining the polyomino class
#NOTE: If you rerun this cell, things may break in weird ways, 
#because an old-Polyomino object and a new-Polyomino object won't be equal.

#mapping from strings to polyominoes
polyomino_dict = {
    'O': [[0,0],[0,1],[0,2],[0,3],[0,4]],
    'P': [[0,0],[0,1],[0,2],[1,1],[1,2]],
    'Q': [[0,0],[0,1],[0,2],[0,3],[1,0]],
    'R': [[0,1],[1,1],[1,2],[2,2],[1,0]],
    'S': [[0,0],[0,1],[0,2],[1,2],[1,3]],
    'T': [[0,0],[1,0],[2,0],[1,1],[1,2]],
    'U': [[0,0],[1,0],[2,0],[0,1],[2,1]],
    'V': [[0,0],[1,0],[2,0],[0,1],[0,2]],
    'W': [[0,0],[0,1],[1,1],[1,2],[2,2]],
    'X': [[0,1],[1,1],[1,0],[1,2],[2,1]],
    'Y': [[0,0],[0,1],[0,2],[0,3],[1,2]],
    'Z': [[0,0],[1,0],[1,1],[1,2],[2,2]]
}
polyomino_dict['I'] = polyomino_dict['O']
polyomino_dict['L'] = polyomino_dict['V']
polyomino_dict['M'] = polyomino_dict['W']
polyomino_dict['F'] = polyomino_dict['R']
class Polyomino:
    def __init__(self,s,lattice='default'):
        if type(s)==str:
            s=polyomino_dict[s.upper()]
        if s==[]:
            self.n=0
            self.lattice='square'
        else:
            if lattice=='default' and type(s[0])==list:
                if len(s[0])==2:
                    lattice='square'
                elif len(s[0])==3:
                    lattice='triangle'
            elif isinstance(s[0],Point):
                lattice=dual(s[0].lattice)
            self.vertex_lattice = lattice
            self.face_lattice = dual(lattice)
            s=pointify(s,self.face_lattice)
            s=normalize(s)
            self.cells=sorted(s)
            self.n = len(s)
    def rotate_clockwise(self):
        return Polyomino([c.rotate() for c in self.cells],lattice=self.vertex_lattice)
    def flip(self):
        return Polyomino([c.flip() for c in self.cells],lattice=self.vertex_lattice)
    def __eq__(self,other):
        #are these polyominoes the same up to translation?
        #since cells are normalized, we just check equality
        c1=self.cells
        c2=other.cells
        if len(c1)!=len(c2):
            return False
        else:
            for i in range(len(c1)):
                if c1[i]!=c2[i]:
                    return False
            return True
    def scale(self,k):
        #TODO: implement for other lattices
        new_cells=[]
        for c in self.cells:
            for i in range(k):
                for j in range(k):
                    new_cells.append(Point([k*c.coords[0]+i,k*c.coords[1]+j]))
        return Polyomino(new_cells)
    def perimeter(self):
        #returns a list of Path objects
        pass
    def face_xrange(self):
        #range of x-values spanned by the faces of the shape
        if self.face_lattice=='square':
            return [0,max([c.coords[0] for c in self.cells])+1]
        elif self.face_lattice=='hexagon':
            mxs=[(c.coords[0]+c.coords[1]/2+c.coords[2]/2) for c in self.cells]
            return [min(mxs),max(mxs)+1]
        elif self.face_lattice=='triangle':
            mxs=[(c.coords[0]+c.coords[1]/2) for c in self.cells]
            return [min(mxs)-0.5,max(mxs)+0.5]
    def face_yrange(self):
        #range of y-values spanned by the faces of the shape
        if self.face_lattice=='square':
            return [0,max([c.coords[1] for c in self.cells])+1]
        elif self.face_lattice=='hexagon':
            mys=[0.866*c.coords[1] for c in self.cells]
            return [0,max(mys)+0.866]
        elif self.face_lattice=='triangle':
            mys=[0.866*c.coords[1] for c in self.cells]
            return [-0.57735,max(mys)+0.57735]
    def axis_spans(self):
        #how far of a range is there along each axis?
        #uses normalization
        result = [max([p.coords[i] for p in self.cells])+1 for i in range(len(self.cells[0].coords))]
        if self.face_lattice=='triangle':
            #this is a hack to make hexiamonds have the right border
            result[0]=max(result[0],max([p.coords[0]+p.coords[1] for p in self.cells]))
        if self.face_lattice=='hexagon':
            result[2]=2
        return result
    def draw(self,draw,translate,scale,fill='red',outline='white'):
        for p in self.cells:
            p.draw(draw,translate,scale,fill,outline)
    def render(self,img=None,N=200,M=-1,border=0.1,color='red'):
        if M==-1:
            M=N
        xrange = self.face_xrange()
        yrange = self.face_yrange()
        if img==None:
            img=Image.new('RGB',(N,M))
        draw = ImageDraw.Draw(img)
        s=img.size
        biggest_ratio=max((xrange[1]-xrange[0])/s[0],(yrange[1]-yrange[0])/s[1])
        scale = (1-2*border)/biggest_ratio
        translate = [border*s[0]-scale*xrange[0],s[1]+scale*yrange[0]-border*s[1]]
        self.draw(draw,translate,scale,color)
        return img
    def key(self):
        return tuple(self.cells)
    def __hash__(self):
        return hash(self.key())
        
#pentominoes = [Polyomino(s) for s in 'OPQRSTUVWXYZ']
# %%
def show(x):
    if isinstance(x,Polyomino):
        show(x.render())
    else:
        IPython.display.display(x)
# %%
def isometric_copies(p,reflection=True,only_translate=False):
    if only_translate:
        return [p]
    q=p.flip()
    copies=[]
    for i in range([1,2][reflection]):
        x=[p,q][i]
        if not x in copies:
            copies.append(x)
            for j in range(lattice_rotations(p.face_lattice)):
                x=x.rotate_clockwise()
                if not x in copies:
                    copies.append(x)
    return copies
# %%
def adjacent_cells(P):
    cells=P.cells
    adj = set()
    for c in cells:
        for c2 in c.neighbors():
            if not c2 in cells:
                adj.add(c2)
    return adj
def get_extensions(polys):
    new_polys = set()
    for poly in polys:
        for c in adjacent_cells(poly):
            new_cells = poly.cells+[c]
            new_poly = Polyomino(new_cells)
            add=True
            for variant in isometric_copies(new_poly):
                if variant in new_polys:
                    add=False
                    break
            if add:
                new_polys.add(new_poly)
    return new_polys
# %%
monominoes=[Polyomino([[0,0]])]
ominoes = [[],monominoes]
for i in range(2,11):
    ominoes.append(get_extensions(ominoes[-1]))
dominoes = ominoes[2]
trominoes = ominoes[3]
tetrominoes = ominoes[4]
pentominoes = ominoes[5]
hexominoes = ominoes[6]
heptominoes = ominoes[7]
octominoes = ominoes[8]
nonominoes = ominoes[9]
decominoes = ominoes[10]
print([len(x) for x in ominoes])
# %%
#optional: generate some more polyominoes, this takes a minute or two because the code is slow
for i in range(len(ominoes),12+1):
    print(i)
    ominoes.append(get_extensions(ominoes[-1]))
    print('\t',len(ominoes[-1]))
    print('\t',now())
# %%
def adjacency(c,cells):
    i,j=c
    coords = [(i,j+1),(i,j-1),(i-1,j),(i+1,j)]
    return tuple([int(c in cells) for c in coords])
 
def duplicates(P):
    cells = [tuple(c.coords) for c in P.cells]
    adjs = [adjacency(c,cells) for c in cells]
    return len(adjs)-len(set(adjs))
# %%
#some helper functions to make boards of arbitrary dimensions
#made to work on points for easier calling
def make_list(tup,elt=1):
    if len(tup)==1:
        return [elt]*tup[0]
    else:
        return [make_list(tup[1:],elt) for i in range(tup[0])]
def get_elt(l,path):
    if isinstance(path,Point):
        return get_elt(l,path.coords)
    if len(path)==1:
        return l[path[0]]
    else:
        return get_elt(l[path[0]],path[1:])
def set_elt(l,path,v):
    if isinstance(path,Point):
        set_elt(l,path.coords,v)
    else:
        if len(path)==1:
            l[path[0]]=v
        else:
            set_elt(l[path[0]],path[1:],v)
def full_copy(l):
    if type(l)==list:
        return [full_copy(e) for e in l]
    else:
        return l
def eltwise_max(l):
    maxes=l[0]
    for e in l[1:]:
        for i in range(len(e)):
            if e[i]>maxes[i]:
                maxes[i]=e[i]
    return maxes
# %%
class Count:
    #Class to store possible numbers of an object
    #currently very simple, just stores a range
    #but could modify to e.g. be all primes or something
    def __init__(self,low,high):
        self.low = low
        self.high = high
    def empty(self):
        return self.high<=0
    def remove(self,k=1):
        self.low-=k
        self.high-=k
    def copy(self):
        return Count(self.low,self.high)
    def sat(self,n=0):
        #does n satisfy the object's requirements?
        return self.low<=n<=self.high
    def __str__(self):
        return ("Count(%.1f,%.1f)" % (self.low,self.high))
class Collection:
    #Class for a collection of objects with restrictions on how many you can have
    #e.g "0-2 A's, and exactly 3 B's, and any nonnegative number of C's"
    def __init__(self,objects=None,counts=None,obj_map=None,costs=None,max_cost=None,style='standard'):
        self.style = style
        if style=='standard':
            self.n = len(objects)
            if counts is None:
                counts=[Count(0,float('inf')) for i in range(self.n)]
            if obj_map is None:
                obj_map = [i for i in range(self.n)]#just a mapping from objects to counts
            self.objects = objects
            self.counts = counts
            self.obj_map=obj_map
        elif self.style == 'costs':
            self.objects = objects
            self.points = max_cost
            self.costs = costs
    def remove(self,i):
        if self.style=='standard':
            self.counts[self.obj_map[i]].remove()
        elif self.style=='costs':
            self.points -= self.costs[i]
    def copy(self):
        if self.style=='standard':
            return Collection(self.objects,[c.copy() for c in self.counts],self.obj_map[:])
        elif self.style=='costs':
            return Collection(self.objects,costs=self.costs,max_cost=self.points,style='costs')
    def valid_indices(self):
        if self.style=='standard':
            return [i for i in range(self.n) if not self.counts[self.obj_map[i]].empty()]
        elif self.style=='costs':
            return [i for i in range(len(self.costs)) if self.costs[i]<=self.points]

# %%
def backtracking_cover(board,neighbors,targets,shapes,method,depth=0):
    #will attempt to cover targets via backtracking
    #board: cells, 1 if open and 0 if restricted. Indexed by 
    #targets: a list of cells, in order of covering priority
    #shapes: a Collection of allowable polyominoes via translation
    #returns a list of (list of cell)s covering the shape, or False if no such exists
    
    #Method determines the search process
    #standard: search among targets in the order provided
    #neighbors: prioritize those cells in target with the fewest permissible neighbors
    #random: at every stage, select a random element of target
    if len(targets)==0:
        return []
    t=None
    if method=='standard':
        t=targets[0]
    elif method=='neighbors':
        best_t = None
        fewest_neighbors=float("inf")
        for t in targets:
            nn=get_elt(neighbors,t)
            if nn<fewest_neighbors:
                fewest_neighbors = nn
                best_t=t
            if fewest_neighbors==0: #save time
                break
        t=best_t
    elif method=='random':
        t=random.choice(targets)
    #print('\t'*depth+"t=",t)
    #print("backtracking_cover called")
    #print([[c.low,c.high] for c in shapes.counts])
    #print(shapes.valid_indices())
    for i in shapes.valid_indices():
        poly = shapes.objects[i] #get the corresponding object from the collection
        #put shape S on t somehow
        s = poly.cells
        for c in s:
            #translate so cell c covers t
            vec = t-c #might fail, only proceed if not
            if vec:
                new_s = translate(s,vec)
                attempt_works = True #does the new shape fit?
                for c2 in new_s:
                    if not get_elt(board,c2):
                        attempt_works=False
                        break
                if attempt_works: #shape fits
                    new_targets = [c2 for c2 in targets if not(c2 in new_s)]
                    new_board=full_copy(board)
                    new_neighbors=full_copy(neighbors)
                    for c2 in new_s:
                        set_elt(new_board,c2,0)
                        if method=='neighbors': #save time if we're not consulting the array
                            for c3 in c2.neighbors():
                                one_lower = get_elt(new_neighbors,c3)-1 #decrement neighbor count by 1
                                set_elt(new_neighbors,c3,one_lower)
                    new_shapes = shapes.copy()
                    new_shapes.remove(i)
                    cover_attempt = backtracking_cover(new_board, new_neighbors, new_targets,new_shapes,method,depth+1)
                    if cover_attempt!=False:
                        return [new_s]+cover_attempt
    return False
def get_border(s):
    #find maximum x/y span of a polyomino
    return max([max(c) for c in s.cells])+1
def tile(target,shapes,counts=None,reflection=True,tile_costs=None,max_tile_cost=0,only_translate=False,give_args=False,method='standard'):
    return cover(target,shapes,counts,reflection,True,tile_costs,max_tile_cost,only_translate,give_args,method)
def cover(target,shapes,counts=None,reflection=True,tile=False,tile_costs=None,max_tile_cost=0,only_translate=False,give_args=False,method='standard'):
    #Determine if the target can be covered by nonoverlapping copies of the shape
    if isinstance(shapes,Polyomino):
        shapes = [shapes]
    if isinstance(target,Polyomino):
        target = target.cells
    if isinstance(target,str):
        target = str_to_cells(target)
    if len(shapes)==1 and tile:
        if len(target)%shapes[0].n !=0: #make sure the number of cells matches
            return False
    lattice = shapes[0].face_lattice
    target = pointify(target,lattice)

    border = max([max(s.axis_spans()) for s in shapes])
    
    if not isinstance(target[0],Point):
        target = [Point(e) for e in target]
    t = normalize(target) #give things nonnegative coordinates
    t_axes = Polyomino(t,lattice).axis_spans()
    
    if method=='taxicab': #sort points by sum of coords
        target = sorted(target,key=(lambda p: sum(p.coords)))
        method='standard'

    v=(1-tile) #value to fill board with by default
    board=make_list([t_axes[i]+2*border*(i!=2) for i in range(len(t_axes))],v)
    
    neighbor_v = [neighbor_count(lattice),0][tile] #default number of neighbors
    neighbors = make_list([t_axes[i]+2*border*(i!=2) for i in range(len(t_axes))],neighbor_v)
    
    k=2 #TODO: maake k vary based on the lattice (happens to be 2 for everything so far)
    translation = Point([border]*k,translation_lattice(lattice))
    t = translate(t,translation) #now we have a buffer zone around the target
    
    if tile: #board is zeroed, set the target cells to 1
        for p in t:
            set_elt(board,p,1)
        #after board is configured, specify neighbors
        for p in t:
            num_neighbors=sum([get_elt(board,q) for q in p.neighbors()])
            set_elt(neighbors,p,num_neighbors)
    #print("board,target,shapes",len(board),t,shapes)
    full_shapes = []
    if counts is None:
        counts = [Count(0,float('inf')) for s in shapes] #full range
    obj_map = [] #map from shapes to counts
    full_costs=[]
    for i in range(len(shapes)):
        s = shapes[i]
        iso_copies = isometric_copies(s,reflection,only_translate)
        full_shapes+=iso_copies
        obj_map+=[i]*len(iso_copies)
        if tile_costs is not None:
            full_costs+=[tile_costs[i]]*len(iso_copies)
    if tile_costs is None:
        full_shape_collection = Collection(full_shapes,counts,obj_map)
    else:
        full_shape_collection = Collection(full_shapes,costs=full_costs,max_cost = max_tile_cost,style='costs')
    full_cells = [s.cells for s in full_shapes] #our goals up to translation
    if method!='neighbors':
        neighbors=[]#save time and memory in calls to backtracking_cover
    if give_args:
        return (board,neighbors,t,full_shape_collection,method)
    sol = backtracking_cover(board,neighbors,t,full_shape_collection,method)
    if sol==False:
        return False
    else:
        #translate solution back to match original tiling shape, so the buffer doesn't mess things up
        return [translate(bit,-translation) for bit in sol]
#todo: add subtile_board, see if that helps
def cover_board(Q,P,verbose=False,N=200,only_translate=False,method='standard'):
    #Cover a board Q with polyominoes P
    #Q is given as a multi-line string of 2s (must cover), 1s (may cover) and 0s (must not cover)
    #then is padded with 1s as need be
    if type(Q)==str:
        Q = Q.split('\n')
    if type(P)==Polyomino:
        P=[P]
    lattice = P[0].face_lattice
    board = {} #dict coord -> value
    for i in range(len(Q)):
        for j in range(len(Q[i])):
            board[(j,-i)]=int(Q[i][j])

    two_cells=[Point(coord,lattice) for coord in board if board[coord]==2]
    zero_cells = [coord for coord in board if board[coord]==0]
    bc = min(two_cells)
    args = cover(two_cells,P,only_translate=only_translate,method=method,give_args = True)
    #args[0] is the board
    #args[2] are the targets
    p = (min(args[2])-bc).coords #translation vector
    for z in zero_cells:
        args[0][p[0]+z[0]][p[1]+z[1]]=0
    attempt = backtracking_cover(args[0],args[1],args[2],args[3],args[4])

    
    if verbose:
        show(tiling_render(attempt,N=N,scale=15,border=0))
        return attempt
    else:
        return attempt
# %%
def tiling_render(tiles,scale,border=0.1,N=200,M=-1,target=[],k=0,r=0.2):
    if tiles==False:
        return "Error: no tiling!"
    lattice = tiles[0][0].lattice
    if M==-1:
        M=N
    img = Image.new('RGB',(N,M))
    draw = ImageDraw.Draw(img)
    s=img.size
    min_dim = min(s)
    for t in tiles:
        color=(random.randint(0,255),random.randint(0,255),random.randint(0,255))
        translate = [border*min_dim,s[1]-border*min_dim]
        for c in t:
            c.draw(draw,translate,scale,fill=color,outline='white')
    for c in target:
        color=(255,255,255)
        coords=[c[0]+0.5-r,-c[1]-0.5-r,c[0]+0.5+r,-c[1]-0.5+r]
        translate = [k*scale+border*s[0],s[1]-border*s[1]-k*scale]
        real_coords = [coords[i]*scale+translate[i%2] for i in range(4)]
        #print(real_coords)
        draw.ellipse(real_coords,fill=color)
    return img
# %%
def str_to_cells(q):
    q=q.split('\n')
    cells=[]
    for i in range(len(q)):
        for j in range(len(q[i])):
            if q[i][j]=='1':
                cells.append([j,-i])
    return cells
def str_to_poly(q,shape='square'):
    return Polyomino(str_to_cells(q),lattice=shape)
# %%
#example to demonstrate functionality for checking candidate polyominoes

#this board describes a tiling condition, where 2s must be covered, 1s may be covered, and 0s may not be covered
B = """
121
20
121
"""
P = str_to_poly("""
...1
...1
1111111
1.....1
1.....1
1.....1
1.....1
1111111
...1
...1
""")
P = str_to_poly("""
..1111
..1..1
111...11
1......1
1......1
111..111
..1..1
..1111
""")   


P = str_to_poly("""
..11111111
..11111111
......1.....
11....1...11
1.....1....1
1.....1....1
1.....1....1
111111111111
1.....1....1
1.....1....1
1.....1....1
11....1...11
......1.....
..11111111
..11111111
""")   

print(P.n)            
for P in [P]:
    c = cover_board(B,P)
    show(tiling_render(c,N=500,scale=15))