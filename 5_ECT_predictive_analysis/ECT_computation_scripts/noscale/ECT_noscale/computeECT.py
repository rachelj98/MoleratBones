# This computes ECT for a given OFF file
# Optional array of directions (coordinate vectors on S^2 as a subset of R^3:
# A row per direction:
# x_1 y_1 z_1
# x_2 y_2 z_2
#
# Returns the EC map parameters for the given directions (I.e. a device for computing exact EC for the given directions)
import numpy as np
from sets import Set
from itertools import combinations
#import PHTzero
import operator

def Euler_Curve(vertex_height_dict, edge_heights, face_heights):
    sort = sorted(vertex_height_dict.items(), key=operator.itemgetter(1))
    V=0
    E=0
    F=0
    levelset=[]
    EC_old=0
    EC_dict=dict()
    for pair in sort:
        levelset.append(pair[0])
        E=sum(i<=pair[1] for i in edge_heights)
        F=sum(i<=pair[1] for i in face_heights)
        #sub_edges=Set(list(combinations(levelset,2))) &edges
        #sub_faces=Set(list(combinations(levelset,3))) &faceSet
        EC=len(levelset)-E+F#-len(sub_edges)#+len(sub_faces)
        if (EC != EC_old):
            EC_dict[pair[1]]=EC
            EC_old=EC
    return EC_dict

def computeDiameter(filename):
    file=open(filename,"r")
    A=file.readline().split()
    (V,F,E)=map(int,file.readline().strip().split(' '))
    vertices=np.empty([V,3])
    #faces=np.empty([F,3])
    # Read in the vertices
    for i in range(0,V):
        vertices[i]=map(float,file.readline().strip().split(' '))
    V=np.transpose(vertices)
    xscale=max(V[0])-min(V[0])
    yscale=max(V[1])-min(V[1])
    zscale=max(V[2])-min(V[2])
    scale=max(xscale,yscale,zscale)
    return(scale)

def preparecomplex(filename):
    file=open(filename,"r")
    # Checking we have valid headers
    A=file.readline().split()
        #if(A != 'OFF')
        #TODO error handling    
    #Reading in the number of vertices, faces and edges, and formatting their arrays
    (V,F,E)=map(int,file.readline().strip().split(' '))
    vertices=np.empty([V,3])
    faces=np.empty([F,3])
    # Read in the vertices
    for i in range(0,V):
        vertices[i]=map(float,file.readline().strip().split(' '))
    #Read in the faces
    for i in range(0,F):
        line=map(int,file.readline().strip().split(' '))
        #TODO Sanity check: Insist we have a triangulated mesh
        #if line[0]!=3
        # Warning: the mesh contains non-triangular faces, holes might have been created
        faces[i]=line[1:4]
    tvertices=np.transpose(vertices)
    tvertices[0]=tvertices[0]-np.mean(tvertices[0])
    tvertices[1]=tvertices[1]-np.mean(tvertices[1])
    tvertices[2]=tvertices[2]-np.mean(tvertices[2])
    vertices=np.transpose(tvertices)

    #TODO Creating set of edges from the faces

    edges=Set()
    faceSet=Set()
    for i in faces:
        edge=combinations(i,2)
        edges.update(edge)
        face=combinations(i,3)
        faceSet.update(face)
    list_of_edges=list(edges)
    list_of_faces=list(faceSet)

    ordered_list_of_edges=list()
    for t in list_of_edges:
        r=min(t)
        R=max(t)
        tu=(r,R)
        ordered_list_of_edges.append(tu)
    edges=set(ordered_list_of_edges)
    list_of_edges=list(edges)
    dict_of_vertices=dict()    
    k=0
    for i in vertices:
        dict_of_vertices.update({k:tuple(i)})
        k=k+1
    vertex_dict=dict()
    for j in range(len(vertices)):
        vertex_dict[j]=vertices[j]
    return(vertex_dict,list_of_edges,list_of_faces)

# Takes in: output from preparecomplex, and a direction
# Returns: the EC curve parameters that can be fed to evaluate_EC
def computeECT(vertex_dict,list_of_edges,list_of_faces,direction):
    vertex_height_dict=dict()
    for vertex in vertex_dict:
        height=np.dot(direction,vertex_dict[vertex])
        vertex_height_dict[vertex]=height

    #edge_height_dict=dict()
    edge_heights=list()
    for edge in list_of_edges:
        h1=vertex_height_dict[edge[0]]
        h2=vertex_height_dict[edge[1]]
        height=max(h1,h2)
        #edge_height_dict[edge]=height
        edge_heights.append(height)

    #face_height_dict=dict()
    face_heights=list()
    for face in list_of_faces:
        h1=vertex_height_dict[face[0]]
        h2=vertex_height_dict[face[1]]
        h3=vertex_height_dict[face[2]]
        height=max(h1,h2,h3)
        #face_height_dict[face]=height
        face_heights.append(height)
    EC=Euler_Curve(vertex_height_dict, edge_heights, face_heights)
    return(EC)



# A function for evaluating an Euler Curve at x 
def evaluate_EC(EC,x):
    if x< min(EC):
        return(0)
    threshold=max(element for element in EC.keys() if element <= x)
    return(EC[threshold])
