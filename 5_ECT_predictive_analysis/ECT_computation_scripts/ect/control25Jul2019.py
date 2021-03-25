# This is the orchestration script
import computeECT
import os
import math
import numpy


def normalize(v):
    s=1.0/math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) 
    w= (v[0]*s, v[1]*s, v[2]*s)
    return w
def find_directions():
    #   Since we are performing the PHT for a simplicial complex in R^3 we need to have a set of directions
    #   to approximate S^2.
    #   We are effectivly using a barycentric subdivision of an icosahedron (the plationic solid with 20 faces).
    phi=(1 + math.sqrt(5))/2
    list_of_keys=[]

    icos_vertices ={1:(0,1, phi), 2:(0,1,-phi), 3:(0,-1,phi), 4:(0,-1,-phi),5:(1,phi, 0),
    6:(1,-phi,0), 7:(-1,phi,0), 8:(-1,-phi,0), 9:(phi,0,1), 10:(phi,0,-1), 
    11:(-phi,0,1), 12:(-phi,0,-1)}
    # These are the 12 vertices of an icosahedron. We will normalize later to get unit vectors.


    dict_of_adj={1:set(),2:set(),3:set(),4:set(),5:set(),6:set(),7:set(),
    8:set(),9:set(),10:set(),11:set(),12:set()}

    dict_of_faces={}

    # vertices from first subdivision - here we split each edge. 
    # vertices form second subdivision we further spit the edges and have three in each face

    sub1 ={}
    #The midpoint of each of the edges in the icosahedron

    sub2edge={}
    #The points point the vertices of the icosahedron and the midpoints in the edges of the icosahedron

    sub2face={}
    #The points in the face (we have 3 in each face?)



    length1=math.sqrt(1/(phi*phi) +1)/2
    length2=math.sqrt(1/(4*phi*phi) +1)/4

    normalization=1.0/math.sqrt(math.pow(phi,2) + 1)

    for i in range(1,13):
        for j in range(i+1,13):
            a=icos_vertices[i]
            b=icos_vertices[j]
            if (a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2])==4:

                v1=(length2*(3*a[0]+b[0]),length2*(3*a[1]+b[1]),length2*(3*a[2]+b[2]))
                v2=(length1*(a[0]+b[0]),length1*(a[1]+b[1]),length1*(a[2]+b[2]))
                #			print math.pow(normalization*length1*(a[0]+b[0]),2) + math.pow(normalization*length1*(a[1]+b[1]),2)+math.pow(normalization*length1*(a[2]+b[2]),2) 
                v3=(length2*(a[0]+3*b[0]),length2*(a[1]+3*b[1]),length2*(a[2]+3*b[2]))
                sub1.update({(i,i,j,j):v2})

                list_of_keys.append((i,i,j,j))
                list_of_keys.append((i,i,i,j))
                list_of_keys.append((j,j,j,i))

                sub2edge.update({(i,i,i,j):v1, (j,j,j,i):v3})

                dict_of_adj[i].add(j)
                dict_of_adj[j].add(i)

    a=phi*phi+1
    b=8*phi*phi + 8*phi+4
    length3=math.sqrt(a/b)

    for i in range(1,13): 	
        for j in dict_of_adj[i]:
            a=icos_vertices[i]
            b=icos_vertices[j]
            for k in dict_of_adj[i] & dict_of_adj[j]:
                c=icos_vertices[k]
                det=a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
                if det>0:	
                    dict_of_faces.update({(i,j):k})
                    v=(length3*(2*a[0]+b[0] + c[0]),length3*(2*a[1]+b[1]+c[1]),length3*(2*a[2]+b[2]+c[2]))
                    sub2face.update({(i,i,j,k):v})
                    list_of_keys.append((i,i,j,k))

    dict_of_directions=dict(sub2edge.items() + sub1.items()  + sub2face.items())
    #  	dict_of_directions=dict()

    for v in icos_vertices:
        dict_of_directions.update({(v,v,v,v):icos_vertices[v]})
        list_of_keys.append((v,v,v,v))


    for key in list_of_keys:
        dict_of_directions[key]=normalize(dict_of_directions[key])
    # normalises so vectors of unit length.

    return (list_of_keys, dict_of_directions)

# Get directions ala Kate
keys,dict_of_dirs=find_directions()
directions=dict_of_dirs.values()


#Get the function handle (the expensive step!)
cwd=os.getcwd()
path=cwd+'/feminput/'
filename='FILEINFO_FEM_100K_aligned_aligned.off'
V,E,F=computeECT.preparecomplex(path+filename)


#parallelization goes here
ECS=list()
params=list()
for t in directions:
    params.append((V,E,F,t))

def multi_run_wrapper(args):
   return computeECT.computeECT(*args)
if __name__ == "__main__":
    from multiprocessing import Pool
    pool = Pool(48)
    ECS = pool.map(multi_run_wrapper,params)
    #print results
#for t in directions:
   #EC=computeECT.computeECT(V,E,F,t)
   #ECS.append(EC)

#for t in os.listdir(path):

# Evaluating the function handle (A quick script for turning the above function handle to a comma separated matrix
outpath=cwd+'/output25Jul2019_fem/'
nsublevelset=5000
xset=[3.0**(0.5)*(2*j/float(nsublevelset-1)-1) for j in range(nsublevelset)]
file =open(outpath+filename,'w')
for j in range(len(directions)):
    EC=list()
    for t in xset:
        EC.append(computeECT.evaluate_EC(ECS[j],t))
    for k in EC:
        file.write(str(k)+',')
    file.write('\n')
file.close()


