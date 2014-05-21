from igraph import *
import networkx as nx
import numpy as np
import random as r
import scipy.sparse.linalg

def cut(Adjacency,c1,c2):
	cut=0
	return cut

def vol(c):
	return 0

def ncut(Adjacency,c1,c2):
	ncut=cut(Adjacency,c1,c2)*((1/vol(c1))+(1/vol(c2)))
	return ncut
	
	


G=nx.read_gml("karate.gml")
 
A=nx.to_numpy_matrix(G)
print A
N=A.shape[0]
degree_dict=G.degree()
degree_list_values=degree_dict.values()
D=np.diag(degree_list_values)
Dinv=np.linalg.inv(D)
P=Dinv*A
w,v=scipy.sparse.linalg.eigs(P,k=2)
v_2=v[:,1]
v_2_sorted=np.sort(v_2)


