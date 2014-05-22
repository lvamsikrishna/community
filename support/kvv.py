from igraph import *
import networkx as nx
import numpy as np
import random as r
import scipy.sparse.linalg

def cut(A,groups):	
	cut_size=0
	for i in range(A.shape[0]):
		for j in range(A.shape[1]):
			if(A.item(i,j)==1 and groups[i]!=groups[j]):
				cut_size=cut_size+1
	#print "cut_size"
	#print cut_size
	return cut_size

def vol(D,groups):
	volA=0
	volB=0
	for i in range(groups.shape[0]):
		if(groups[i]==0):
			volA=volA+D.item(i,i)
		else:
			volB=volB+D.item(i,i)
	#print "volA"
	#print volA
	#print "volB"
	#print volB
	return (volA,volB)

def ncut(A,D,groups):
	vols=vol(D,groups)
	volA=vols[0]
	volB=vols[1]
	n_cut=cut(A,groups)*((1.0/volA)+(1.0/volB))
	return n_cut
	
	


G=nx.read_gml("karate.gml") 
A=nx.to_numpy_matrix(G)
N=A.shape[0]
degree_dict=G.degree()
degree_list_values=degree_dict.values()
D=np.diag(degree_list_values)
Dinv=np.linalg.inv(D)
P=Dinv*A
w,v=scipy.sparse.linalg.eigs(P,k=2)
v_2=v[:,1]

k = int(raw_input("Enter the number of clusters: "))

v_2_sorted=np.sort(v_2)
indices=np.argsort(v_2)
groups=np.zeros(N)
ncuts=np.empty(N-1)
minimum_ncut=np.inf
finalGrouping=np.empty(N)
for i in range(N-1):	
	groups[indices[i]]=1
	#print "i  "+str(i)
	#print groups
	ncuts[i]=ncut(A,D,groups)
	if(ncuts[i]<minimum_ncut):
		minimum_ncut=ncuts[i]
		finalGrouping=groups.copy()

mask1=finalGrouping==0
mask2=finalGrouping==1
A_1=A[mask1][:,mask1]
A_2=A[mask2][:,mask2]
D_1=D[mask1][:,mask1]
D_2=D[mask2][:,mask2]
D_inv_1=np.linalg.inv(D_1)
D_inv_2=np.linalg.inv(D_2)
P_1=D_inv_1*A_1
P_2=D_inv_2*A_2
w1,v1=scipy.sparse.linalg.eigs(P_1,k=2)
w2,v2=scipy.sparse.linalg.eigs(P_2,k=2)
if(w1[1]<w2[1]):
	v_2=v1[:,1]
else:
	v_2=v2[:,1]

#print A_new.shape
'''print ncuts
print finalGrouping
print minimum_ncut
print min(ncuts)
print np.argmin(ncuts)
'''



