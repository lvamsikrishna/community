from igraph import *
import networkx as nx
import numpy as np
import random as r
import scipy.sparse.linalg

GROUPS_FORMED=2

def cut(A,groups,formed):	
	cut_size=0
	for i in range(A.shape[0]):
		for j in range(A.shape[1]):
			if(A.item(i,j)==1 and groups[i]!=groups[j]):
				cut_size=cut_size+1
	#print "cut_size"
	#print cut_size
	return cut_size

def vol(D,groups,formed):
	volA=0
	volB=0
	for i in range(groups.shape[0]):
		if(groups[i]==formed):
			volA=volA+D.item(i,i)
		else:
			volB=volB+D.item(i,i)
	#print "volA"
	#print volA
	#print "volB"
	#print volB
	return (volA,volB)

def ncut(A,D,groups,formed):
	vols=vol(D,groups,formed)
	volA=vols[0]
	volB=vols[1]
	n_cut=cut(A,groups,formed)*((1.0/volA)+(1.0/volB))
	return n_cut
	
	
def SM(A,D,k,formed,grpNumber):
	N=A.shape[0]
	Dinv=np.linalg.inv(D)
	P=Dinv*A
	w,v=scipy.sparse.linalg.eigs(P,k=2)
	v_2=v[:,1]
	v_2_sorted=np.sort(v_2)
	indices=np.argsort(v_2)
	groups=np.zeros(N)
	groups.fill(grpNumber)
	ncuts=np.empty(N-1)
	minimum_ncut=np.inf
	finalGrouping=np.empty(N)
	for i in range(N-1):	
		groups[indices[i]]=formed
		ncuts[i]=ncut(A,D,groups,formed)
		if(ncuts[i]<minimum_ncut):
			minimum_ncut=ncuts[i]
			finalGrouping=groups.copy()
	#print finalGrouping
	mask1=finalGrouping==grpNumber
	mask2=finalGrouping==formed
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
	mask=np.empty(N)
	if(w1[1]<w2[1]):
		A=A_2
		D=D_2
		mask=mask2
		grpNumber=formed
	else:
		A=A_1
		D=D_1
		mask=mask1
	indices=np.where(mask==True)
	#print indices
	if(formed<k):
		ret=SM(A,D,k,formed+1,grpNumber)
		indices=np.where(mask==True)
		indices=indices[0]
		for i in range(indices.shape[0]):
			finalGrouping[indices[i]]=ret[i]	
	return finalGrouping


G=nx.read_gml("karate.gml") 
A=nx.to_numpy_matrix(G)
N=A.shape[0]
degree_dict=G.degree()
degree_list_values=degree_dict.values()
D=np.diag(degree_list_values)
k = int(raw_input("Enter the number of clusters: "))
result=SM(A,D,k,2,1)
print result

#print A_new.shape
'''print ncuts
print finalGrouping
print minimum_ncut
print min(ncuts)
print np.argmin(ncuts)
'''



