import networkx as nx
import numpy as np
import random as r
import math
import operator
import pprint

MAX_ITERATIONS=1000

def shouldStop(oldCentroids, centroids,iterations):
	if iterations > MAX_ITERATIONS:
		return True
	return (centroids==oldCentroids).all()

def get_labels(T,centroids):
	labels=np.array([min(enumerate(np.sqrt(((T[i]-centroids)**2).sum(axis=1))), key=operator.itemgetter(1))[0] for i in range(T.shape[0])])
	return labels

def get_centroids(T,labels,k):
	indices=np.array([[i for i, x in enumerate(labels) if x == j]for j in range(k)])
	list_values=np.array([[T[j] for j in indices[i]]for i in range(k)])
	centroids=np.array([np.mean(list_values[i],axis=0) for i in range(k)])
	return centroids
	

def k_means(T, k):
	centroids= np.array(r.sample(T, k))
	print centroids
	iterations=0
	oldCentroids=np.empty(shape=(k,k))
	while not (shouldStop(oldCentroids,centroids,iterations)):
		oldCentroids = centroids
		iterations += 1
		labels=get_labels(T,centroids)
		centroids=get_centroids(T,labels,k)
	return centroids


def NJW(A,D,D_1_2,k):
	L=D-A
	L_sym=D_1_2*L*D_1_2
	eig_vals, eig_vecs = np.linalg.eig(L_sym)
	eig_vals_sorted = np.sort(eig_vals)
	eig_vecs_sorted = eig_vecs[eig_vals.argsort()]
	U=eig_vecs_sorted[:,:k]
	norms=np.apply_along_axis(np.linalg.norm, 1, U)
	T = np.array([[U.item(i,j)/norms[i] for j in range(U.shape[1])] for i in range(U.shape[0])])
	centroids=k_means(T,k)
	labels=get_labels(T,centroids)
	clusters=np.array([[i for i, x in enumerate(labels) if x == j]for j in range(k)])
	return (centroids,clusters)



G=nx.read_gml("karate.gml") 
A=nx.to_numpy_matrix(G)
degree_dict=G.degree()
degree_list_values=degree_dict.values()
degree_list_values_1_2=np.power(degree_list_values,-0.5)
D=np.diag(degree_list_values)
D_1_2=np.diag(degree_list_values_1_2)
k = int(raw_input("Enter the number of clusters: "))
centroids,clusters=NJW(A,D,D_1_2,k)
print "results"
print centroids
print clusters

