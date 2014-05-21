from igraph import *
import networkx as nx
import numpy as np
import random as r
import math
import operator
import pprint
from scipy.cluster.vq import vq, kmeans, whiten

MAX_ITERATIONS=1000

def shouldStop(oldCentroids, centroids,iterations):
	if iterations > MAX_ITERATIONS:
		return True
	return (oldCentroids == centroids).all()

def get_labels(T,centroids):
	labels=np.array([min(enumerate(np.sqrt(((T[i]-centroids)**2).sum(axis=0))), key=operator.itemgetter(1))[0] for i in range(T.shape[0])])
	return labels

def get_centroids(T,labels,k):
	indices=np.array([[i for i, x in enumerate(labels) if x == j]for j in range(k)])
	list_values=np.array([[T[j] for j in indices[i]]for i in range(k)])
	centroids=np.array([np.mean(list_values[i],axis=0) for i in range(k)])
	return centroids
	

def k_means(T, k):
	centroids=centroids = r.sample(T, k)
	iterations=0
	oldCentroids=np.empty(shape=(k,k))
	while not (shouldStop(oldCentroids,centroids,iterations)):
		oldCentroids = centroids
		iterations += 1
		#print iterations
		labels=get_labels(T,centroids)
		centroids=get_centroids(T,labels,k)
	#print centroids	
	return centroids



#input similarity matrix

G=nx.read_gml("astro-ph.gml") 
A=nx.to_numpy_matrix(G)
#print A
degree_dict=G.degree()
degree_list_values=degree_dict.values()
degree_list_values_1_2=np.power(degree_list_values,-0.5)
D=np.diag(degree_list_values)
L=D-A

D_1_2=np.diag(degree_list_values_1_2)
L_sym=D_1_2*L*D_1_2

k = int(raw_input("Enter the number of clusters: "))


eig_vals, eig_vecs = np.linalg.eig(L_sym)
eig_vals_sorted = np.sort(eig_vals)
eig_vecs_sorted = eig_vecs[eig_vals.argsort()]

#print eig_vecs_sorted

k_eig_vecs_sorted=eig_vecs_sorted[:k,:]

U=k_eig_vecs_sorted.T

#print U
#print U.shape

norms=np.apply_along_axis(np.linalg.norm, 1, U)

#print "norms"
#print norms

T = np.array([[U.item(i,j)/norms[i] for j in range(U.shape[1])] for i in range(U.shape[0])])
#print T

centroids=k_means(T,k)
print centroids


#whitened = whiten(T)
book,distortion=kmeans(T,k)
print book
labels=get_labels(T,centroids)
labels1=get_labels(T,book)
clusters=np.array([[i for i, x in enumerate(labels) if x == j]for j in range(k)])
#pprint.pprint(clusters,width=1)
for i in range(k):
	print "Cluster "+str(i)
	for j in range(len(clusters[i])):
		print clusters[i][j]+1
	
clusters1=np.array([[i for i, x in enumerate(labels1) if x == j]for j in range(k)])

for i in range(k):
	print "Cluster "+str(i)
	for j in range(len(clusters1[i])):
		print clusters1[i][j]+1


#print clusters1
#pprint.pprint(clusters,width=1)



