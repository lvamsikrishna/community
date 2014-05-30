import networkx as nx
import numpy as np
import random as r
import math
import operator
import pprint

#sets the maximum number of iterations for the k_means algorithm
MAX_ITERATIONS=1000

'''
Checks for the stopping conditon of the k_means algorithm 
k-means algorithm stops either when the number of iterations has exceeded
the maximum number of iterations
or 
When centroids in two consecutive iterations do not change.
input:
oldCentroids -	centroids in previous iteration
centroids	 -	centroids in current itearation
iterations   - 	current number of iterations
returns:
true - if iterations should Stop
false - otherwise
'''


def shouldStop(oldCentroids, centroids,iterations):
	if iterations > MAX_ITERATIONS:
		return True
	return (centroids==oldCentroids).all()

'''
Assigns a group label to each point. Given centroids, 
the cluster label is the cluster number for which the 
eucledian distance from that point to the cluster centroid
is minimum.
input:
T - Array of points for which labels are to be assigned
centroids - the centroids of the clusters
returns:
labels - group labels corresponding to each point to be clustered
'''


def get_labels(T,centroids):
	labels=np.array([min(enumerate(((T[i]-centroids)**2).sum(axis=1)), key=operator.itemgetter(1))[0] for i in range(T.shape[0])])
	return labels


'''
Given group labels for each point(node) clusters finds the centroid for each cluster. Centroid of a particular group or cluster is the average of all the points in that cluster.
Input:
T - array of points
labels - group label for each point 
k - number of clusters to be formed
returns:
centroids - an array of centroids, each element corresponding to the centroid of a cluster
'''


def get_centroids(T,labels,k):
	indices=np.array([[i for i, x in enumerate(labels) if x == j]for j in range(k)])
	list_values=np.array([[T[j] for j in indices[i]]for i in range(k)])
	centroids=np.array([np.mean(list_values[i],axis=0) for i in range(k)])
	return centroids
	

'''
Clusters the points in array T into k clusters using k means clustering algorithm.
input:
T - array of points to be clustered
k - number of clusters to be formed.
returns:
centroids - centroids of the k clusters formed
'''

def k_means(T, k):
	np.random.shuffle(T)
	centroids=T[:k,:]
	iterations=0
	oldCentroids=np.empty(shape=(k,k))
	while not (shouldStop(oldCentroids,centroids,iterations)):
		oldCentroids = centroids
		iterations += 1
		labels=get_labels(T,centroids)
		centroids=get_centroids(T,labels,k)
	return centroids

'''
Clustering the nodes of network represented by adjacency matrix A into k clusters.
input:
A -  	Adjacency matrix
D -  	degree matrix, represented as a diagonal matrix
k -		number of clusters to be formed
'''

def MS(A,D,k):
	L=np.subtract(D,A)
	L_rw=np.linalg.inv(D)*L
	eig_vals, eig_vecs = np.linalg.eig(L_rw)
	eig_vals_sorted = np.sort(eig_vals)
	eig_vecs_sorted = eig_vecs[eig_vals.argsort()]
	T=eig_vecs_sorted[:,:k]
	T=np.around(T,decimals=12)
	centroids=k_means(T,k)
	labels=get_labels(T,centroids)
	clusters=np.array([[i for i, x in enumerate(labels) if x == j]for j in range(k)])
	return (centroids,clusters)


#reading the dataset from the file in gml format
G=nx.read_gml("football.gml") 

A=nx.to_numpy_matrix(G)
#obtaining the dictionary of (node_key : degree) of the graph
degree_dict=G.degree()
#obtaining the list of degrees only
degree_list_values=degree_dict.values()
#D diagonal matrix, diagonal elements representing the degrees of the nodes
D=np.diag(degree_list_values)
k = int(raw_input("Enter the number of clusters: "))
centroids,clusters=MS(A,D,k)
print "results"
print centroids
print clusters

