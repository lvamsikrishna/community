from igraph import *
import networkx as nx
import numpy as np
import random as r
import math
import pylab
import scipy.sparse.linalg

def is_positive(x):
	return x>(1*pow(10,-6))

def power_method(A, x, maxit):
	"""Does maxit iterations of the power method
	on the matrix mat starting at start.
	Returns an approximation of the largest
	eigenvector of the matrix mat."""
	eigenvalue=0.0
	tolerence = 1 * pow(10,-9)
	for i in xrange(maxit):
		oldx = x 
		y = A*x
		x=y/np.linalg.norm(y)
		oldeigenvalue = eigenvalue 
		eigenvalue = np.linalg.norm(A*x) 
		if abs(eigenvalue - oldeigenvalue) < tolerence:
			break
	if i==maxit:
		return []
	else:
		return x 
		

def compute_eigenvalue(A,x): 
	Axx=x.T*A*x
	xx=x.T*x
	print xx
	result= Axx/xx
	return result

def delta(i,j,groups):
	#print groups[i]==groups[j]
	if(groups[i]==groups[j]):
		return 1
	else:
		return 0


def recalculate_lead_eigen(modularity_mat,eigenvalue):
	I=np.identity(modularity_mat.shape[0])
	myList=[r.randint(-1,1) for i in range(modularity_mat.shape[0])]
	start=np.matrix(myList)
	start=start.T
	eI=eigenvalue*I	
	modularity_mat_new= np.array([[modularity_mat.item((i,j))-eI.item((i,j)) for j in range(modularity_mat.shape[1])] for i in range(modularity_mat.shape[0])])
	result=power_method(modularity_mat_new,start,10000)
	#print result.T

	largest_eigen_value=compute_eigenvalue(modularity_mat_new,result)
	#print largest_eigen_value+eigenvalue
	return result


def leadEigenVectorMethod(G,B,k,two_m,groups,group_number):
	print G.nodes()	
	A=nx.to_numpy_matrix(G)
	if(A.shape[0]!=0):
		d_g=sum(k[np.array(G.nodes())-1]) #total degree sum. k total degrees
		k_g=G.degree().values()
		print '0000000----------------0000000000'
		print k_g[0]-k[0]*(float(d_g)/float(two_m))
		print '0000000----------------0000000000'
		ind=np.array(np.array(G.nodes())-1)
		B_g=np.array([[ B.item(i,j)-(delta(i,j,groups)*sum(B[i][l] for l in ind)) for j in range(ind.shape[0])] for i in range(ind.shape[0])])
		#B_g=np.array([[ B.item(i,j)-sum(B[i][l] for l in G) for j in range(G.shape[0])] for i in range(G.shape[0])])
					
		#B_g=np.array([[B.item((i,j))-(delta(i,j,groups)*(k_g[i]-((k[i]*float(d_g))/float(two_m)))) for j in range(A.shape[1])]for i in range(A.shape[0])])
		print '-------------------Mod Matrix-------------------'
		print B_g
		print '-------------------Mod Matrix-------------------'
		a=A.shape[0]
		'''start=np.array([r.randint(-1,1) for i in range(a)])
		start=np.reshape(start,(a,1))
		dominant_eigenvector=power_method(B_g,start,1000)
		print dominant_eigenvector
		dominant_eigenvalue=compute_eigenvalue(B_g,dominant_eigenvector)
		dominant_eigenvalue=dominant_eigenvalue.item((0,0))
		if(not(is_positive(dominant_eigenvalue))):
			print "eigenvalue obtained to be negative....recalculating"
			lead_eigenvector=recalculate_lead_eigen(B_g,dominant_eigenvalue)
			lead_eigenvalue=compute_eigenvalue(B_g,lead_eigenvector)
		else:
			lead_eigenvector=dominant_eigenvector
			lead_eigenvalue=dominant_eigenvalue	
		'''
		#lead_eigenvalue,lead_eigenvector=scipy.sparse.linalg.eigs(B_g,k=1)
		eigenvalue,eigenvector=np.linalg.eigh(B_g)
		#print '---------------eigen values---------------'
		#print eigenvalue
		#print '---------------eigen values---------------'
		#print eigenvector
		lead_eigenvalue=max(eigenvalue)
		lead_eigenvector=eigenvector[:,np.argmax(eigenvalue)]
		lead_eigenvalue=np.real_if_close(lead_eigenvalue)
		lead_eigenvector=np.real_if_close(lead_eigenvector)
		print "Lead Eigen Vector"
		print lead_eigenvector
		print "Lead Eigen Value"
		print lead_eigenvalue	
		if(is_positive(lead_eigenvalue)):
			index_vector_list=[]

			for i in range(lead_eigenvector.shape[0]):
				if(is_positive(lead_eigenvector[i])):
					index_vector_list.append(1)
				else:
					index_vector_list.append(0)
		
			s=np.array(index_vector_list)
			
			s=np.reshape(s,(a,1))
			indices=np.array([[i for i, x in enumerate(s) if x == j]for j in (0,1)])
			print "-------------------indices----------"
			print indices
			print "-------------------indices----------"
			temp=np.dot(B_g,s)
			modularity=np.dot(s.T,temp)
			print "----------------modularity----------------"
			print modularity
			print "----------------modularity----------------"
			if(is_positive(modularity)):
				
				nodeList=G.nodes()
				nodeList=np.array(nodeList)
				grp1=nodeList[indices[0]]
				grp2=nodeList[indices[1]]
				#print grp1
				#print grp2
				#print "properly working"
				#print group_number
				G1=G.subgraph(grp1)
				G2=G.subgraph(grp2)
				group_number+=1
				groups[grp1]=group_number
				group_number+=1
				groups[grp2]=group_number
				#print G1.number_of_nodes()
				#print G2.number_of_nodes()
				(groups,group_number)=leadEigenVectorMethod(G1,B,k,two_m, groups,group_number)
				(groups,group_number)=leadEigenVectorMethod(G2,B,k,two_m, groups,group_number)
				print groups
				print "-------------------return 0----------"
				return (groups,group_number)
			else:
				#print groups
				print "-----------------return 1-------------"
				return (groups,group_number)
		else:
			#print groups
			print "return 2"
			return (groups,group_number)
	else:
		print "return 3"
		return (groups,group_number)		
			
			
			
#reading the dataset from the file in gml format
G=nx.read_gml("adjnoun.gml")
A=nx.to_numpy_matrix(G)
print A.shape

#obtaining the dictionary of (node_key : degree) of the graph
degree_dict=G.degree()

#obtaining the list of degrees only
degree_list_values=np.array(degree_dict.values())

#sum of the degrees of all the nodes in the graph		
two_m=sum(degree_list_values)

#jth element of the ith column in degree_mat stores the product of the degrees of ith node and jth node
degree_mat = np.array([[degree_list_values[i]*degree_list_values[j] for j in range(len(degree_list_values))] for i in range(len(degree_list_values))])


#elements of degree_mat_2m stores the corresponding elements of degree_mat each element divided by 2m (sum of the degree of all nodes)
degree_mat_2m=degree_mat/float(two_m)


modularity_mat = np.array([[A.item((i,j))-degree_mat_2m.item((i,j)) for j in range(A.shape[1])] for i in range(A.shape[0])])

groups=np.zeros(G.number_of_nodes()+1) 
groups=leadEigenVectorMethod(G,modularity_mat,degree_list_values,two_m, groups,1)
gnumbers=np.unique(groups[0])
#n=gnumbers.size[0]
final_grouping=np.array([[i for i, x in enumerate(groups[0]) if x == j]for j in gnumbers])
ngroups=groups[1]
print groups
print '------------------------final--------------------------'
print final_grouping
print gnumbers


