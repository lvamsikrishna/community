from igraph import *
import networkx as nx
import numpy as np
import random as r
import math
import pylab





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
#		print x
#		print eigenvalue 
#		print compute_eigenvalue(A,x)
		return x 
		

def compute_eigenvalue(A,x): 
	Axx=x.T*A*x
	xx=x.T*x
	result= Axx/xx
	return result

def delta(i,j,groups):
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

def leadEigenVectorMethod(G,B,degrees,degree_sum,groups,group_number):
	#converting the graph (in networkx format) to numpy matrix format
	#This is equivalent to obtaining the adjacency matrix of the graph
	A=nx.to_numpy_matrix(G)
	print "A.shape"
	print A.shape
	if(A.shape[0]!=0):
		print "G nodes"
		print G.nodes()
		total_degrees=degrees[np.array(G.nodes())-1]
		d_g=sum(total_degrees)
	
		#obtaining the dictionary of (node_key : degree) of the graph
		degree_dict=G.degree()

		#obtaining the list of degrees only
		degree_list_values=degree_dict.values()
		print "degree_list_values"
		print degree_list_values
		#sum of the degrees of all the nodes in the graph
		print "two_m"		
		two_m=sum(degree_list_values)
		print two_m
		#jth element of the ith column in degree_mat stores the product of the degrees of ith node and jth node
		degree_mat = np.array([[degree_list_values[i]*degree_list_values[j] for j in range(len(degree_list_values))] for i in range(len(degree_list_values))])


		#elements of degree_mat_2m stores the corresponding elements of degree_mat each element divided by 2m (sum of the degree of all nodes)
		degree_mat_2m=degree_mat/float(two_m)

	
		#modularity_mat = np.array([[A.item((i,j))-degree_mat_2m.item((i,j)) for j in range(A.shape[1])] for i in range(A.shape[0])])
		modularity_mat=np.array([[B.item((i,j))-(delta(i,j,groups)*(degree_list_values[i]-(degrees[i]*(d_g/degree_sum)))) for j in range(A.shape[1])]for i in range(A.shape[0])])
	
		#print "Modularity Matrix"
		#print modularity_mat

		'''u,v=np.linalg.eigh(modularity_mat)
		print "Eigen Values"
		print u'''


		#A.shape[0] gives the number of rows in the matrix A, storing it in a
		a=A.shape[0]

		#generating a random list whose size is equal to the number of rows, where each element ranges from -1 t0 1
		myList=[r.randint(-1,1) for i in range(a)]

		#generating a row vector from the elements in the list, number of rows is equal to the number of rows in the network. This row vector serves as the initial vector for power iteration method
		start=np.matrix(myList)

		#start now stores the transpose of the previous start
		start=start.T

		#calling the power_method which returns the eigenvector corresponding to the eigenvalue of largest magnitude
		dominant_eigenvector=power_method(modularity_mat,start,1000)
		#dominant_eigenvector stores the eigenvector of the modularity matrix corresponding to the eigenvalue of highest magnitude

		#compute_eigenvalue takes as arguments, a square matrix and an eigenvector of the matrix and returns eigenvalue corresponding to that particular eigenvector of the matrix as a 1X1 matrix
		dominant_eigenvalue=compute_eigenvalue(modularity_mat,dominant_eigenvector)

		#dominant_eigvalue store the eigenvalue corresponding to the dominant_eigenvector


		#obtaining the eigenvalue as a number instead of a matrix,
		#matrix.item(i,j) returns the jth element in the ith row of the matrix
		dominant_eigenvalue=dominant_eigenvalue.item((0,0))


		if(not(is_positive(dominant_eigenvalue))):
			print "negative....recalculating the lead eigenvector using power method"
			lead_eigenvector=recalculate_lead_eigen(modularity_mat,dominant_eigenvalue)
			lead_eigenvalue=compute_eigenvalue(modularity_mat,lead_eigenvector)
		else:
			lead_eigenvector=dominant_eigenvector
			lead_eigenvalue=dominant_eigenvalue	
		print "Lead eigen vector "
		print lead_eigenvector
		print "Lead Eigen Value"
		print lead_eigenvalue
		#print lead_eigenvector
		if(is_positive(lead_eigenvalue)):
			index_vector_list=[]

			for i in range(lead_eigenvector.shape[0]):
				for j in range(lead_eigenvector.shape[1]):
					if(is_positive(lead_eigenvector.item(i,j))):
						index_vector_list.append(1)
					else:
						index_vector_list.append(-1)
		
			#print index_vector_list
			index_vector=np.matrix(index_vector_list)
			index_vector=index_vector.T
			indices=np.array([[i for i, x in enumerate(index_vector) if x == j]for j in (-1,1)])		
			#print indices
		
			'''for i in range(2):
				print "Cluster "+str(i)
				for j in range(len(indices[i])):
					print indices[i][j]+1'''
			modularity=index_vector.T*modularity_mat*index_vector
			print "Modularity"
			print modularity
			if(is_positive(modularity)):
			
				'''grp1=indices[0]
				grp1=np.array(grp1)+1
				grp1.shape
				grp2=indices[1]
				grp2=np.array(grp2)+1
				grp2.shape'''
				nodeList=G.nodes()
				nodeList=np.array(nodeList)
				grp1=nodeList[indices[0]]
				grp2=nodeList[indices[1]]
				print grp1
				print grp2
				group_number+=1
				groups[grp1]=group_number
				group_number+=1
				groups[grp2]=group_number
				print "properly working"
				print group_number
				G1=G.subgraph(grp1)
				G2=G.subgraph(grp2)
				print G1.number_of_nodes()
				print G2.number_of_nodes()
				(groups,group_number)=leadEigenVectorMethod(G1,B,degrees,degree_sum, groups,group_number)
				(groups,group_number)=leadEigenVectorMethod(G2,B,degrees,degree_sum, groups,group_number)
				print groups
				return (groups,group_number)
			else:
				print groups
				return (groups,group_number)
		else:
			print groups
			return (groups,group_number)
	else:
		print groups
		return (groups,group_number)




#reading the dataset from the file in gml format
G=nx.read_gml("football.gml")
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
groups=leadEigenVectorMethod(G,modularity_mat,degree_list_values,two_m, groups,0)
print "----------------final-----------------"
print groups
print groups[0].shape




'''
	x=np.arange(1,A.shape[0]+1)
	x=np.reshape(x,A.shape[0])
	y=lead_eigenvector.T
	y=np.array(lead_eigenvector)
	y=y.reshape(A.shape[0])
	print x.shape
	print y.shape
	pylab.plot(x,y,'ro')
	pylab.xlim(0,A.shape[0]+1)
	#pylab.ylim(-1,1)
	pylab.xlabel('Node Numbers')
	pylab.show('Eigen Vectors')
'''




#Community structure detection method inspired from Kernighan Lin's method for graph partioning
'''In this method, an initial random division of network is chosen. Though an obvious choice would be to make the graph into two groups such that, all the elements lies in group and none lies in the other, here we start with the division of network that has been obtained as a result of the above algorithm. Then the following steps are performed as a part of the algorithm. A node is selected at random and '''



"""L=np.matrix([[1,2,0],[0,2,0],[0,0,3]])
print L
start=np.matrix([1,1,1]).T
print start
"""

"""for i in xrange(result.shape[0]): 
	result[i]=round(result[i])
print result"""
#print start.shape
#res=power_method(L,)
#astro=Graph.Load("karate.net")
#summary(astro)
#L=Graph.get_adjacency(astro,True)
#Graph.write_adjacency(astro,"op1",' ',';',True)
#print L
#w,v=LA.eigh(L)
#print("Largest eigenvalue:", max(w))
#print("Smallest eigenvalue:", min(w))


