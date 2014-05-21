from igraph import *
karate=Graph.Read_Pajek("karate.net")
cg1=karate.community_leading_eigenvector()

for component in cg1:
	print component

print(cg1.modularity)	
print "Dividing into two communities"

cg2=karate.community_leading_eigenvector(2)
for comp in cg2:
	print comp
print(cg2.modularity)	
	
	
