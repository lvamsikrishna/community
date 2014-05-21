from igraph import *

karate=Graph.Read_GML("astro-ph.gml")
cg1=karate.community_leading_eigenvector()

print 

for component in cg1:
	print component
print "vamsi krishna"
s="aeigop"
i=0
for component in cg1:
	#print "component"
	i=i+1
	for g in component:
		f=open(s+str(i),'w')
		f.write(str(g)+'\n')
print i
	
