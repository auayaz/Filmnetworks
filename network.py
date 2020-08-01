import networkx as nx
import os
import scipy.io
import numpy as np

"""
This class takes the bridgenetwork found with the matlab routine capillarypressure_bridges and 
utilizes the class networkx to calculate measures related to graph theory.
"""

class network:

	def __init__(self,experiment):

		self.experiment = experiment
		self.path = '/home/monem/Dropbox/code/MATLAB/osloexperiments/' + experiment
		self.filename = 'bridgecoordinates.mat'
		self.cwd = os.getcwd()
		self.latexfile = '/home/monem/Dropbox/code/latex/article_filmflow/'
		self.figure_graph = '/home/monem/Dropbox/code/python/filmnetworks/figure_graphs/'
		self.data_network_size = '/home/monem/Dropbox/code/python/filmnetworks/data_network_size/'
		self.data_angledistribution = '/home/monem/Dropbox/code/python/filmnetworks/angle_distribution/'

	def create_graph(self):

		"""
		This function takes in Vertex ID and its respective spatial coordinates, and spits
		out a graph object. G is the Graph object connecting all the links, and H gives connectivity
		to the trapped clusters, F is the superimposed Graph.
		"""
		
		os.chdir(self.path)
		mat = scipy.io.loadmat(self.filename)
		coordinates = mat['E_coordinates']						#Coordinates of the vertex ID
		E =mat['E']												#Vertex ID
		bridge_length = mat['bridge_l'][0]
		cluster_E_coordinate = mat['cluster_E_coordinates']
		cluster_E = mat['cluster_E']
		snap_E = mat['snap_E']
		snap_E_coordinates = mat['snap_E_coordinates']

		G = nx.Graph(); i = 0;

		#flip the graph with respect to y.
		maks_y =  np.max(coordinates)
		coordinates[:,1] = np.abs(coordinates[:,1]-(maks_y-1))
		coordinates[:,3] = np.abs(coordinates[:,3]-(maks_y-1))
		cluster_E_coordinate[:,1] = np.abs(cluster_E_coordinate[:,1] - (maks_y-1))
		cluster_E_coordinate[:,3] = np.abs(cluster_E_coordinate[:,3] - (maks_y-1))
		snap_E_coordinates[:,1] = np.abs(snap_E_coordinates[:,1] - (maks_y-1))
		snap_E_coordinates[:,3] = np.abs(snap_E_coordinates[:,3] - (maks_y-1))


		while i<len(E):
			
			tmp = coordinates[i]
			vertex_N = E[i]
			l = np.sqrt((tmp[0] - tmp[2])**2 + (tmp[1] - tmp[3])**2)
			G.add_node(vertex_N[0],pos=(tmp[0],tmp[1]))
			G.add_node(vertex_N[1],pos=(tmp[2],tmp[3]))
			G.add_edge(vertex_N[0],vertex_N[1], weight=1, length =l, resistance=1)
			i +=1

		#Here we add the connectivity with the trapped clusters.

		H = nx.Graph()
		k = 0;
		while k<len(cluster_E):

			tmp = cluster_E_coordinate[k]
			vertex_N = cluster_E[k]
			#l = np.sqrt((tmp[0] - tmp[2])**2 + (tmp[1] - tmp[3])**2)
			H.add_node(vertex_N[0],pos=(tmp[0],tmp[1]))
			H.add_node(vertex_N[1],pos=(tmp[2],tmp[3]))
			H.add_edge(vertex_N[0],vertex_N[1], weight=1, length = l, resistance=1)
			k +=1
		
		#Here we add the connectivity of the snap-off with the trapped clusters
		S = nx.Graph()
		s = 0;
		while s<len(snap_E):

			tmp = snap_E_coordinates[s]
			vertex_N = snap_E[s]
			S.add_node(vertex_N[0],pos=(tmp[0],tmp[1]))
			S.add_node(vertex_N[1],pos=(tmp[2],tmp[3]))
			S.add_edge(vertex_N[0],vertex_N[1], weight=1, length = l, resistance=1)
			s+=1

		#combining the two different graphs.
		F = nx.compose(G,H) #connecting capillary bridges to clusters
		J = nx.compose(F,S) #connecting capillary bridges to clusters and snap-off
		os.chdir(self.cwd)    	

		print "graph F has %d nodes with %d edges" % (nx.number_of_nodes(F), nx.number_of_edges(F))
		print "graph G has %d nodes with %d edges" % (nx.number_of_nodes(G), nx.number_of_edges(G))
		print "graph J has %d nodes with %d edges" % (nx.number_of_nodes(J), nx.number_of_edges(J))
		print "graph S has %d nodes with %d edges" % (nx.number_of_nodes(S), nx.number_of_edges(S))

		print "G has %d connected components"%(nx.number_connected_components(G))
		print "F has %d connected components"%(nx.number_connected_components(F))
		print "J has %d connected components"%(nx.number_connected_components(J))
		print "S has %d connected components"%(nx.number_connected_components(S))

		return G,F

	def network_node_size_distribution(self):

		G,F,J = self.create_graph()

		"""
		This function finds the number of nodes in subsets of the graphs of G and F
		"""

		from networkx.algorithms.components.connected import connected_components

		s_G =  [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]
		l_G = sorted(nx.connected_components(G), key=len, reverse=True)
		s_F =  [len(c) for c in sorted(nx.connected_components(F), key=len, reverse=True)]
		l_F = sorted(nx.connected_components(F), key=len, reverse=True)

		tmp = nx.connected_components(G)
		import collections
		counter = collections.Counter(s_G)
		v_G = counter.values()
		k_G = counter.keys()
		counter = collections.Counter(s_F)
		v_F = counter.values()
		k_F = counter.keys()
		
		#Plotting the edge size distribution.
		import matplotlib.pyplot as plt

		plt.figure(1)
		plt.loglog(k_F,v_F,'o-')
		plt.ylabel('log[N(s)]')
		plt.xlabel('log[s]')
		plt.savefig(self.experiment+"Distributionlog",bbox_inches= 'tight',format = 'eps', dpi=300)

		plt.figure(2)
		plt.loglog(k_G,v_G,'o-')
		plt.ylabel('N(s)')
		plt.xlabel('s')
		plt.savefig(self.experiment+ "Distribution",bbox_inches= 'tight',format = 'eps',dpi= 300)
		plt.show()

		#Write to file
		np.savetxt(self.data_network_size + self.experiment + "_node_distribution_G"+".txt",np.c_[k_G,v_G])
		np.savetxt(self.data_network_size + self.experiment + "_node_distribution_F"+".txt",np.c_[k_F,v_F])

	def network_edge_size_distribution(self):

		"""
		This function finds the edge_size distribution.
		"""		
		from collections import Counter
		G,F,J = self.create_graph()
		
		# F #
		for graph in ['F','G','J']:

			tmp = eval(graph)
			graphs = list(nx.connected_component_subgraphs(tmp))
			edge_array = np.zeros(len(graphs))

			for i in range(len(graphs)):
				tmp =  list(graphs[i])
				sub = F.subgraph(tmp)
				edge_array[i] = sub.number_of_edges()
		
			node_edges = Counter(edge_array)
			F_k =  node_edges.keys()
			F_v =  node_edges.values()

			np.savetxt(self.data_network_size + self.experiment + "_edge_distribution_"+graph+".txt",np.c_[F_k,F_v])                                                

		
	def Bondingbox(self):
		"""
		Redundant function.
		"""
		os.chdir(self.path)
		mat = scipy.io.loadmat(self.filename)
		xbox = mat['xbox'][0]
		zbox = mat['zbox'][0]

		import matplotlib.pyplot as plt
		import pylab as pl
	
		hist, bins = np.histogram(xbox, bins='auto')
		width = 0.7 * (bins[1] - bins[0])
		center = (bins[:-1] + bins[1:]) / 2

		plt.figure(2)
		#plt.bar(center, hist, align='center', width=width)
		print len(center)
		print len(hist)
		plt.loglog(center,hist,'.')
		plt.show()

		np.savetxt(self.experiment + "_boundingbox"+".txt",np.c_[xbox,zbox])

		#plt.savefig(self.experiment+"boundingbox",bbox_inches= 'tight',format = 'eps', dpi=300)

	def isomorphism(self):

		#Under progress.
		G,F = self.create_graph()		
		s = [sorted(nx.connected_components(G), key=len, reverse=True)]
		print s

		print "isomorphism"

	def node_degree_distribution(self):

		"""
		This function gives us the node-degree distribution of the Graphs F and G
		"""

		G,F,J = self.create_graph()
		from collections import Counter

		for graph_name in ['G','F','J']:
			tmp = eval(graph_name)
			nd = list(tmp.degree())
			nd_array = np.zeros(len(nd))				

			for k in range(len(nd)):
				nd_array[k] = nd[k][1]

			node_degree = Counter(nd_array)
			n_d_k =  node_degree.keys()
			n_d_v =  node_degree.values()
					
			np.savetxt(self.data_network_size + self.experiment + "_node_degree_distribution"+graph_name+".txt",np.c_[n_d_k, n_d_v])
		
		#Finding the mean value of k (node degree)

			print n_d_v
			print n_d_k

		n_d_v = np.array(n_d_v)
		n_d_k = np.array(n_d_k)
		value = n_d_v*n_d_k
		v= sum(value)
		tmp = sum(n_d_v)
		print v/tmp
		"extracing the degree from the vertices"
		#Saving to textfiles
		#np.savetxt(self.data_network_size + self.experiment + "_node_degree_distribution"+".txt",np.c_[n_d_k, n_d_v])

	def Average_path_length(self):

		G,F = self.create_graph()
		tmp = nx.average_shortest_path_length(F)
		print tmp


	def spectralanalysis(self):
		
		"""
		Use in progress
		In graph theory and computer science, an adjacency matrix is a square matrix used to
		represent a finite graph. The elements of the matrix indicate if pairs of vertices
		are adjacent or not in the graph.
		"""

		G,F = self.create_graph()
		#a = nx.adjacency_matrix(G)
		eigenvalues = nx.adjacency_spectrum(G)
		for i in eigenvalues:
			print i

		print len(eigenvalues)

	def clustering(self):

		import random
		G,F = self.create_graph()

		#A = nx.average_clustering(G,1000)
		#G2 = nx.G.to_undirected()
		#G = nx.G.to_undirected()
		C_local = nx.clustering(F)					#Gives the local clusering coefficient 

		print "clustering ="
		print C_local
		print nx.average_clustering(F)

	def drawing_the_graph(self):

		import random
		import collections
		from matplotlib import pyplot as plt

		#gui_env = ['TKAgg','GTKAgg','Qt4Agg','WXAgg']
		#matplotlib.use(gui_env[0],warn=False, force=True)
		#from matplotlib import pyplot as plt
		#print "Switched to:",matplotlib.get_backend()

		G,F = self.create_graph()
		count = 1;
		for graph_name in ['F','G']:
			C=nx.connected_component_subgraphs(eval(graph_name))
			s = nx.number_connected_components(eval(graph_name))

			fig = plt.figure(count,figsize=(2.5*3,3.5*3))

			for i in C:        
				c=[random.random()]*nx.number_of_nodes(i) # random color...
				pos = nx.get_node_attributes(i,'pos')
				nx.draw_networkx_edges(i,pos, edge_color='k')
				nodes=nx.draw_networkx_nodes(i,pos, node_color = c, cmap = 'jet',vmin= 0, vmax =1,  node_size = 10,linewidths = 0.7, alpha =1)
				nodes.set_edgecolor('k')
			plt.axis('off')		
			plt.savefig(self.figure_graph+self.experiment+ 'Graph'+'_'+graph_name+'.eps',bbox_inches='tight',dpi=fig.dpi, pad_inches=0.0,
				frameon=None ,format = 'eps')

			count += 1

	def angle_distribution(self):

		"""
		This function finds the angle distribution of the capillary bridges.
		"""

		os.chdir(self.path)
		mat = scipy.io.loadmat(self.filename)
		coordinates = mat['E_coordinates']	#Coordinates of the vertex ID [y0,x0,y1,x1]
		E =mat['E'] 						#Vertex ID
		angles_alongx = np.zeros(len(coordinates))
		angles_alongy = np.zeros(len(coordinates))

		import math
		
		for i in range(len(coordinates)):

			x0 = coordinates[i][1]; y0 = coordinates[i][0]
			x1 = coordinates[i][3]; y1 = coordinates[i][2]

			if y1>y0:
				uy = y1-y0
			else:	 
				uy = y0-y1
			
			if x1>x0:
				ux = x1-x0
			else:
				ux = x0-x1
			
			#uy  = y1 - y0
			#ux = x1-x0
			unorm = np.sqrt(ux**2 + uy**2)
			vnorm = 1
			tmpy = uy/(vnorm*unorm)
			tmpx = ux/(vnorm*unorm)
			
			angles_alongx[i] = (math.acos(tmpx)*180)/math.pi
			angles_alongy[i] = (math.acos(tmpy)*180)/math.pi


		import matplotlib.pyplot as plt
		plt.figure()
		plt.hist(angles_alongx.ravel(), bins=30, fc='k', ec='k')
		plt.show()
		np.savetxt(self.data_angledistribution + self.experiment + "_angledistribution"+".txt",np.c_[angles_alongx,angles_alongy])


	def test(self):

		import matplotlib.pyplot as plt
		"Attemping to find ways to iterate over the graph object"


	def test2(self):
		G,F = self.create_graph()
		#s = [sorted(nx.connected_components(F), key=len, reverse=True)]
		value = nx.global_efficiency(F)
		print value

def main():

	experiment = ['experiment2_15angle','experiment3_30angle','experiment2_45angle','experiment1_60angle']
	N = network(experiment[3])
	#N.Average_path_length()
	#N.create_graph()
	#N.network_edge_size_distribution()
	#N.network_node_size_distribution()
	#N.angle_distribution()
	#N.Bondingbox()
	N.drawing_the_graph()
	#N.clustering()
	#N.spectralanalysis()
	#N.test()
	#N.node_degree_distribution()
	#N.Assortivity()
	#N.test2()
	#for i in range(len(experiment)):
	#	N = network(experiment[i])
	#	N.network_size_distribution()

if __name__ == '__main__':
	main()