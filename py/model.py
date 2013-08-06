'''
MODEL 
by Rick Page
Stores the parameters, domain information and boundary conditions that define a one-dimensional advection diffusion transient storage model.
Used by Approximation to represent the one-dimensional PDE. 
'''
import numpy as np
from subprocess import call
from pylab import plot,show,savefig,clf,xlabel,ylabel
	
class Model:
	def __init__(self):
		# upstream boundary condition vector C(x,t=0)
		self._upstream = None
		self.upstreamFilename = None

	def __str__(self):
		'''
		Prints domain discretization and model parameters. Also indicates if upstream conditions have been properly set.
		'''		
		return self.info()

	def info(self):
		'''
		Returns a string containing useful information on the state of the model, including discretization information and model parameters.
		'''
		string = str()
		
		# print model parameters				
		string += "Model parameters:\t"
		try:
			string += "U %e\tD %e\tk1 %e\tk2 %e\n" % (self.u, self.d, self.k1, self.k2)
		except AttributeError:
			string += "Not set.\n"		

		# print domain info, and depending if self._upstream is None, print upstream info
		# if not set yet, exceptions will be thrown; catch them and report Not Set
		string += "Domain information:\t"
		try:	
			# temporary flag to store wether an exception occured
			noselfm = False	
			string += "M %d\tN %d\tdt %e\tdx %e\n" % (self.m, self.n, self.dt, self.dx)
		except AttributeError:
			string += "Not set.\n"
			noselfm = True
			if(self._upstream != None):			
				string += "Upstream values:\t Set, but note that domain information is not set." 
			else:
				string 	+= "Upstream values:\t Not set."	
		
		# print status on upstream data only if we know that self.m is set and upstream is not None		
		# at this point, we know self.m is set but we don't know if _upstream is set
		if(noselfm == False and self._upstream != None):	
			if(len(self._upstream) == self.m and self.upstreamFilename == None): 	
				string += "Upstream values:\t Set correctly.\n"
			elif(len(self._upstream) == self.m and self.upstreamFilename != None): 	
				string += "Upstream values:\t Set correctly from %s\n"  % self.upstreamFilename
			elif(len(self._upstream) != self.m): 	
				string += "Upstream values: WARNING number of time-steps M does not equal number of size of upstream vector \n"
		
		return string

	def domain(self,d=None):
		'''
		Gets or sets the values M,N,dt,dx as a length 4 list of floats.
		With no arguments, params returns a length 4 list of floating point values: [M,N,dt,dx]. 
		With params as the argument, where d is a length 4 list of floating point values the model variables 
		are set as M = int(d[0]),N = int(d[1]),dt = float(d[2]),dx = float(d[3])
		'''
		if(d == None):
			try:
				return [self.m,self.n,self.dt,self.dx]
			except AttributeError:
				raise StandardError, "domain() - Domain information has not been set. Please use domain(x) where x is a 4 element list to set M,N,dt,dx."
		elif(type(d) == list and len(d) == 4):
			self.m = int(d[0]);
			self.n = int(d[1]);
			self.dt = float(d[2]);
			self.dx = float(d[3]);
		else:
			raise TypeError, "domain() - Domain information must be input as a 4 element list for U,D,k1,k2"
	def params(self,p=None):
		'''
		Gets or sets the values U,D,k1,k2 as a length 4 list of floats.
		With no arguments, params returns a length 4 list of floating point values: [U,D,k1,k2]. 
		With params as the argument, where p is a length 4 list of floating point values the model variables 
		are set as U = p[0], D=p[1], k1=p[2],k2=p[3].
		'''
		if(p == None):
			try:
				return [self.u,self.d,self.k1,self.k2]
			except AttributeError:
				raise StandardError, "params - Model parameters have not yet been set. Please use params(x) where x is a 4 element list to set U,D,k1,k2."
		elif(type(p) == list and len(p) == 4):
			self.u = float(p[0]);
			self.d = float(p[1]);
			self.k1 = float(p[2]);
			self.k2 = float(p[3]);
		else:
			raise TypeError, "Error: setParams - parameters must be a 4 element list for U,D,k1,k2"
		return

	def upstream(self,x=None):
		'''
		Gets or sets the length M numpy array of upstream values depending on the input.
		If the input to upstream is None, it returns the upstream vector as a numpy array size 1 by M.
		If the input is a list or numpy ndarray, the values are copied to the upstream vector. 
		If the input argument is a string, it is treated as a file name and the method readUpstreamFile
		is used to load the upstream data.
		'''
		if(type(x) == np.ndarray):
			self._upstream = x.copy()
		elif(type(x) == list):
			self._upstream = np.array(x)
		elif(type(x) == str):
			self.readUpstreamFile(x)
		elif(x == None):
			return self._upstream
		else:
			raise TypeError, "Error: upstream(x) - accepts either a list or a numpy.ndarray, or a filename containing concentration values.\nTo get instead of set the values, simply use upstream() with no arguments."
		
		return

	def readUpstreamFile(self,filename):
		'''
		Read vector of upstream boudary values from file as a numpy array of floats and store internally. Also sets the
		internal string upstreamFilename to the filename .
		'''
		# try to open the file and read the contents; on failure print helpful error message
		try:		
			F = open(filename,'r')
			contents = F.read().split()
			F.close()
		except IOError:
			print "Error: readUpstreamFile - could not read", filename
			return
		vector = np.array([float(contents[i]) for i in range(len(contents))])
		# inform user of success, and then set the filename and vector values
		print '\nSuccessfully read ',len(vector), ' elements from ', filename, '\n' 
		self.upstreamFilename = filename
		self.upstream(vector)
		return

	def copy(self,aModel):
		'''
		Sets attributes equal to attributes of the Model instance aModel. If input is not of type Model, an error is generated.
		'''		
		try:
			# copy parameters
			self.params(aModel.params())
			# copy domain info
			self.domain(aModel.domain())
			# copy upstream info
			# since upstream is an numpy array use np.copy
			self.upstream(np.copy(aModel.upstream()))
			# also get filename
			self.upstreamFilename = aModel.upstreamFilename
		except AttributeError:
			raise TypeError, "copy() - argument aModel must be an instance of class Model"
		return

