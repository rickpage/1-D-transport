'''	
Approximation
by Rick Page
	The base class for the OTIS and DISCUS solvers. The solver, an external executable, 
generates an approximation to the downstream values and writes them to a file rb.dat. To run 
the solver, lb.dat must conatin upstream values. Downstream values can be loaded with the 
routine upstream(). Approximations are generated with the routine run(), which in the base 
class simply reports an error message.
	In order to calculate the norm, the observed concentration must be loaded. To load
a two-column file, where the left column is upstream data and the (observed) downstream data 
in the right column, use loadTwoColumn. Files specially formatted for the module have M,N,dx
and dt as the first four lines, upstream data in the left column and (observed) downstream
data in the right column, and can be loaded with loadFile(filename).
'''
from model import *

class Approximation:
	def __init__(self):
		# make instance of the river or stream being modelled
		self.stream = Model()
		# observed and approximate values 
		self._observed = np.array([])
		self._approx = np.array([])
		self.norm = -1
		self.relnorm = -1
		# define a flag that is set to true when requisite lb.dat has been written.
		# this is to prevent unintentionally running a simulation with a previous set of data,
		# otherwise simply checking for the file would suffice
		self.lbWritten = False
		#  assign a flag to see if the simulation has been run
		self.simRun = False
		# This variable is set to the path of the executable file
		self.binPath = '../bin'
		self.binName = None

	def __str__(self):
		return self.info()

	def info(self):
		'''
		Returns a string containing model and solver information, including the status of the
		simulation.
		'''
		# add path, name, stream.info, observed status and simulation status to a string
		string = str()
		string += "Path:\t\t%s\n" % self.binPath
		string += "Execcutable:\t%s\n" % self.binName
		string += self.stream.info()
		string += '\n'
		if(len(self._observed) > 0):
			if(len(self._observed) == self.stream.m):
				string += "Observed:\tSet correctly.\n"
			else:
				string += "Observed: WARNING number of time-steps M does not equal number of size of observed vector"
		else:
			string += "Observed:\tNot set (optional).\n"		
		if(self.ready()):
			string += "Status:\tSimulation may be run.\n"
		elif(self.ready() == False and self.lbWritten==False):
			string += "Status:\tSimulation is not ready because lb.dat has not been written (use Approximation.upstream() to write it).\n"
		if(self.simRun):
			string += "Note: Results have been generated\n"	
		return string
	
	def approximated(self):
		'''
		Returns approximation to downstream vector if calculated. Returns None if approximation has not been run.
		'''
		if(len(self._approx) == 0):
			return None
		return self._approx

	def observed(self,x=None):
		'''
		Gets or sets the length M numpy array of observed values. 
		If the input to observed is None, it returns the observed vector as a numpy array size 1 by M.
		If the input is a list or numpy array, the values are copied to the observed vector. 
		If the input argument is a string, it is treated as a file name containing the data.
		'''
		if(type(x) == np.ndarray):
			self._observed = x
		elif(type(x) == list):
			self._observed = np.array(x)
		elif(x == None):
			return self._observed
		elif(type(x) == str):
			# attempt to open the file and read its contents into a list; if the file cannot be found, display a custom error message
			try:		
				F = open(filename,'r')
				contents = F.read().split()
				F.close()
			except IOError:
				raise IOError, "observed() - could not read ", filename
			# assign the elements
			vector = np.array([float(contents[i]) for i in range(len(contents))])
			# confirm to the user that the elements have been read			
			print '\nSuccessfully read ',len(vector), ' elements from ', filename, '\n' 
			# set the vector
			self.observed(vector) # Recursion is nice
		else:
			raise TypeError, "observed() - accepts either a list or a numpy.ndarray, or a filename containing concentration values.\n\
			To get the values instead of set the values, simply use observed() with no arguments."
		return

	def copyStream(self,stream):
		try:
			if(stream == None):
				raise ValueError, "copyStream() - argument must be instance of Model class\n\
					 (try spam.stream where spam is your Approximation instance)"
			else:
				self.stream.copy(stream)
				us = self.stream.upstream()				
				# the C binaries depends on the file lb.dat; write upstream values to this		
				outfile = open('lb.dat','w');
				[outfile.write('%e ' % us[i]) for i in range( int(self.stream.m) )]
				outfile.close()
				self.lbWritten = True
		except AttributeError:
			raise ValueError, "copyStream() - argument must be instance of Model class\n\
					 (try spam.stream where spam is your Approximation instance)"
		return

	def upstream(self,x=None):
		'''
		Gets or sets the length M numpy array of upstream values depending on the input.
		If the input to upstream is None, it returns the upstream vector as a numpy array size 1 by M.
		If the input is a list or numpy ndarray, the values are copied to the upstream vector. 
		If the input argument is a string, it is treated as a file name and the method readUpstreamFile
		is used to load the upstream data.
		'''
		# call Model method to set upstream
		if(x == None):		
			return self.stream.upstream()
		else:	
		# attempt to set upstream
			self.stream.upstream(x)	
			# if successful write it to a file	
			us = self.stream.upstream()				
			# the C binaries depends on the file lb.dat; write upstream values to this		
			outfile = open('lb.dat','w');
			[outfile.write('%e ' % us[i]) for i in range( int(self.stream.m) )]
			outfile.close()
			self.lbWritten = True
		return
		

	def setBinName(self,string):
		'''
		Sets the directory that the external executables (solvers) reside in. Can be either an 
		absolute path or a path relative to the present working directory.
		'''
		# make sure input is string, then set it
		if(type(string) != str):
			raise ValueError,"Error: setBinName() - input must be a string"
		self.binName = string
		return

	def setBinDir(self,string):
		'''
		Sets the directory that the external executables (solvers) reside in. Can be either an 
		absolute path or a path relative to the present working directory. By default, it is set 
		to ../bin/. 
		'''
		# make sure input is string, then set it
		if(type(string) != str):
			raise ValueError,"setBinDir() - input must be a string"
		self.binPath = string 
		return

	def ready(self):
		'''
		Returns True only if the lb.dat file has been written, the binary name and path have been set,
		domain information has been set, and model parameters have been set. Otherwise returns False.
		The solver will not run if ready returns False.
		'''
		# an easy way to check params and domain programatically is to call them and catch the Error
		try:
			self.stream.params()
			self.stream.domain()
		except:
			return False	
		return self.lbWritten and self.binName != None and self.binPath != None

	
	def getNorm(self):
		'''
		Returns the absolute error norm and relative error norm as a tuple. Sets attributes Approximation.norm and Approximation.relnorm. 
		'''
		# make sure observed and approx are both set
		# if observed is still an empty vector, length of vector will be zero, so raise exception
		# also raise exception when lengths of vectors differ or simulation has not been run
		if(len(self.observed())==0):
			raise StandardError, "Error: Norm cannot be calculated because the observed data has not been loaded. Use observed(filename),\
				 loadTwoColumnFile(filename), or loadJRMFile(filename) to load the observed data."
		elif(len(self._approx) != len(self._observed)):
			raise ValueError, "Error: Cannot calculate norm; approximation and observed vectors differ in length."
		elif (self._approx.size == 0 or self.simRun == False):
			raise ValueError, "Error: Cannot calulcate norm; approximated vector not calculcated." 
		elif (self._observed.size == 0):	
			raise ValueError, "Error: Cannot calulcate norm; observations not set."
		self.norm = np.linalg.norm(self.approximated()-self._observed,2)
		self.relnorm = np.linalg.norm((self.approximated()-self._observed),2)/np.linalg.norm(self._observed,2)
		return (self.norm,self.relnorm)

	def loadTwoColumnFile(self,filename):
		'''
		Used for file format:
		upstream1	downstream1
		...		...
		upstreamM	downstreamM
		
		Parse the file given by filename, stores the upstream and
		the observed downstream values. Also writes lb.dat and sets the lbWritten flag 
		to True. This file is needed in order to run the Otis or Discus solvers.
		'''
		# open the two column file and put its lines in a list of strings
		infile = open(filename,'r')
		lines = infile.read().splitlines()
		infile.close()
		# the number of lines is the number of observations so set the Model domain information
		self.stream.m = int(len(lines))
		# add upstream values to the Model instance and known downstream values to the observed vector		
		self.stream.upstream(np.array([float(lines[j].split()[0]) for j in range(len(lines))]))
		self.stream.upstreamFilename = filename
		self.observed(np.array( [float(lines[j].split()[1]) for j in range(len(lines))] ))
		
		return
	
	def __call__(self):
		'''
		Returns the result of run().
		'''
		return self.run()
 
	def run(self):
		'''
		This function calls the external solver if ready() returns True. Returns False if the solver did not
		succeed, and returns True if results were recorded. Uses the solver output file rb.dat to load data.
		Note that if ready() is True, the file rb.dat will be removed (if it exists) when run() is called, 
		and a new rb.dat will be written only if the solver does not terminate computations due to instability
		or incorrect arguments. info() or "print spam()" where spam is the name of the class instance
		can be helpful when run() fails.
		'''
		# Exit immediately if the solver isn't fully initialized		
		if not self.ready():
			return False
		# prepare to call the solver by removing rb.dat; this way, if the solver fails, we know because there is no results file
		try:
			call('rm rb.dat',shell=True);
		except:
			pass		
		# make sure binary path ends with a forward slash
		if (self.binPath[-1] != '/'):
			self.binPath += '/'
		# call solver
		call([self.binPath + self.binName,str(int(self.stream.m)),str(int(self.stream.n)),\
			str(float(self.stream.u)),str(float(self.stream.d)),str(float(self.stream.k1)),str(float(self.stream.k2)),\
			str(float(self.stream.dt)),str(float(self.stream.dx))])
		# if rb.dat was not found then the solver surely failed
		# try opening for read and on success return True; else return false on exception		
		try:
			infile = open("rb.dat","r"); contents = infile.read().split(); infile.close()
			self._approx = np.array([float(contents[i]) for i in range(len(contents))]);
			self.simRun = True
			return True
		# only handle IOErrors in case
		except IOError:
			print "Warning: Results not available - solver failed to complete with given parameters."
			return False

	def loadFile(self,string):
		'''
		Loads M,N,dt, and dx, as well as upstream boundary values and downstream observations.
		Used for the following file format:
		M
		N
		dt
		dx
		upstream1 	downstream1
		...		...
		upstreamM 	downstreamM

		Parses a file designed for the module. The domain info is loaded from the first
		four lines, and a two column format follows. loadFile() stores the upstream and
		the observed downstream values. Also writes lb.dat and sets the lbWritten flag 
		to True. This file is needed in order to run the Otis or Discus solvers.
		'''
		# open file
		F = open(string,'r')
		# get rows of file and put into list
		lines = F.read().splitlines()
		F.close()
		# get dt,M,L,N from first four lines
		# use list comprehension to make compact loop and grab the first number
		values = [float(lines[i].split()[0]) for i in range(4)]
		# assign values
		self.stream.domain(values)		
		# get the concentrations by using split on the remaining lines				
		self.stream.upstream( np.array([float(lines[4+i].split()[0]) for i in range(self.stream.m)]) )
		self.stream.upstreamFilename = string
		self.observed( np.array([float(lines[4+i].split()[1]) for i in range(int(self.stream.m))]) )
		# write lb.dat binaries will find the data
		outfile = open('lb.dat','w');
		us = self.stream.upstream()
		[outfile.write('%e ' % us[i]) for i in range( int(self.stream.m) )]
		outfile.close()
		# set the flag so the module knows it is okay to run the solver (solver requires lb.dat in pwd to run)		
		self.lbWritten = True
		
	def showApprox(self,compare=None):
		'''
		Plots the concentration approximation versus time as triangles; optionally 
		plots another equal length vector as squares. For example, to compare the 
		observed values to the approximation:
		>	approximationInstance.showApprox(approximationInstance.observed())
		
		It is also useful to compare different results with the same parameters:
		>	otisInstance.run(); discusInstance.run()
		>	otisInstance.showApprox(discusInstance.approximated())

		Plots are of concentration versus time for M values. Times are derived from 
		M and dt of the domain.
		'''
		if(self.simRun == False or len(self._approx) == 0):
			raise StandardError, "showApprox() - please run the simulation so that an approximation can be generated"
		M,dt = self.stream.m,self.stream.dt
		times = np.linspace(dt,M*dt,M)
				
		# if compare=None, plot and return
		if(compare == None):
			# clear figure and plot
			clf(); plot(times,self.approximated(),'r^'); show()
		elif(type(compare) != np.ndarray):
			raise TypeError,"showApprox() - argument compare must be an np.array."
		# if the length of approx and observed is the same		
		elif(len(self._approx) == M and len(compare) == M):
			# clear figure and plot
			clf();	plot(times,self._approx,'r^',times,compare,'bs'); show()
		else:
			raise ValueError,"showApprox() - Approximation length not equal to argument length."
		return
			
	def saveFigure(self,filename=None,format='eps'):
		'''
		Saves the current figure as an image file. By default, saves the plot as an .eps file.
		'''
		if(type(filename) != str):
			raise TypeError, "saveFigure() - filename must be a string."
		savefig(filename,format)	
		return

class Otis(Approximation):
	'''
	Otis - The traditional central-difference approximation, solved with a tridiagonal 
	matrix solver. The name Otis comes from the FORTRAN code named OTIS by R. L. Runkel
	(1998)
	'''
	def __init__(self):
		# ensure that parent's init routine is called
		Approximation.__init__(self);
		# set the executable name string	
		self.setBinName('CN')

class Discus(Approximation):	
	'''
	Discus - A semi-Lagrangian method by J Russell Manson. See (Manson, et al 2011)
	for usage in transient storage applications.
	'''
	def __init__(self):
		# ensure that parent's init routine is called
		Approximation.__init__(self);
		# set the executable name string
		self.setBinName('DISCUS')
	
