import approximation

# Note: It is most convenient to run this file in the same directory as the approximation.py module. See README.
# First, define an instance of Discus, the solver that will be used in the example.

solver = approximation.Discus()

# Next, define the model parameters.
# run00 stream5 folder
U=2.81e-01; D=5.89e-02; k1=1.87e-01; k2=1.51e-02
solver.stream.params([U,D,k1,k2])
#You can check and see if you forgot anything with:

# print solver.info()

# Or equivalently:

print solver

# 5.3 Load the domain information along with upstream and observed values.

solver.loadFile("../example/scotland.dat")

# You can check and see that everything is loaded with "print solver", as in 5.2.

print solver

# Run the solver.

solver()

# Or:

# solver.run()

# Plot the approximation. To just see the values calculated:

# solver.showApprox()

# Or, to visually compare it with the observed:

solver.showApprox(solver.observed())

# Get the error norm. Because we have supplied observed data, the solver can calculate the norm.

print solver.getNorm()

'''
Extended Example: This takes off where the above example leaves off, and demonstrates comparing the two solvers. First, the Otis class is instantiated. Next, the stream from the Discus run is copied using copyStream() - this inserts the upstream, domain, and parameters into the new solver. Then the observed data is copied with the observed() routine. Finally, the simulation is run and the simulations can be compared.
'''
'''
solver2 = approximation.Otis()
solver2.copyStream(solver.stream)
solver2.observed(solver.observed())
print solver2
solver2.run()
solver2.showApprox(solver.approximated())
'''
