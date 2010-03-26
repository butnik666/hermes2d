# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_dof_newt.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="Newton solve on coarse mesh")
data = numpy.loadtxt("conv_dof_proj.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="no Newton solve on coarse mesh")

legend()

# finalize
show()
