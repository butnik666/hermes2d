# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_dof_m.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="error (multi)")
data = numpy.loadtxt("conv_dof_s.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="error (single)")

legend()

# initialize new window
pylab.figure()

axis('equal')
pylab.title("Error convergence")
pylab.xlabel("CPU time (s)")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_cpu_m.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="error (multi)")
data = numpy.loadtxt("conv_cpu_s.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="error (single)")
legend()


# finalize
show()
