from pylab import *

def g(x):
  result = zeros(x.shape)
  pos = abs(x) < 1.
  result[pos] = 1. - x[pos]**2
  return result

# Make a cartesian grid
x = linspace(-5., 5., 256)
[x, y] = meshgrid(x, x)
y = y[::-1, :]

# Calculate the plane polar coordinates
r = sqrt(x**2 + y**2)
phi = arctan2(y, x + 1E-6)
phi[phi < 0.] += 2*pi

# Initialise image to zero
f = zeros(x.shape)

imshow(g(r**2 - phi**2.))
show()

