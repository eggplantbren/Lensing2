from pylab import *

b = 0.8
q = 0.7
rc = 0.

def alpha(x, y):
  psi = sqrt(q**2*(x**2 + rc**2) + y**2)
  ax = b/sqrt(1. - q**2)*arctan(sqrt(1. - q**2)*x/(psi + rc))
  ay = b/sqrt(1. - q**2)*arctanh(sqrt(1. - q**2)*y/(psi + q**2*rc))
  return [ax, ay]

def deriv(f, h):
  dfdx = zeros(f.shape)
  dfdy = zeros(f.shape)
  dfdx[:,1:-1:] =  (f[:,2:] - f[:,0:-2])/(2*h)
  dfdy[1:-1:, :] = (f[0:-2, :] - f[2:, :])/(2*h)
  return [dfdx, dfdy]

x = linspace(-2., 2., 1001) + 1E-6
y = linspace(-2., 2., 1001) + 1E-6
[x, y] = meshgrid(x, y)
y = y[::-1, :]
h = x[0, 1] - x[0, 0]

[ax, ay] = alpha(x, y)

# Inside critical curve (only valid for rc=0)
inside = (x/(b/q))**2 + (y/b)**2 < 1

# Terms in divergence of alpha
term1 = deriv(ax, h)[0]
term2 = deriv(ay, h)[1]

# Integral of divergence/(2pi) (inside critical curve) = b^2/q
divergence = term1 + term2
divergence[divergence == 0] = min(divergence[divergence > 0])
print((divergence[inside]).sum()*h**2/(2.*pi), b**2/q)
imshow(log(divergence + 1E-6))
show()

# Compute magnification
Dax = deriv(x - ax, h)
Day = deriv(y - ay, h)
J = Dax[0]*Day[1] - Dax[1]*Day[0]
imshow(-log(abs(J) + 1E-6))
title('Magnification')
show()

