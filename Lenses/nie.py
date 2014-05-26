from pylab import *

b = 1.5
q = 0.6
rc = 0.

def alpha(x, y):
  psi = sqrt(q**2*(x**2 + rc**2) + y**2)
  ax = b/sqrt(1. - q**2)*arctan(sqrt(1. - q**2)*x/(psi + rc))
  ay = b/sqrt(1. - q**2)*arctanh(sqrt(1. - q**2)*y/(psi + q**2*rc))
  return [ax, ay]

x = linspace(-2., 2., 1001) + 1E-6
y = linspace(-2., 2., 1001) + 1E-6
[x, y] = meshgrid(x, y)
y = y[::-1, :]
h = x[0, 1] - x[0, 0]

[ax, ay] = alpha(x, y)
psi = sqrt(q**2*(x**2 + rc**2) + y**2)

# Terms in divergence of alpha
term1 = zeros(x.shape)
term2 = zeros(x.shape)
term1[:,1:-1:] = (ax[:,2:] - ax[:,0:-2])/(2*h)
term2[1:-1:, :] = (ay[0:-2, :] - ay[2:, :])/(2*h)

# Integral of divergence/(2pi) (inside radius of psi < b*q) = b^2
divergence = term1 + term2
divergence[divergence == 0] = min(divergence[divergence > 0])
print((divergence*(psi < b*q)).sum()*h**2/(2.*pi))
imshow(log(divergence + 1E-6))
show()

