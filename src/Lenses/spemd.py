"""
Experimenting to figure out the best parameterisation for the SPEMD
"""

from pylab import *

logr = linspace(-20., 1., 3001)
phi = linspace(0., 2.*pi, 3002)
phi = phi[0:-1]
du = logr[1] - logr[0]
dphi = phi[1] - phi[0]
[logr, phi] = meshgrid(logr, phi)
phi = phi[::-1, :]

# u = log(r) --> du = dr/r
# dx dy = r dr dphi = r^2 du dphi = e^(2u) du dphi 

def spemd_density(x, y, params):
	"""
	params are (b, q, gamma)
	b = intermediate axis
	"""
	b, q, gamma = params

	b_major = b/sqrt(q)

	# kappa(x, y) = q[x1^2+x2^2/q^2 + rc^2]&(-gam)
	coeff = (1.0 - gamma)*b_major**(2.0*gamma)

	return coeff*(x**2 + y**2/q**2 + 1E-7**2)**(-gamma)

x = exp(logr)*cos(phi)
y = exp(logr)*sin(phi)
dA = exp(2*logr)*du*dphi

b = 0.8
q = 0.7

sie = spemd_density(x, y, [b, q, 0.3])
inside = q*x**2 + y**2/q < b**2
print((sie*dA)[inside].sum(), pi*b**2)
imshow(log(sie))
show()

