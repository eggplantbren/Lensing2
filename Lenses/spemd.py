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
	params are (b0, q, gamma)
	b0 = semi-major axis
	"""
	b, q, gamma = params

	# kappa(x, y) = q[x1^2+x2^2/q^2 + rc^2]&(-gam)
	coeff = 0.5*b**(2*gamma)*(2.-2.*gamma)

	return coeff*(x**2 + y**2/q**2 + 1E-7**2)**(-gamma)

x = exp(logr)*cos(phi)
y = exp(logr)*sin(phi)
dA = exp(2*logr)*du*dphi

sie = spemd_density(x, y, [1.3, 1., 0.5])
within = x**2 + y**2 < 1.3**2
print((sie*dA)[within].sum())
imshow(log(sie))
show()

