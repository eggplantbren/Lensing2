from pylab import *

def alpha(x, y, shear=0.1, theta_shear=0.):
	cos_theta_shear = cos(theta_shear)
	sin_theta_shear = sin(theta_shear)

	# Go into shear coordinate system
	xx =  x*cos_theta_shear + y*sin_theta_shear;
	yy = -x*sin_theta_shear + y*cos_theta_shear;

	# Calculate external shear
	alphax =  shear*xx;
	alphay = -shear*yy;

	# Add external shear
	ax = alphax*cos_theta_shear - alphay*sin_theta_shear;
	ay = alphax*sin_theta_shear + alphay*cos_theta_shear;

	return [ax, ay]



x = -0.5
y = 0.7

h = 1E-3

[ax1, ay1] = alpha(x + h, y, theta_shear=0.6)
[ax2, ay2] = alpha(x - h, y, theta_shear=0.6)
[ax3, ay3] = alpha(x, y + h, theta_shear=0.6)
[ax4, ay4] = alpha(x, y - h, theta_shear=0.6)

divergence = (ax2 - ax1)/(2*h) + (ay4 - ay3)/(2*h)
print(divergence)

ion()
hold(False)
for theta_shear in linspace(0, pi, 101):
	[ax, ay] = alpha(x, y, theta_shear=theta_shear)

	plot(ax, ay, 'bo')
	hold(True)
	axis('equal')
	title(theta_shear)

	draw()

ioff()
show()

