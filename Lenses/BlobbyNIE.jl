# Function to compute the deflection angles
# at position (x, y), for the lens specified by
# the parameters array (a row from sample.txt or posterior_sample.txt)
function alpha(parameters, x, y, cos_theta=cos(parameters[8]),
								 sin_theta=sin(parameters[8]),
								 cos_theta_shear=cos(parameters[10]),
								 sin_theta_shear=sin(parameters[10]))
	# Get NIE parameters
	b = posterior_sample[3]
	q = posterior_sample[4]
	rc = posterior_sample[5]
	xc = posterior_sample[6]
	yc = posterior_sample[7]
	theta = posterior_sample[8]
	shear = posterior_sample[9]
	theta_shear = posterior_sample[10]

	# Stuff derived from b and q
	qq = q;
	if(qq == 1.)
		qq = 0.99999
	end
	q_term = sqrt(1. - qq*qq)
	bb = b*sqrt(qq) # Minor axis

	# Rotate and center
	xx =  (x - xc)*cos_theta + (y - yc)*sin_theta
	yy = -(x - xc)*sin_theta + (y - yc)*cos_theta

	psi = sqrt(qq*qq*(xx*xx + rc*rc) + yy*yy)
	alphax = bb/q_term*atan(q_term*xx/(psi + rc))
	alphay = bb/q_term*atanh(q_term*yy/(psi + qq*qq*rc))

	# Rotate back
	ax = alphax*cos_theta - alphay*sin_theta
	ay = alphax*sin_theta + alphay*cos_theta

	# Go into shear coordinate system
	xx =  x*cos_theta_shear + y*sin_theta_shear
	yy = -x*sin_theta_shear + y*cos_theta_shear

	# Calculate external shear
	alphax = -shear*xx
	alphay = shear*yy

	# Add external shear
	ax += alphax*cos_theta_shear - alphay*sin_theta_shear
	ay += alphax*sin_theta_shear + alphay*cos_theta_shear

	# Add the blobs
	num_blobs = int64(parameters[19])
	for(i in 1:num_blobs)
		rsq = (x - parameters[20 + 4*(i-1)])^2 + (y - parameters[30 + 4*(i-1)])^2
		widthsq = parameters[50 + 4*(i-1)]^2
		Menc = parameters[40 + 4*(i-1)]

		if(rsq < widthsq)
			Menc = 4.*Menc/widthsq*(0.5*rsq - rsq*rsq/(4*widthsq));
		end
		ax += Menc*(x - parameters[20 + 4*(i-1)])/(pi*rsq);
		ay += Menc*(y - parameters[30 + 4*(i-1)])/(pi*rsq);
	end

	return [ax, ay]
end

function fire_ray(parameters, x, y)
	(ax, ay) = alpha(parameters, x, y)
	xs = x - ax
	ys = y - ay
	return [xs, ys]
end

# Compute the magnification at position (x, y).
# h is the spacing for numerical differentiation
function magnification(parameters, x, y, h=1E-6)
	(xs1, ys1) = fire_ray(parameters, x + h, y)
	(xs2, ys2) = fire_ray(parameters, x - h, y)

	J11 = (xs1 - xs2)/(2*h)
	J12 = (ys1 - ys2)/(2*h)

	(xs1, ys1) = fire_ray(parameters, x, y + h)
	(xs2, ys2) = fire_ray(parameters, x, y - h)

	J21 = (xs1 - xs2)/(2*h)
	J22 = (ys1 - ys2)/(2*h)

	return -2.5*log10(abs(J11*J22 - J12*J21))
end

