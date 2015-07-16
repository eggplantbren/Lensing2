# Function to compute the deflection angles
# at position (x, y), for the lens specified by
# the parameters array (a row from sample.txt or posterior_sample.txt)
function alpha(parameters, x, y, cos_theta=cos(parameters[8]),
								 sin_theta=sin(parameters[8]),
								 cos_theta_shear=cos(parameters[10]),
								 sin_theta_shear=sin(parameters[10]))
	# Get NIE parameters
	b = parameters[3]
	q = parameters[4]
	rc = parameters[5]
	xc = parameters[6]
	yc = parameters[7]
	theta = parameters[8]
	shear = parameters[9]
	theta_shear = parameters[10]

	# Stuff derived from b and q
	qq = q;
	if(qq == 1.0)
		qq = 0.99999
	end
	q_term = sqrt(1.0 - qq*qq)
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

# A version that keeps magnification less than 10
function magnification2(parameters, x, y, h=1E-6)
	mag = magnification(parameters, x, y, h)
	if mag > 10.0
		mag = 10.0
	end
	return mag
end

# Current terrible implementation: numerically differentiate
# magnification (which itself involves numerical derivatives).
# Better: actually do second derivatives properly
function magnification_derivs(parameters, x, y, h=1E-6)
	dmdx = (magnification2(parameters, x+h, y) -
				magnification2(parameters, x-h, y))/(2*h)
	dmdy = (magnification2(parameters, x, y+h) -
				magnification2(parameters, x, y-h))/(2*h)
	return [dmdx, dmdy]
end

# Do one step of a leapfrog integrator
# on the position (x, y) in the lens plane,
# treating the magnification as a potential
function update(parameters, pos, vel, dt=1E-3)
	pos = pos + 0.5*dt*vel
	accel = -magnification_derivs(parameters, pos[1], pos[2])
	vel = vel + dt*accel
	pos = pos + 0.5*dt*vel
	return (pos, vel)
end

function dynamics(parameters, pos=[0.0, 0.0], vel = [0.0, 0.0],
						steps=10000, skip=10, dt=1E-3)
	keep = Array(Float64, div(steps, skip), 2)

	for(i in 1:steps)
		(pos, vel) = update(parameters, pos, vel, dt)

		if rem(i, skip) == 0
			keep[div(i, skip), :] = pos
			println(i)
		end
	end

	return keep
end

