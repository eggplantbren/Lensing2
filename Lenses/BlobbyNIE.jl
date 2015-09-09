# Function to compute the deflection angles
# at position (x, y), for the lens specified by
# the parameters array (a row from sample.txt or posterior_sample.txt)
function alpha(parameters::Array{Float64, 2}, x::Float64, y::Float64, cos_theta::Float64=cos(parameters[8]),
								 sin_theta::Float64=sin(parameters[8]),
								 cos_theta_shear::Float64=cos(parameters[10]),
								 sin_theta_shear::Float64=sin(parameters[10]))
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

function fire_ray(parameters::Array{Float64, 2}, x::Float64, y::Float64)
	(ax, ay) = alpha(parameters, x, y)
	xs = x - ax
	ys = y - ay
	return [xs, ys]
end

# Compute the jacobian of the lens-->source plane mapping at position (x, y).
# h is the spacing for numerical differentiation
function jacobian(parameters::Array{Float64, 2}, x::Float64, y::Float64, h::Float64=1E-7)
	(xs1, ys1) = fire_ray(parameters, x + h, y)
	(xs2, ys2) = fire_ray(parameters, x - h, y)

	J11 = (xs1 - xs2)/(2*h)
	J12 = (ys1 - ys2)/(2*h)

	(xs1, ys1) = fire_ray(parameters, x, y + h)
	(xs2, ys2) = fire_ray(parameters, x, y - h)

	J21 = (xs1 - xs2)/(2*h)
	J22 = (ys1 - ys2)/(2*h)

	return(J11*J22 - J12*J21)
end


# Compute the magnification (in magnitudes) at position (x, y).
# h is the spacing for numerical differentiation
function magnification(parameters::Array{Float64, 2}, x::Float64, y::Float64, h::Float64=1E-6)
	return -2.5*log10(abs(jacobian(parameters, x, y, h)))
end

function density(parameters::Array{Float64, 2},
					x_min, x_max, y_min, y_max; n=1001)
	assert(size(parameters)[1] == 1)

	x = linspace(x_min, x_max, n)
	y = linspace(y_max, y_min, n)
	f = zeros(n, n)

	for(k in 1:parameters[19])
		xc, yc, mass, width = parameters[19+k], parameters[29+k], parameters[39+k], parameters[49+k]
		widthsq = width^2
		C = mass*2.0/pi/widthsq

		for(j in 1:n)
			for(i in 1:n)
				rsq = (x[i] - xc)^2 + (y[j] - yc)^2
				if(rsq <= widthsq)
					f[i, j] += C*(1.0 - rsq/widthsq)
				end
			end
		end
	end

	return f
end

# Input: parameters of a blob (length four)
function blob_mass_within(parameters::Array{Float64, 2},
							x_min, x_max, y_min, y_max; n=1001)
	assert(size(parameters) == (1, 4))
	xc, yc, mass, width = parameters
	widthsq = width^2
	C = mass*2.0/pi/widthsq

	x = linspace(x_min, x_max, n)
	y = linspace(y_max, y_min, n)
	f = zeros(n, n)
	for(j in 1:n)
		for(i in 1:n)
			rsq = (x[i] - xc)^2 + (y[j] - yc)^2
			if(rsq <= widthsq)
				f[i, j] += C*(1.0 - rsq/widthsq)
			end
		end
	end

	return (x[2] - x[1])*(y[1] - y[2])*sum(f)
end

params = zeros(1, 4)
params[1], params[2], params[3], params[4] = 0.0, 0.0, 1.0, 0.5
println(blob_mass_within(params, -0.5, 0.5, -0.5, 0.5))


