# Function to compute the deflection angles
# at position (x, y), for the lens specified by
# the parameters array (a row from sample.txt or posterior_sample.txt)
function alpha(parameters, x, y)
	# Get NIE parameters
	b = posterior_sample[i, 3]
	q = posterior_sample[i, 4]
	rc = posterior_sample[i, 5]
	xc = posterior_sample[i, 6]
	yc = posterior_sample[i, 7]
	theta = posterior_sample[i, 8]
	shear = posterior_sample[i, 9]
	theta_shear = posterior_sample[i, 10]

	cos_theta = cos(theta)
	sin_theta = sin(theta)
	cos_theta_shear = cos(theta_shear)
	sin_theta_shear = sin(theta_shear)

	# Rotate and center
	xx =  (x - xc)*cos_theta + (y - yc)*sin_theta
	yy = -(x - xc)*sin_theta + (y - yc)*cos_theta

	psi = sqrt(qq*qq*(xx*xx + rc*rc) + yy*yy)
	alphax = bb/q_term*atan(q_term*xx/(psi + rc))
	alphay = bb/q_term*atanh(q_term*yy/(psi + qq*qq*rc))

	# Rotate back
	ax = alphax*cos_theta - alphay*sin_theta;
	ay = alphax*sin_theta + alphay*cos_theta;

	# Go into shear coordinate system
	xx =  x*cos_theta_shear + y*sin_theta_shear
	yy = -x*sin_theta_shear + y*cos_theta_shear;

	# Calculate external shear
	alphax = -shear*xx;
	alphay = shear*yy;

	# Add external shear
	ax += alphax*cos_theta_shear - alphay*sin_theta_shear;
	ay += alphax*sin_theta_shear + alphay*cos_theta_shear;

	return [ax, ay]
end

#void BlobbyNIE::alpha(double x, double y, double& ax, double& ay) const
#{

#	// Add blobs
#	if(BlobbyNIE::disable_blobs)
#		return;

#	const vector< vector<double> >& components = blobs.get_components();
#	double rsq, widthsq, Menc;
#	for(size_t i=0; i<components.size(); i++)
#	{
#		rsq = pow(x - components[i][0], 2)
#				+ pow(y - components[i][1], 2);
#		widthsq = pow(components[i][3], 2);

#		if(rsq < widthsq)
#		{
#			Menc = 4.*components[i][2]/widthsq*(0.5*rsq -
#				rsq*rsq/(4*widthsq));
#		}
#		else
#			Menc = components[i][2];
#		ax += Menc*(x - components[i][0])/(M_PI*rsq);
#		ay += Menc*(y - components[i][1])/(M_PI*rsq);
#	}
#}

