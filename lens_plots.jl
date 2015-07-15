# So we can use matplotlib
using PyCall
@pyimport matplotlib.pyplot as plt

# Load the posterior samples
posterior_sample = readdlm("posterior_sample.txt")
# Size of the 2D array
M = size(posterior_sample)[1]
N = size(posterior_sample)[2]

# Loop over samples
for(i in 1:M)
	# Get NIE parameters
	b = posterior_sample[i, 3]
	q = posterior_sample[i, 4]
	rc = posterior_sample[i, 5]
	xc = posterior_sample[i, 6]
	yc = posterior_sample[i, 7]
	theta = posterior_sample[i, 8]
	shear = posterior_sample[i, 9]
	theta_shear = posterior_sample[i, 10]
end

