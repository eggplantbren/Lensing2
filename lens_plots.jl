# So we can use matplotlib
using PyCall
@pyimport matplotlib.pyplot as plt

include("Lenses/BlobbyNIE.jl")

# Load the posterior samples
posterior_sample = readdlm("posterior_sample.txt")
# Size of the 2D array
M = size(posterior_sample)[1]
N = size(posterior_sample)[2]

# Make a grid
x = linspace(-1, 1, 501)
y = linspace(-1, 1, 501)

# Function to calculate magnification on a grid
function magnification_image(parameters, x, y)
	mag = zeros(length(y), length(x))
	for(j in 1:length(y))
		for(i in 1:length(x))
			mag[i, j] = magnification(parameters, x[i], y[j])
		end
	end
	return mag
end

# Loop over samples
for(i in 1:M)
	mag = magnification_image(posterior_sample[i, :], x, y)
	plt.imshow(mag, interpolation="nearest")
	plt.show()
end

