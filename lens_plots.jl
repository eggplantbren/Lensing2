# So we can use matplotlib
using PyCall
@pyimport matplotlib.pyplot as plt

include("Lenses/BlobbyNIE.jl")

# Load the posterior samples
posterior_sample = readdlm("posterior_sample.txt")
# Size of the 2D array
M = size(posterior_sample)[1]
N = size(posterior_sample)[2]

# Loop over samples
for(i in 1:M)
	println(alpha(posterior_sample[i, :], 0., 0.,))
end

