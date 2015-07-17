# So we can use matplotlib
using PyCall
@pyimport matplotlib.pyplot as plt

include("Lenses/BlobbyNIE.jl")

# Load the posterior samples
posterior_sample = readdlm("posterior_sample.txt")

# Size of the 2D array
M = size(posterior_sample)[1]
N = size(posterior_sample)[2]

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

x = linspace(-10,  10, 1001)
y = linspace( 10, -10, 1001)
mag = magnification_image(posterior_sample[1, :], x, y)
plt.imshow(mag, interpolation="nearest", extent=[-10, 10, -10, 10],
					vmin=-10.0, vmax=10.0, cmap="coolwarm")
plt.show()

# http://stackoverflow.com/questions/5666056/matplotlib-extracting-data-from-contour-lines
contour = plt.contour(x, y, mag, [3.0, 4.0])
contour = contour[:collections][:1]
contour = contour[:get_paths]()[1][:vertices]
plt.plot(contour[:,1], contour[:,2])
plt.show()

