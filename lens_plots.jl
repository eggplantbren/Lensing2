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
plt.ion()
plt.hold(true)
plt.imshow(mag, interpolation="nearest", extent=[-10, 10, -10, 10],
					vmin=-10.0, vmax=10.0)
pos = [5.0, 7.0]
for(i in 1:1000)
	(pos[1], pos[2]) = move_along_contour(posterior_sample[1, :],
							pos[1], pos[2], 0.03)
	plt.plot(pos[1], pos[2], "ro")
	plt.draw()
end
plt.show()



## Loop over samples
## Make a grid
#x = linspace(-10,  10, 501)
#y = linspace( 10, -10, 501)
#plt.ion()
#plt.hold(false)
#for(i in 1:M)
#	mag = magnification_image(posterior_sample[i, :], x, y)
#	plt.imshow(mag, interpolation="nearest", vmin=-10.0, vmax=10.0)
#	plt.title(i)
#	plt.draw()
#end
#plt.ioff()
#plt.show()

