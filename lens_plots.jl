# So we can use matplotlib
using PyCall
@pyimport matplotlib.pyplot as plt

include("Lenses/BlobbyNIE.jl")

# Load the posterior samples
posterior_sample = readdlm("posterior_sample.txt")

# Specify a row
which = 1

# Size of the 2D array
M = size(posterior_sample)[1]
N = size(posterior_sample)[2]

# Function to calculate magnification on a grid
function magnification_image(parameters::Array{Float64, 2}, x::Array{Float64, 1}, y::Array{Float64, 1})
	mag = zeros(length(y), length(x))
	for(j in 1:length(y))
		for(i in 1:length(x))
			mag[i, j] = magnification(parameters, x[j], y[i])
		end
	end
	return mag
end

# Function to calculate jacobian on a grid
function jacobian_image(parameters::Array{Float64, 2}, x::Array{Float64, 1}, y::Array{Float64, 1})
	mag = zeros(length(y), length(x))
	for(j in 1:length(y))
		for(i in 1:length(x))
			mag[i, j] = jacobian(parameters, x[j], y[i])
		end
	end
	return mag
end

x = linspace(-10,  10, 3001)
y = linspace( 10, -10, 3001)
mag = magnification_image(posterior_sample[which, :], x, y)

plt.imshow(mag, interpolation="nearest", extent=[-10, 10, -10, 10],
					vmin=-5.0, vmax=5.0, cmap="coolwarm")


params = posterior_sample[which, :]
rays = zeros(100000, 4)
n_rays = 0
for(j in 1:length(x))
	for(i in 1:(length(y)-1))
		f1 = jacobian(params, x[j], y[i])
		f2 = jacobian(params, x[j], y[i+1])
		if(sign(f1) != sign(f2))
			(xs, ys) = fire_ray(params, x[j], y[i])
			rays[n_rays+1, :] = [[x[j], y[i], xs, ys]]
			n_rays += 1
		end
	end
end
rays = rays[1:n_rays, :]

plt.hold(true)
plt.plot(rays[:,1], rays[:,2], "w.", markersize=1)
plt.plot(rays[:,3], rays[:,4], "b.", markersize=1)
plt.axis([-10.0, 10.0, -10.0, 10.0])
plt.show()
writedlm("rays.txt", rays)
plt.show()

