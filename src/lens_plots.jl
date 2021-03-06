# So we can use matplotlib
using PyCall
@pyimport matplotlib.pyplot as plt

include("Lenses/BlobbyNIE.jl")

# Load the posterior samples
posterior_sample = readdlm("posterior_sample.txt")
metadata = readdlm("Data/mock_metadata.txt")

plt.figure(figsize=(8, 8))
for(which in 1:3)
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

#	x = linspace(metadata[3], metadata[4], 1001)
#	y = linspace(metadata[5], metadata[6], 1001)
#	mag = magnification_image(posterior_sample[which, :], x, y)
#	plt.imshow(mag, interpolation="nearest", extent=[-10, 10, -10, 10],
#						vmin=-5.0, vmax=5.0, cmap="coolwarm")

	params = posterior_sample[which, :]
	x_substructures = params[20:29]
	y_substructures = params[30:39]
	m_substructures = params[40:49]
	x_substructures = x_substructures[m_substructures .!= 0.]
	y_substructures = y_substructures[m_substructures .!= 0.]
	x_nie, y_nie = params[6], params[7]

#	rays = zeros(100000, 4)
#	n_rays = 0
#	for(j in 1:length(x))
#		for(i in 1:(length(y)-1))
#			f1 = jacobian(params, x[j], y[i])
#			f2 = jacobian(params, x[j], y[i+1])
#			if(sign(f1) != sign(f2))
#				(xs, ys) = fire_ray(params, x[j], y[i])
#				rays[n_rays+1, :] = [[x[j], y[i], xs, ys]]
#				n_rays += 1
#			end
#		end
#	end
#	for(i in 1:length(y))
#		for(j in 1:(length(x)-1))
#			f1 = jacobian(params, x[j], y[i])
#			f2 = jacobian(params, x[j+1], y[i])
#			if(sign(f1) != sign(f2))
#				(xs, ys) = fire_ray(params, x[j], y[i])
#				rays[n_rays+1, :] = [[x[j], y[i], xs, ys]]
#				n_rays += 1
#			end
#		end
#	end

#	rays = rays[1:n_rays, :]

##	plt.hold(true)
##	plt.plot(rays[:,1], rays[:,2], "w.", markersize=1)
##	plt.plot(rays[:,3], rays[:,4], "b.", markersize=1)
##	plt.axis([-10.0, 10.0, -10.0, 10.0])
##	writedlm("rays.txt", rays)
##	plt.show()

	src = vec(posterior_sample[which, 469:468 + metadata[1]*metadata[2]*metadata[8]^2])
	img = vec(posterior_sample[which, 469 + 2*metadata[1]*metadata[2]*metadata[8]^2:(size(posterior_sample)[2]-2)])
	src = transpose(reshape(src, int(metadata[2]*metadata[8]), int(metadata[1]*metadata[8])))
	img = transpose(reshape(img, int(metadata[2]), int(metadata[1])))

	plt.rc("font", size=14, family="serif", serif="Computer Sans")
	plt.rc("text", usetex=true)
	plt.subplot(3, 3, 3*which - 2)
	plt.imshow(src, interpolation="nearest", cmap="Oranges", extent=metadata[3:6])
	plt.hold(true)
#	plt.plot(rays[:,3], rays[:,4], "k.", markersize=1)
	if(which == 1)
		plt.title("Source")
	end
	plt.axis(0.5*metadata[3:6])
	plt.gca()[:set_xticks]([])
	plt.gca()[:set_yticks]([])

	plt.subplot(3, 3, 3*which - 1)
	f = density(posterior_sample[which, :], metadata[3], metadata[4], metadata[5], metadata[6])
	plt.imshow(f, interpolation="nearest", cmap="Oranges", extent=metadata[3:6])
	plt.gca()[:set_xticks]([])
	plt.gca()[:set_yticks]([])
	if(which == 1)
		plt.title("Density in substructures")
	end

	plt.subplot(3, 3, 3*which)
	plt.imshow(img, interpolation="nearest", cmap="Oranges", extent=metadata[3:6])
	plt.hold(true)
#	plt.plot(rays[:,1], rays[:,2], "k.", markersize=1)

	# Plot center of NIE
	plt.plot(x_nie, y_nie, "wo")
	# Substructures
	plt.plot(x_substructures, y_substructures, "g*", alpha=0.5)
	plt.axis(metadata[3:6])

	if(which == 1)
		plt.title("Image")
	end
	plt.gca()[:set_xticks]([])
	plt.gca()[:set_yticks]([])
end

plt.savefig("sources.pdf", bbox_inches="tight")
plt.show()

