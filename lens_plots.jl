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
			mag[i, j] = magnification(parameters, x[i], y[j])
		end
	end
	return mag
end

x = linspace(-10,  10, 3001)
y = linspace( 10, -10, 3001)
mag = magnification_image(posterior_sample[which, :], x, y)

## Use matplotlib's contour function to get the critical curve
## http://stackoverflow.com/questions/5666056/matplotlib-extracting-data-from-contour-lines
#contour = plt.contour(x, y, mag, [5.5, 5.51])
#plt.clf()
#plt.imshow(mag, interpolation="nearest", extent=[-10, 10, -10, 10],
#					vmin=-7.0, vmax=7.0, cmap="coolwarm")
#contour = contour[:collections][:1]
#contour = contour[:get_paths]()[1][:vertices]
#plt.plot(contour[:,1], contour[:,2], "k", linewidth=2)
#plt.show()

## Fire the contour through to the source plane
#xs = contour[:,1]
#ys = contour[:,2]
#for(i in 1:size(xs)[1])
#	(xs[i], ys[i]) = fire_ray(posterior_sample[which, :],
#									contour[i, 1], contour[i, 2])
#end

plt.imshow(mag, interpolation="nearest", extent=[-10, 10, -10, 10],
					vmin=-7.0, vmax=7.0, cmap="coolwarm")
plt.show()

xx = zeros(size(y)[1], size(x)[1])
yy = zeros(size(y)[1], size(x)[1])
for(j in 1:size(yy)[1])
	for(i in 1:size(xx)[1])
		xx[i, j] = x[j]
		yy[i, j] = y[i]
	end
end

xx = xx[mag .> 5.5]
yy = yy[mag .> 5.5]
xs = xx
ys = yy
for(i in 1:size(xs)[1])
	(xs[i], ys[i]) = fire_ray(posterior_sample[which, :],
									xx[i], yy[i])
end



plt.figure(figsize=(8, 8))
plt.plot(xs, ys, "k.", markersize=1)
plt.xlabel("x_s")
plt.ylabel("y_s")
rays = zeros(size(xs)[1], 4)
rays[:,1] = xx
rays[:,2] = yy
rays[:,3] = xs
rays[:,4] = ys
writedlm("rays.txt", rays)
plt.show()

