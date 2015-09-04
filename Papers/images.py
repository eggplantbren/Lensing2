from pylab import *
rc("font", size=16, family="serif", serif="Computer Sans")
rc("text", usetex=True)

img = loadtxt("../Data/mock_image.txt")
truth = loadtxt("../Data/mock_truth.txt")
metadata = loadtxt("../Data/mock_metadata.txt")

# Convert x and y to pixel coordinates for overplotting
dx = (metadata[3] - metadata[2])/metadata[1]
dy = (metadata[5] - metadata[4])/metadata[0]
x_substructures = truth[19:29]
y_substructures = truth[29:39]
x_substructures = x_substructures[x_substructures != 0.]
y_substructures = y_substructures[y_substructures != 0.]
x_substructures = (x_substructures - metadata[2])/dx - 0.5
y_substructures = (metadata[5] - y_substructures)/dy - 0.5
x_nie = (truth[5] - metadata[2])/dx - 0.5
y_nie = (metadata[5] - truth[6])/dy - 0.5

imshow(img, interpolation="nearest", cmap="Oranges")
hold(True)
# Plot center of NIE
plot(x_nie, y_nie, 'wo', markersize=15)
# Substructures
plot(x_substructures, y_substructures, 'g*', markersize=15)
xlim([-0.5, metadata[1] - 0.5])
ylim([metadata[0] - 0.5, -0.5])

gca().set_xticklabels("")
gca().set_yticklabels("")
title("Simulated Dataset")
savefig("image1.pdf", bbox_inches="tight")
show()

img = loadtxt("../Data/harder_image.txt")
truth = loadtxt("../Data/harder_truth.txt")
metadata = loadtxt("../Data/mock_metadata.txt")

# Convert x and y to pixel coordinates for overplotting
dx = (metadata[3] - metadata[2])/metadata[1]
dy = (metadata[5] - metadata[4])/metadata[0]
x_substructures = truth[19:29]
y_substructures = truth[29:39]
x_substructures = x_substructures[x_substructures != 0.]
y_substructures = y_substructures[y_substructures != 0.]
x_substructures = (x_substructures - metadata[2])/dx - 0.5
y_substructures = (metadata[5] - y_substructures)/dy - 0.5
x_nie = (truth[5] - metadata[2])/dx - 0.5
y_nie = (metadata[5] - truth[6])/dy - 0.5

imshow(img, interpolation="nearest", cmap="Oranges")
hold(True)
# Plot center of NIE
plot(x_nie, y_nie, 'wo', markersize=15)
# Substructures
plot(x_substructures, y_substructures, 'g*', markersize=15)
xlim([-0.5, metadata[1] - 0.5])
ylim([metadata[0] - 0.5, -0.5])
gca().set_xticklabels("")
gca().set_yticklabels("")
title("Simulated Dataset")
savefig("image2.pdf", bbox_inches="tight")
show()

img = loadtxt("../Data/horseshoe_image.txt")
imshow(img, interpolation="nearest", cmap="Oranges")
gca().set_xticks([-0.5, 0.5*img.shape[1], img.shape[1] - 0.5])
gca().set_yticks([-0.5, 0.5*img.shape[0], img.shape[0] - 0.5])
gca().set_xticklabels([-7.66, 0, 7.66])
gca().set_yticklabels([-7.66, 0, 7.66])
xlabel("$x$ (arcseconds)")
ylabel("$y$ (arcseconds)")
title("Cosmic Horseshoe")
savefig("image3.pdf", bbox_inches="tight")
show()

