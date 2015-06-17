from pylab import *
rc("font", size=16, family="serif", serif="Computer Sans")
rc("text", usetex=True)

img = loadtxt("../Data/mock_image.txt")
imshow(img, interpolation="nearest", cmap="Oranges")
gca().set_xticklabels("")
gca().set_yticklabels("")
title("Simulated Dataset")
savefig("simulated_image.pdf", bbox_inches="tight")
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
savefig("horseshoe_image.pdf", bbox_inches="tight")
show()

