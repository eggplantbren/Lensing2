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

