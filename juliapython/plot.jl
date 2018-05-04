# PLOT example
using PyPlot
using PyCall
@pyimport numpy as np

x = np.random[:randn](200)
xm = mean(x)

plot(x)
title("μ = $xm")
xlabel("time")
