using NetCDF
using PyPlot
using PyCall
using Printf

@pyimport cmocean.cm as cm
@pyimport matplotlib.colors as mc
@pyimport matplotlib.cm as cmaps

#path = "/local/kloewer/julsdata/"
path = "/network/aopp/chaos/pred/kloewer/julsdata/ssthr/"

runs = [2]

# LOAD DATA
runpath = path*"run"*Printf.@sprintf("%04d",runs[1])
ncu = NetCDF.open(runpath*"/u.nc")
ncη = NetCDF.open(runpath*"/eta.nc")

u = ncu.vars["u"][:,:,1000:1400]
η = ncη.vars["eta"][:,:,1000:1400]

# permute
u = permutedims(u,[2,1,3])
η = permutedims(η,[2,1,3])

# flip up down
u = cumsum(u[end:-1:1,:,:],dims=1)
η = η[end:-1:1,:,:]


# LOAD OPERATORS
include(runpath*"/scripts/parameters.jl")
include(runpath*"/scripts/src/grid.jl")

## PLOTTING
levs = Array(LinRange(-2.5,2.5,64))
ioff()

for i = 1:400
    fig,axs = subplots(1,1,figsize=(6,3))
    tight_layout(rect=[-0.15,-0.05,1.03,0.98])

    pos1 = axs[:get_position]()
    cax = fig[:add_axes]([pos1[:x1]-0.1,pos1[:y0],0.02,pos1[:y1]-pos1[:y0]])

    norm = mc.Normalize(-2.5,2.5)
    lsource = mc.LightSource(360,5)
    rgb = lsource[:shade_rgb](cmaps.viridis(norm(η[:,:,i])),0.5*u[:,:,i])

    #hand = axs[:pcolormesh](x_T,y_T,η[:,:,i]',vmin=levs[1],vmax=levs[end],cmap="viridis",shading="gouraud")
    hand = axs[:imshow](rgb,interpolation="bilinear",vmin=-2.5,vmax=2.5)

    cb = colorbar(hand,cax=cax,orientation="vertical",extend="both")
    cb[:set_label]("Sea surface height [m]")
    cb[:set_ticks](-2:2)

    axs[:set_xlabel](L"$x$")
    axs[:set_ylabel](L"$y$")

    axs[:set_xticks]([])
    axs[:set_yticks]([])

    axs[:set_title]("Turbulence topography")

    savefig("/home/kloewer/julia/vid/frames3/frame"*Printf.@sprintf("%04d",i)*".png",dpi=300)
    close(fig)
end
