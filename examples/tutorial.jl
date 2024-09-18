#=
# Tutorial

Here we show two use cases for this package. Corner plots can be produced
from either a dictionary or an MCMCChains instance. When using a dictionary,
each value can be a vector (representing a single chain) or an array, where
each column represents a chain. Below we create some mock values uniformly
samples between zero and one, and add a weight. If the weight is not included,
the package will weigh all samples equally.
=#

using CornerPlotting
using CairoMakie
using Distributions

results = Dict(:a=>rand(500_000,3), :b=>rand(500_000,3), :c=>rand(500_000,3))
results[:weights] = results[:a].*results[:b]

##
#=
We can the produce the corner plot. We use Makie as the plotting backend,
and we provide a default theme.
=#

set_theme!(CornerPlotting.default_theme())
corner_plot = CornerPlotting.CornerPlot(results,[:a, :b, :c])
corner_plot.fig

##
#=
The plotter also works with MCMCChains produced with Turing.jl.
The example below samples a few variables (two of them correlated)
and then includes in the corned plot the sampled variables to verify
we obtain the expected result. Model is sampled with 4 chains here.
=#
using Turing
using Distributions

@model function test_MCMC()
    a ~ Normal(0.5,3)
    b ~ Normal(0.1, 5)
    x ~ MvNormal([1.0, 2.0], [1.0 0.5;0.5 1.0])
end

num_chains = 4
chain = sample(test_MCMC(), NUTS(), MCMCThreads(), 1000, num_chains)

##
corner_plot = CornerPlotting.CornerPlot(chain,[:a, :b, Symbol("x[1]"), Symbol("x[2]")])
CornerPlotting.plot_extra_1D_distribution(corner_plot, :a, Normal(0.5,3))
CornerPlotting.plot_extra_1D_distribution(corner_plot, :b, Normal(0.1,5))
CornerPlotting.plot_extra_1D_distribution(corner_plot, Symbol("x[1]"), Normal(1.0,1.0))
CornerPlotting.plot_extra_1D_distribution(corner_plot, Symbol("x[2]"), Normal(2.0,1.0))
corner_plot.fig

##
#=
The same can be achieved when there is only one chain available
=#
chain = sample(test_MCMC(), NUTS(), 5000)

##
corner_plot = CornerPlotting.CornerPlot(chain,[:a, :b, Symbol("x[1]"), Symbol("x[2]")])
CornerPlotting.plot_extra_1D_distribution(corner_plot, :a, Normal(0.5,3))
CornerPlotting.plot_extra_1D_distribution(corner_plot, :b, Normal(0.1,5))
CornerPlotting.plot_extra_1D_distribution(corner_plot, Symbol("x[1]"), Normal(1.0,1.0))
CornerPlotting.plot_extra_1D_distribution(corner_plot, Symbol("x[2]"), Normal(2.0,1.0))
corner_plot.fig

##
#=
Finally, all axes of the plot can be accesed from the CornerPlot struct.
This allows us to add arbitrary content.
=#
xvals = LinRange(-5.0, 5.0, 100)
yvals = LinRange(0.0, 0.1, 100)

axis_1d = corner_plot.distributions_1d[:a]
lines!(axis_1d, xvals, yvals)
axis_2d = corner_plot.distributions_2d[:a][Symbol("x[1]")]
scatter!(axis_2d, [0.0], [1.0], color=:red, markersize=30)

corner_plot.fig

##
#=
Beware that in the example above you can access `corner_plot.distributions_2d[:a][:b]`,
but `corner_plot.distributions_2d[:b][:a]` is not defined. The first reference corresponds
to the x-axis of the 2D marginalized distribution, while the second one corresponds to
the y-axis, and a matching plot must be present in the figure.
=#