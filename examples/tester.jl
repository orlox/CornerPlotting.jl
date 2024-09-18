using CornerPlotting
using CairoMakie
using Distributions

results = Dict(:a=>rand(500_000,3), :b=>rand(500_000,3), :c=>rand(500_000,3))
results[:weights] = results[:a].*results[:b]

##
set_theme!(CornerPlotting.default_theme())
corner_plot = CornerPlotting.CornerPlot(results,[:a, :b, :c], fraction_1D=0.90)
corner_plot.fig

##
using Turing
using Distributions

@model function test_MCMC()
    a ~ Normal(0.5,3)
    b ~ Normal(0.1, 5)
    x ~ MvNormal([1.0, 2.0], [1.0 0.5;0.5 1.0])
end

num_chains = 5
chain = sample(test_MCMC(), NUTS(), MCMCThreads(), 1000, num_chains)

##
corner_plot = CornerPlotting.CornerPlot(chain,[:a, :b, Symbol("x[1]"), Symbol("x[2]")])
CornerPlotting.plot_extra_1D_distribution(corner_plot, :a, Normal(0.5,3))
CornerPlotting.plot_extra_1D_distribution(corner_plot, :b, Normal(0.1,5))
CornerPlotting.plot_extra_1D_distribution(corner_plot, Symbol("x[1]"), Normal(1.0,1.0))
CornerPlotting.plot_extra_1D_distribution(corner_plot, Symbol("x[2]"), Normal(2.0,1.0))
corner_plot.fig

##
chain = sample(test_MCMC(), NUTS(), 5000)

##
corner_plot = CornerPlotting.CornerPlot(chain,[:a, :b, Symbol("x[1]"), Symbol("x[2]")])
CornerPlotting.plot_extra_1D_distribution(corner_plot, :a, Normal(0.5,3))
CornerPlotting.plot_extra_1D_distribution(corner_plot, :b, Normal(0.1,5))
CornerPlotting.plot_extra_1D_distribution(corner_plot, Symbol("x[1]"), Normal(1.0,1.0))
CornerPlotting.plot_extra_1D_distribution(corner_plot, Symbol("x[2]"), Normal(2.0,1.0))
corner_plot.fig
