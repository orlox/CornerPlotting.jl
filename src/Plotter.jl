using StatsBase
using CairoMakie
using StyledStrings
using Distributions



"""
    struct CornerPlot(fig, ranges, distributions_1d, distributions_2d)

Stores the elements of a CornerPlot to allow easy access to figure elements.
"fig" contains the Makie figure, "ranges" corresponds to the ranges used for each
variable in each plot, as well as to determine credible intervals. The Axes with
1D marginalized distributions are stored in the dictionary "distributions_1d",
using as keys the name of each variable. Similarly, "distributions_2d" provides
the Axes with 2D marginalized distributions as a dictionary of dictionaries. 

"""
struct CornerPlot
    fig
    ranges
    distributions_1d
    distributions_2d
end

"""
    CornerPlot(results, names; ...)

Constructs a corner plot from the provided "results", using the variables given
in the "names" vector. 

# Arguments:
- results: Array containing the samples to be plotted. This can be either of type MCMCChains
or a Dictionary containing vectors with values for an individual chain, or arrays for multiple
chains (in this case each column represents a chain). 
- names: Vector of symbols containing the key needed to access each result from `results`.
The corner plot will only include the values specified in `names`.
- labels: Dictionary of strings containing the labels that should be used for each variable.
- ranges: Dictionary of two element vectors, containing the ranges that will be used for each
variable of the plot. If ranges are not provided for a variable these are determined based on the
`quantile_for_range` option.
- scaling: Dictionary containing scaling factors for variables. For any `name` in `names` that
is also a key of `scaling`, all values are divided by `scaling[name]`.
- fig: The Makie figure used for the plot. If not provided it is created.
- quantile_for_range: If ranges are not specified for an axis, then they are set to be between
the quantiles `quantile_for_range` and `1- quantile_for_range`. This is done using weighted
quantiles if `use_weights=true`
- use_weights: If true, then `results[:weights]` is expected to be defined to provide weights
for each sample.
- fraction_1d: Fraction of samples contained in the shown 1D credible intervals. Credible intervals
are determined using highest density intervals. By default 90% credible intervals are shown.
- fractions_2d: Similar to `fraction_1d`, but used to determine the contours in the 2D marginalized
distributions. Values are provided as a Vector of fractions.
- show_CIs: If true, credible intervals are shown in the corner plot.
- nbins: Number of bins in each axis used to plot the heatmaps and the 1D marginalized distributions
- nbins_contour: Number of bins used to plot the contours in the 2D marginalized distributions.
using `nbins_contour<nbins` allows for smoother contous.
- axis_size: The Makie axis will be set to have width and height equal to this value.
- supertitle: The title of the plot
- supertitlefontsize: Fontsize of the title 

# Output:
Returns an instance of CornerPlot
"""
function CornerPlot(results, names::Vector{Symbol};
        labels=nothing, ranges=Dict(), scaling=Dict(),
        fig=Figure(), quantile_for_range=0.01,
        use_weights = true, fraction_1D=0.9, fractions_2D=[0.9], 
        show_CIs=true, nbins=100, nbins_contour=20,
        axis_size=100, supertitle=nothing, supertitlefontsize=25)

    num_col = length(names)
    
    distributions_2d = Dict() # this will store the 2D distributions in a dictionary of dictionaries
    for name in names
        if name ∉ keys(results)
            throw(ArgumentError("$name is not a valid key"))
        end
        distributions_2d[name] = Dict()
    end
    if isnothing(labels)
        labels = Dict()
        for name in names
            labels[name] = String(name)
        end
    end
    if :weights ∈ keys(results) && use_weights
        sample_weights = results[:weights]
    else
        # just create an array of ones with the needed size
        # not the most efficient but avoids having to write special cases later
        sample_weights = ones(size(results[names[1]]))
    end
    # detect ranges
    for name in names
        if name ∉ keys(ranges)
            ranges[name] = quantile(vec(results[name]),weights(sample_weights),
                                    [quantile_for_range, 1-quantile_for_range])
            if name ∈ keys(scaling)
                ranges[name] /= scaling[name]
            end
        end
    end
 
    # Create 2D density plots
    for ii in 1:num_col-1         # ii is the x-coord param
        name_x = names[ii]
        for  jj in ii+1:num_col   # jj is the y-coord param
            name_y = names[jj]
            axis = Axis(fig[jj+1,ii], xlabel=labels[name_x], ylabel=labels[name_y], height=axis_size, width=axis_size)
            distributions_2d[name_x][name_y] = axis
            if name_x ∉ keys(scaling)
                values_x = results[name_x]
            else
                values_x = results[name_x]./scaling[name_x]
            end
            if name_y ∉ keys(scaling)
                values_y = results[name_y]
            else
                values_y = results[name_y]./scaling[name_y]
            end

            plot_2D_density(axis, name_x, name_y, values_x, ranges[name_x], values_y, ranges[name_y],
                sample_weights, fractions_2D, nbins_heatmap=nbins, nbins_contour=nbins_contour)
            if ii>1
                hideydecorations!(axis, ticks=false, minorticks=false)
            end
            if jj!=num_col
                hidexdecorations!(axis,ticks=false, minorticks=false)
            end         
            xlims!(axis, (ranges[name_x][1], ranges[name_x][2]))
            ylims!(axis, (ranges[name_y][1], ranges[name_y][2]))
        end  
    end
    # Create 1D PDFs along the diagonal
    distributions_1d = Dict() # we will add the distributions with their names here
    latex_bounds_array = Array{AbstractString}(undef, num_col)
    for ii in 1:num_col
        name_x = names[ii]
        axis = Axis(fig[ii+1,ii], xlabel=labels[name_x], height=axis_size, width=axis_size)
        distributions_1d[name_x] = axis
        if name_x ∉ keys(scaling)
            values_x = results[name_x]
        else
            values_x = results[name_x]./scaling[name_x]
        end
        (xmin, xmode, xmax), x, h, y, dx, frac_lost =
            plot_compound_1D_density(axis, name_x, values_x, ranges[name_x], sample_weights, fraction_1D, nbins)

        # Remove labels on the diagonals, except xlabel on the bottom right
        hideydecorations!(axis)
        if ii !=num_col
            hidexdecorations!(axis,ticks=false, minorticks=false)
        end
        # Configure confidence intervals
        str_xmode = round(xmode, sigdigits=3)
        str_upp = round(xmax-xmode, sigdigits=3)
        str_low = round(xmode-xmin, sigdigits=3)
        latex_bounds = L"%$(str_xmode)^{+%$(str_upp)}_{-%$(str_low)}"
        latex_bounds_array[ii] = latex_bounds
        println(labels[name_x]*"="*latex_bounds)
        if show_CIs
            Label(fig[ii,ii], latex_bounds, valign=:bottom)
            axis = Axis(fig[ii,ii])
            hidedecorations!(axis)
            hidespines!(axis)
        end
    end     

    if !isnothing(supertitle)
        Label(fig, text=supertitle, fontsize=supertitlefontsize)
    end
    resize_to_layout!(fig)
    
    return CornerPlot(fig, ranges, distributions_1d, distributions_2d)
end

"""
    get_bounds_for_fractions(h, fractions)

Calculate the probability values that corresponds to a given set of highest density intervals.

# Arguments:
- h: `h` is expected to be a vector or array containing the values of a PDF
within a regular grid. This can be either a 1D vector for 1D marginalized distributions,
or a 2D vector for the 2D marginalized distributions. Can also be used for higher dimensional
cases though.
- fractions: Vector containing the individual fractions that are contained within the HDIs.
for example, `fractions=[0.5,0.9]` means that the value of the PDF corresponding to the 50%
and 90% HDIs will be computed.

# Output:
- bounds: Vector containing the values of the PDF corresponding to the HDIs that contain
the fraction of the samples given by `fractions`.
"""
function get_bounds_for_fractions(h, fractions)
    integral = sum(h.weights)
    bounds =zeros(length(fractions))
    for jj in eachindex(fractions)
        fraction = fractions[jj]
        minbound = 0
        maxbound = maximum(h.weights)
        newbound = 0
        # find HDI 
        for ii in 1:15
            newbound = 0.5*(minbound+maxbound)
            integral2 = sum(h.weights[h.weights.>newbound])
            newfraction = integral2/integral
            if newfraction>fraction
                minbound = newbound
            else
                maxbound = newbound
            end
        end
        bounds[jj] = newbound
    end 
    return sort(bounds)
end

"""
    plot_2D_density(axis, values, range, chain_weights, fraction_1D, nbins; color, linewidth)

Creates a 1D marginalized distribution plot in `axis` from the sample values given in `values`.
The plotted line is normalized, such that it corresponds to the PDF followed by the samples.

# Arguments:
- axis: The Axis on which the distribution will be plotted. If set to `nothing`, no plot will be made,
but the values used for the plot will still be returned.
- name_x: Symbol containing the name for the variable on the x-axis. Only used to report if too
many samples are outside the range to determine a good HDI.
- name_y: Same as `name_x` but for the y-axis.
- values_x: Values for each sample corresponding to the variable `name_x`.
Either a vector or a matrix of numbers, corresponding to a single chain or multiple chains
respectively (in the latter case, columns represent individual samples). The heatmap will include
values from all chains.
- values_y: Same as `values_x` but for the y-axis.
- range_x: Vector with two values containing the range of the plot for the x-axis variable. It also corresponds to the range
that will be used to determine the HDI. The HDI does correct for values outside the range, assuming
that these do not fall within the highest density region.
- range_y: Same as `range_x` but for the y-axis.
- sample_weights: weights used for each sample, should have the same size as `values`.
- fractions: Contours will be plotted corresponding the HDIs containing each fraction contained in fractions.
- nbins: Number of bins used for the heatmap.
- nbins_contours: Number of bins used to determine the contours for the HDIs.

# Output:
- x_hm: Bin centers used in the heatmap in the x-axis
- y_hm: Bin centers used in the heatmap in the y-axis
- z_hm: Value of the PDF estimated at each bin of the heatmap
- dx_hm: Size of bins used for heatmap in the x-axis
- dy_hm: Size of bins used for heatmap in the y-axis
- x_ct: Same as `x_hm` but for the contours.
- y_ct: Same as `y_hm` but for the contours.
- z_ct: Same as `z_hm` but for the contours.
- dx_ct: Same as `dx_hm` but for the contours.
- dy_ct: Same as `dy_hm` but for the contours.
- bounds: probablity values corresponding to the HDIs
- frac_lost: Fraction of samples (including weights) that is outside of the ranges
"""
function plot_2D_density(axis, name_x, name_y, values_x, range_x, values_y, range_y, sample_weights, fractions; nbins_heatmap, nbins_contour=-1)

    filter = (values_x .> range_x[1]) .&& (values_x .< range_x[2]) .&& (values_y .> range_y[1]) .&& (values_y .> range_y[1])
    total_weight = sum(sample_weights)
    missing_weight = sum(sample_weights[.!filter])
    frac_lost = missing_weight/total_weight

    edges = (LinRange(range_x[1], range_x[2], nbins_heatmap+1),
             LinRange(range_y[1], range_y[2], nbins_heatmap+1))
    h_hm = fit(Histogram, (values_x[filter], values_y[filter]), weights(sample_weights[filter]), edges)#, nbins=nbins_heatmap) 
    x_hm = (h_hm.edges[1][2:end] .+ h_hm.edges[1][1:end-1])./2
    y_hm = (h_hm.edges[2][2:end] .+ h_hm.edges[2][1:end-1])./2
    dx_hm = x_hm[2]-x_hm[1]
    dy_hm = x_hm[2]-x_hm[1]
    z_hm = h_hm.weights/(total_weight*dx_hm*dy_hm)
    if !isnothing(axis)
        heatmap!(axis, x_hm, y_hm, z_hm, colormap=:dense, colorrange=(0, maximum(z_hm)))
    end

    #we correct for samples outside the range, assuming they are
    #not part of the HDI. A basic check is to verify we dont have too
    #many samples outside.
    if 1 - maximum(fractions) < frac_lost
        println(styled"{red: range used is insufficient to determine 2D HDI of $(String(name_x)) versus $(String(name_y))}")
        println(styled"{red: Ignoring correction for samples outside range}")
        correction = 1
    else
        println(styled"{green: range used is appropriate to determine HDI of $(String(name_x)) versus $(String(name_y))}")
        correction = 1/(1 - frac_lost)
    end

    if nbins_contour != -1
        edges = (LinRange(range_x[1], range_x[2], nbins_contour+1),
                LinRange(range_y[1], range_y[2], nbins_contour+1))
        h_ct = fit(Histogram, (values_x[filter], values_y[filter]), weights(sample_weights[filter]), edges) 
        x_ct = (h_ct.edges[1][2:end] .+ h_ct.edges[1][1:end-1])./2
        y_ct = (h_ct.edges[2][2:end] .+ h_ct.edges[2][1:end-1])./2
        dx_ct = x_ct[2]-x_ct[1]
        dy_ct = x_ct[2]-x_ct[1]
        z_ct = h_ct.weights/(total_weight*dx_ct*dy_ct)
    else
        # use the same values from the heatmap
        h_ct = h_hm
        x_ct = x_hm
        y_ct = y_hm
        dx_ct = dx_hm
        dy_ct = dy_hm
        z_ct = z_hm
    end
    bounds = get_bounds_for_fractions(h_ct, fractions./correction)
    if !isnothing(axis)
        contour!(axis, x_ct, y_ct, h_ct.weights, levels=bounds, color=(:black, 0.5))
    end
    
    return x_hm, y_hm, z_hm, dx_hm, dy_hm, x_ct, y_ct, z_ct, dx_ct, dy_ct, bounds, frac_lost
end  

"""
    plot_1D_density(axis, values, range, chain_weights, fraction_1D, nbins; color, linewidth)

Creates a 1D marginalized distribution plot in `axis` from the sample values given in `values`.
The plotted line is normalized, such that it corresponds to the PDF followed by the samples.

# Arguments:
- axis: The Axis on which the distribution will be plotted. If set to `nothing`, no plot will be made,
but the values used for the plot will still be returned.
- values: Either a vector or a matrix of numbers, corresponding to a single chain or multiple chains
respectively (in the latter case, columns represent individual samples).
- range: Vector with two values containing the range of the plot. It also corresponds to the range
that will be used to determine the HDI. The HDI does correct for values outside the range, assuming
that these do not fall within the highest density region.
- chain_weights: weights used for each sample, should have the same size as `values`.
- fraction_1d: Fraction of samples (including weights) condained within the HDI for which the
credible interval is reported.
- nbins: Number of bins used to build the histogram and determine the HDI
- color: color used for the line plot.
- linewidth: linewidth used for the line plot

# Output:
- x: Bin centers used in the histogram
- h: Output of `Histogram` fit. This is not normalized, so it contains the sum of all weights
(or count of all samples) within each bin.
- y: Value of the PDF estimated at each bin
- frac_lost: Fraction of samples (including weights) that is outside of `range`
"""
function plot_1D_density(axis, values, range, chain_weights, nbins; color, linewidth)

    filter = values .> range[1] .&& values .< range[2]
    values = values[filter]
    total_weight = sum(chain_weights)
    missing_weight = sum(chain_weights[.!filter])
    chain_weights = weights(chain_weights[filter]) # weights is a StatsBase function

    frac_lost = missing_weight/total_weight
    
    edges = LinRange(range[1], range[2], nbins+1)
    h = fit(Histogram, values, chain_weights, edges)
    x =(h.edges[1][2:end] .+ h.edges[1][1:end-1])./2
    dx = x[2]-x[1]
    y = h.weights/(total_weight*dx)
    if !isnothing(axis)
        lines!(axis, x, y, color=color, linewidth=linewidth)
    end
    return x, h, y, dx, frac_lost
end

"""
    plot_compund_1D_density(axis, values, range, chain_weights, fraction_1D, nbins; color, linewidth)

Creates a 1D marginalized distribution plot in `axis` from the sample values given in `values`.
If values contains more than one chain, it will plot the individual chains separetely, as well
as one line plot containing all of the samples from all chains.

# Arguments:
- axis: The Axis on which the distribution will be plotted. If not provided, no plot will be made,
but the values used for the plot will still be returned.
- name: Symbol containing the name of the variable.
- values: Either a vector or a matrix of numbers, corresponding to a single chain or multiple chains
respectively (in the latter case, columns represent individual samples).
- range: Vector with two values containing the range of the plot. It also corresponds to the range
that will be used to determine the HDI. The HDI does correct for values outside the range, assuming
that these do not fall within the highest density region.
- sample_weights: weights used for each sample, should have the same size as `values`.
- fraction_1d: Fraction of samples (including weights) condained within the HDI for which the
credible interval is reported.
- nbins: Number of bins used to build the histogram and determine the HDI

# Output:
- (xmin, xmode, xmax): Credible interval edges and mode.
- x: Bin centers used in the histogram. This will include all samples from all chains.
- h: Output of `Histogram` fit. This is not normalized, so it contains the sum of all weights
(or count of all samples) within each bin.
- y: Value of the PDF estimated at each bin
- frac_lost: Fraction of samples (including weights) that is outside of `range`
"""
function plot_compound_1D_density(axis, name, values_x, range_x, sample_weights, fraction_1D, nbins)

    # If we have multiple chains, iterate over them
    if ndims(values_x)==2 && !isnothing(axis)
        for ii in 1:size(values_x)[2] # nchains
            chain_values = @view values_x[:,ii]
            chain_weights = @view sample_weights[:,ii]
            plot_1D_density(axis, chain_values, range_x, chain_weights, nbins, color=(:gray, 0.25), linewidth=1)
        end
    end

    # Plot once for all the values
    x, h, y, dx, frac_lost = plot_1D_density(axis, values_x, range_x, sample_weights, nbins, color=(:blue, 1.0), linewidth=1)

    #we correct for samples outside the range, assuming they are
    #not part of the HDI. A basic check is to verify we dont have too
    #many samples outside.
    if 1 - fraction_1D < frac_lost
        println(styled"{red: range used is insufficient to determine HDI of $(String(name))}")
        println(styled"{red: Ignoring correction for samples outside range}")
        correction = 1
    else
        println(styled"{green: range used is appropriate to determine HDI of $(String(name))}")
        correction = 1/(1 - frac_lost)
    end

    bound = get_bounds_for_fractions(h, [fraction_1D*correction])[1]
    xmin = minimum(x[h.weights .>= bound]) - dx/2 # get left most value of bin
    xmax = maximum(x[h.weights .>= bound]) + dx/2
    xmode = x[argmax(h.weights)]
    if xmode == x[1]
        xmode = range_x[1]
    elseif xmode == x[end]
        xmode = range_x[2]
    end

    if !isnothing(axis)
        filter = x .>= xmin .&& x .<= xmax
        band!(axis, x[filter], zeros(length(x[filter])), y[filter], color=(:gray, 0.4))
        vlines!(axis, xmode, color=(:black, 1.0), linewidth=1)
        xlims!(axis, range_x[1], range_x[2])
        ylims!(axis, 0, 1.1*maximum(y))
    end

    return (xmin, xmode, xmax), x, h, y, dx, frac_lost
   
end

"""
    plot_extra_1D_distribution()

Plots a distribution in one of the panels containing a marginalized 1D distribution
of a corner plot The Distribution will be plotted within the range determined for the panel.
# Arguments:
- corner_plot: An instance of CornerPlot containing a plotted figure.
- name_x: Identifies the panel in which the distribution will be plotted.
- distribution: 1D distribution to be plotted.
- npoints: Number of points of the line plot.
- linewidth: Property for line plot.
- color: Property for line plot.
- linestyle: Property for line plot.
"""
function plot_extra_1D_distribution(corner_plot, name_x, distribution::Distribution; npoints=100, linewidth=2, color=:red, linestyle=:solid)
    xvals = LinRange(corner_plot.ranges[name_x][1], corner_plot.ranges[name_x][2],npoints)
    yvals = pdf(distribution, xvals)
    lines!(corner_plot.distributions_1d[name_x], xvals, yvals,
            linewidth=linewidth, color=color, linestyle=linestyle)
end
