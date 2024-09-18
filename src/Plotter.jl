using StatsBase
using CairoMakie
using StyledStrings
using Distributions



"""
    struct CornerPlot

Description

# Arguments:

# Output:
"""
struct CornerPlot
    fig
    ranges
    distributions_1d
    distributions_2d
end

"""
    CornerPlot()

Description

# Arguments:

# Output:
"""
function CornerPlot(results, names::Vector{Symbol};
        labels=nothing, ranges=Dict(), scaling=nothing,
        fig=Figure(), quantile_for_range=0.01,
        use_weights = true, fraction_1D=0.9, fractions_2D=[0.9], 
        show_CIs=true, nbins=100, nbins_contour=20,
        axis_size=100)

    num_col = length(names)
    
    distributions_2d = Dict() # this will store the 2D distributions in a dictionary of dictionaries
    for name in names
        if name ∉ keys(results)
            throw(ArgumentError("$name is not a valid key"))
        end
        if name ∉ keys(ranges)
            ranges[name] = quantile(results[name],[quantile_for_range, 1-quantile_for_range])
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
        weights = results[:weights]
    else
        # just create an array of ones with the needed size
        # not the most efficient but avoids having to write special cases later
        weights = ones(size(results[names[1]]))
    end
 
    # Create 2D density plots
    for ii in 1:num_col-1         # ii is the x-coord param
        name_x = names[ii]
        for  jj in ii+1:num_col   # jj is the y-coord param
            name_y = names[jj]
            axis = Axis(fig[jj+1,ii], xlabel=labels[name_x], ylabel=labels[name_y], height=axis_size, width=axis_size)
            distributions_2d[name_x][name_y] = axis
            if isnothing(scaling)
                values_x = results[name_x]
                values_y = results[name_y]
            else
                values_x = results[name_x]./scaling[ii]
                values_y = results[name_y]./scaling[jj]
            end

            plot_2D_density(axis, name_x, name_y, values_x, ranges[name_x], values_y, ranges[name_y],
                weights, fractions_2D, nbins_heatmap=nbins, nbins_contour=nbins_contour)
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
        if isnothing(scaling)
            values_x = results[name_x]
        else
            values_x = results[name_x]./scaling[ii]
        end
        (xmin, xmode, xmax) = plot_compound_1D_density(axis, name_x, values_x, ranges[name_x], weights, fraction_1D, nbins)

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
    resize_to_layout!(fig)
    
    return CornerPlot(fig, ranges, distributions_1d, distributions_2d)
end

"""
    get_bounds_for_fractions(h, fractions)

Description
Calculate the bounds containing the specified fraction(s) of area.

# Arguments:
- h:         the densities contained in the bins
- fractions: the fractional area that should be bounded

# Output:
- bounds:    the limits of the bounding area
"""
function get_bounds_for_fractions(h, fractions)
    integral = sum(h.weights)
    bounds =zeros(length(fractions))
    for (jj,fraction) in enumerate(fractions)
        minbound = 0
        maxbound = maximum(h.weights)
        newbound = 0
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
    create_2D_density()

Description

# Arguments:
"""
function plot_2D_density(axis, name_x, name_y, values_x, ranges_x, values_y, ranges_y, sample_weights, fractions; nbins_heatmap, nbins_contour)

    filter = (values_x .> ranges_x[1]) .&& (values_x .< ranges_x[2]) .&& (values_y .> ranges_y[1]) .&& (values_y .> ranges_y[1])
    total_weight = sum(sample_weights)
    missing_weight = sum(sample_weights[.!filter])
    frac_lost = missing_weight/total_weight

    h_hm = fit(Histogram, (values_x[filter], values_y[filter]), weights(sample_weights[filter]), nbins=nbins_heatmap) 
    x_hm = (h_hm.edges[1][2:end] .+ h_hm.edges[1][1:end-1])./2
    y_hm = (h_hm.edges[2][2:end] .+ h_hm.edges[2][1:end-1])./2
    heatmap!(axis, x_hm, y_hm, h_hm.weights, colormap=:dense)

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

    h_ct = fit(Histogram, (values_x[filter], values_y[filter]), weights(sample_weights[filter]), nbins=nbins_contour) 
    x_ct = (h_ct.edges[1][2:end] .+ h_ct.edges[1][1:end-1])./2
    y_ct = (h_ct.edges[2][2:end] .+ h_ct.edges[2][1:end-1])./2
    bounds = get_bounds_for_fractions(h_ct, fractions./correction)
    contour!(axis, x_ct, y_ct, h_ct.weights, levels=bounds, color=(:black, 0.5))


end  

"""
    plot_1D_density()

Description

# Arguments:

# Output:
"""
function plot_1D_density(axis, values, range, chain_weights, fraction_1D, nbins; color, linewidth)

    filter = values .> range[1] .&& values .< range[2]
    values = values[filter]
    total_weight = sum(chain_weights)
    missing_weight = sum(chain_weights[.!filter])
    chain_weights = weights(chain_weights[filter]) # weights is a StatsBase function

    frac_lost = missing_weight/total_weight
    
    h = fit(Histogram, values, chain_weights, nbins=nbins)
    x =(h.edges[1][2:end] .+ h.edges[1][1:end-1])./2
    dx = 1
    if length(x) > 1
        if x[2]-x[1] > 0
            dx = x[2]-x[1]
        end
    end
    y = h.weights/(total_weight*dx)
    lines!(axis, x, y, color=color, linewidth=linewidth)
    return x, h, y, dx, frac_lost
end

"""
    create_marginalized_1D_densities()

# Arguments:

# Output:
"""
function plot_compound_1D_density(axis, name, values_x, range_x, sample_weights, fraction_1D, nbins)

    # If we have multiple chains, iterate over them
    if ndims(values_x)==2
        for ii in 1:size(values_x)[2] # nchains
            chain_values = @view values_x[:,ii]
            chain_weights = @view sample_weights[:,ii]
            plot_1D_density(axis, chain_values, range_x, chain_weights, fraction_1D, nbins, color=(:gray, 0.25), linewidth=1)
        end
    end

    # Plot once for all the values
    x, h, y, dx, frac_lost = plot_1D_density(axis, values_x, range_x, sample_weights, fraction_1D, nbins, color=(:blue, 1.0), linewidth=1)

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

    filter = x .>= xmin .&& x .<= xmax
    band!(axis, x[filter], zeros(length(x[filter])), y[filter], color=(:gray, 0.4))
    vlines!(axis, xmode, color=(:black, 1.0), linewidth=1)
    xlims!(axis, range_x[1], range_x[2])
    ylims!(axis, 0, 1.1*maximum(y))

    return (xmin, xmode, xmax)
   
end

"""
    plot_extra_1D_distribution()

Description

# Arguments:

# Output:
"""
function plot_extra_1D_distribution(corner_plot, name_x, distribution::Distribution; npoints=100, linewidth=2, color=:red, linestyle=:solid)
    xvals = LinRange(corner_plot.ranges[name_x][1], corner_plot.ranges[name_x][2],npoints)
    yvals = pdf(distribution, xvals)
    lines!(corner_plot.distributions_1d[name_x], xvals, yvals,
            linewidth=linewidth, color=color, linestyle=linestyle)
end