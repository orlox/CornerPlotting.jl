"""
    MultiCornerPlot(results, names; ...)

Constructs a corner plot from the provided "results", using the variables given
in the "names" vector. Results here is taken to contain multiple datasets that will be
overplotted together.

# Arguments:
- results: Vector containing the samples to be plotted. This can be either of type Vector{MCMCChains}
or a vector of Dictionaries containing vectors with values for an individual chain, or arrays for multiple
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
- oneD_lines_kwargs : Named tuple containing the keyword arguments used in the call to `lines!` for individual chains
- oneD_lines_full_kwargs: Named tuple containing keyword arguments used in the call to `lines!` for all chains grouped
- oneD_band_default_kwargs: Named tuple containing keyword arguments used in the call to `band!`
- oneD_vlines_default_kwargs: Named tuple containing keyword arguments used in the call to `vlines!`
- twoD_heatmap_kwargs: named tuple containing named arguments for the call to `heatmap!`
- twoD_contour_kwargs: named tuple containing named arguments for the call to `contour!`

# Output:
Returns an instance of CornerPlot
"""
function MultiCornerPlot(results, names::Vector{Symbol};
        labels=nothing, ranges=Dict(), scaling=Dict(),
        fig=Figure(), quantile_for_range=0.01,
        use_weights = true, fraction_1D=0.9, fractions_2D=[0.9], 
        show_CIs=false, nbins=100, nbins_contour=20,
        axis_size=100,
        )
    corner_plot = nothing

    # verify data exists
    for name in names
        for result in results
            if name ∉ keys(result)
                throw(ArgumentError("$name is not a valid key"))
            end
        end
    end

    # determine ranges for all together
    for name in names
        if name ∉ keys(ranges)
            for (i,result) in enumerate(results)
                if :weights ∈ keys(results) && use_weights
                    sample_weights = results[:weights]
                    subrange = quantile(vec(result[name]),weights(sample_weights),
                                            [quantile_for_range, 1-quantile_for_range])
                else
                    subrange = quantile(vec(result[name]),
                                            [quantile_for_range, 1-quantile_for_range])
                end
                if name ∈ keys(scaling)
                    subrange /= scaling[name]
                end
                if i==1
                    ranges[name] = subrange
                else
                    ranges[name][1] = min(ranges[name][1], subrange[1])
                    ranges[name][2] = max(ranges[name][2], subrange[2])
                end
            end
        end
    end

    for result in results
        @show "caca",ranges
        corner_plot = CornerPlot(result, names; labels=labels, ranges=copy(ranges), scaling=scaling,
        fig=fig, quantile_for_range=quantile_for_range,
        use_weights = use_weights, fraction_1D=fraction_1D, fractions_2D=fractions_2D, 
        show_CIs=show_CIs, show_heatmap=false, nbins=nbins, nbins_contour=nbins_contour,
        axis_size=axis_size, corner_plot=corner_plot)
    end
    return corner_plot
end