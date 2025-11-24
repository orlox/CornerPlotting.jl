using CairoMakie
using MathTeXEngine


"""
    default_theme()

Returns a Makie theme that will produce a nice plot.
"""
function default_theme(;rowcolgap=10, xticklabelrotation=pi/4,labelfontsize=16, credible_interval_fontsize=16, tickfontsize=10)
    return Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                           italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                    fontsize=credible_interval_fontsize, linewidth=2,rowgap=rowcolgap, colgap=rowcolgap,
                 Axis=(xlabelsize=labelfontsize, ylabelsize=labelfontsize, xgridvisible=false, ygridvisible=false,
                       spinewidth=1, xminorticksvisible=false, yminorticksvisible=false, xtickalign=1, ytickalign=1,
                       xminortickalign=1, yminortickalign=1, xticksize=5, xtickwidth=1, yticksize=5,
                       ytickwidth=1, xminorticksize=7, xminortickwidth=1, yminorticksize=5, yminortickwidth=1,
                       xticklabelsize=tickfontsize, yticklabelsize=tickfontsize, xticksmirrored=true, yticksmirrored=true,
                       xticklabelrotation=xticklabelrotation),
                  )
end

"""
    oneD_lines_default_kwargs
    oneD_lines_multi_default_kwargs
Default style options when doing a line in the 1D marginalized plots
for an individual chain. `multi` refers to the `MultiCornerPlot` functionality,
and defines a vector that will be looped over (returning to the beginning)
"""
oneD_lines_default_kwargs = (color=(:gray, 0.25), linewidth=1)
oneD_lines_multi_default_kwargs = [(color=(:gray, 0.25), linewidth=1)]

"""
    oneD_lines_full_default_kwargs
    oneD_lines_full_multi_default_kwargs
Default style options when doing a line in the 1D marginalized plots
for all chains together (or for the single chain if only one is available)
"""
oneD_lines_full_default_kwargs = (color=(:blue, 1.0), linewidth=1)
oneD_lines_full_multi_default_kwargs = [
        (color=(:blue, 1.0), linewidth=1),
        (color=(:orange, 1.0), linewidth=1),
        (color=(:red, 1.0), linewidth=1),
        (color=(:gold, 1.0), linewidth=1),
        (color=(:teal, 1.0), linewidth=1)
    ]

"""
    oneD_band_default_kwargs
    oneD_band_multi_default_kwargs
Default style options when doing a band in the 1D marginalized plots
to show the credible interval
"""
oneD_band_default_kwargs = (color=(:gray, 0.4),)
oneD_band_multi_default_kwargs = [
        (color=(:blue, 0.4),),
        (color=(:orange, 0.4),),
        (color=(:red, 0.4),),
        (color=(:gold, 0.4),),
        (color=(:teal, 0.4),)
    ]

"""
    oneD_vlines_default_kwargs
    oneD_vlines_multi_default_kwargs
Default style options when doing a vline in the 1D marginalized plots
to show the mode
"""
oneD_vlines_default_kwargs = (color=(:black, 1.0), linewidth=1)
oneD_vlines_multi_default_kwargs = [
        (color=(:blue, 1.0), linewidth=1),
        (color=(:orange, 1.0), linewidth=1),
        (color=(:red, 1.0), linewidth=1),
        (color=(:gold, 1.0), linewidth=1),
        (color=(:teal, 1.0), linewidth=1),
    ]

"""
    twoD_heatmap_default_kwargs
Default style options for the call to `heatmap!` in the 2D marginalized plots
"""
twoD_heatmap_default_kwargs = (colormap=:dense,)

"""
    twoD_contour_default_kwargs
    twoD_contour_multi_default_kwargs
Default style for the call to `contour!` in the 2D marginalized plots
"""
twoD_contour_default_kwargs = (color=(:black,0.5),)
twoD_contour_multi_default_kwargs = [
        (color=:blue,),
        (color=:orange,),
        (color=:red,),
        (color=:gold,),
        (color=:teal,),
    ]