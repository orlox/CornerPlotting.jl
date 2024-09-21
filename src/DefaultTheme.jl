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