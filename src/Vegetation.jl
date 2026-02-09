module Vegetation

using EarthSciMLBase
using ModelingToolkit
using ModelingToolkit: t, D
using DynamicQuantities
using DocStringExtensions

include("landis_biomass.jl")
export LANDISBiomass

end
