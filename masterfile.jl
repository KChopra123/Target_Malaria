include("drivematrix.jl")
include("femalefitnessmatrix.jl")
include("linindices.jl")
include("mutationmatrix.jl")
include("newfrequency.jl")
include("TargetSiteResistance.jl")
include("resistancehaplotypes.jl")
include("resistanceprob.jl")
include("TargetSiteResistanceRecovery.jl")
include("recoveryprob.jl")
include("Harry_Anthony_Code\\GaussPoissonHybrid_mnrnd.jl")
include("Harry_Anthony_Code\\Allele_var_Mtx.jl")
include("migration\\fixedmig_TSR.jl")
include("migration\\sqmigration_TSR.jl")
include("migration\\square_deme_migration.jl")
include("migration\\square_deme_size.jl")
include("migration\\sqmignewfrequency.jl")

using Distributions
using LinearAlgebra
using Plots
using PoissonRandom
using StatsBase
using SpecialFunctions
