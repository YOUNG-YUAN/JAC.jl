
"""
`module  JAC.AtomicFeatures`  
... a submodel of JAC that contains all methods for computing the features for a neural-network training. 
    It assumes that the electronic and level structure of atoms and ions can be encoded by (large, real)
    x-vectors (feature vectors) that are based on a proper list of shells, subshells, orbitals as well as
    the configuration-interaction and coupling information of leading CSF.
"""
module AtomicFeatures


using Printf, ..AngularMomentum, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..TableStrings

export  AtomicModel


"""
`struct  AtomicFeatures.AtomicModel`  
    ... defines a type for specifying the atomic model that is used in the training, test and use of an NN. It contains 
        all necessary data but not the underlying "recipi/scheme", how these data are generated and utilized in dealing 
        with the NN. At present, it is assumed that each NN is based on a set of shells (subshells) and associated
        orbitals, from which are features (descriptors) can be generated/extracted.
        

    + nMax              ::Int64                      ... maximum n-shell in the feature generation/extraction.
    + lMax              ::Int64                      ... maximum l-shell in the feature generation/extraction.
    + grid              ::Radial.Grid                ... The radial grid to be used for the computation.
    + asfSettings       ::ManyElectron.AsfSettings   ... AsfSettings used for simplified atomic-structure computations.
    + nuclearModel      ::Nuclear.Model              ... nuclear model used to define the orbitals.
    + shells            ::Array{Shell,1}             ... List of shells of model considered.
    + subshells         ::Array{Subshell,1}          ... List of subshells; this is redundant but simplifies the feature extraction. 
    + orbitals01        ::Dict{Subshell, Orbital}    ... Set of associated orbitals for charge state q=1.
    + orbitals02        ::Dict{Subshell, Orbital}    ... Set of associated orbitals for charge state q=2.
    + orbitals03        ::Dict{Subshell, Orbital}    ... Set of associated orbitals for charge state q=3.
"""
struct  AtomicModel
    nMax                ::Int64
    lMax                ::Int64  
    grid                ::Radial.Grid 
    asfSettings         ::ManyElectron.AsfSettings 
    nuclearModel        ::Nuclear.Model 
    shells              ::Array{Shell,1}  
    subshells           ::Array{Subshell,1} 
    orbitals01          ::Dict{Subshell, Orbital}
    orbitals02          ::Dict{Subshell, Orbital}
    orbitals03          ::Dict{Subshell, Orbital}
end 


# `Base.show(io::IO, am::AtomicFeatures.AtomicModel)`  ... prepares a proper printout of the variable am::AtomicFeatures.AtomicModel
function Base.show(io::IO, am::AtomicFeatures.AtomicModel) 
    println(io, "nMax:              $(am.nMax)  ")
    println(io, "lMax:              $(am.lMax)  ")
    println(io, "grid:              $(am.grid)  ")
    println(io, "asfSettings:       $(am.asfSettings)  ")
    println(io, "nuclearModel:      $(am.nuclearModel)  ")
    println(io, "shells:            $(am.shells)  ")
    println(io, "subshells:         $(am.subshells)  ")
    println(io, "orbitals01:        $(am.orbitals01)  ")
    println(io, "orbitals02:        $(am.orbitals02)  ")
    println(io, "orbitals03:        $(am.orbitals03)  ")
end


"""
`struct  AtomicFeatures.LevelXVector`  
    ... defines a type for specifying a particular LSJ-coupled atomic level in terms of its LSJ-quantum
        numbers, a computed energy as well as a generate feature x-vector.        

    + L                 ::AngularJ64          ... total L of the LSJ-level
    + S                 ::AngularJ64          ... total S of the LSJ-level
    + J                 ::AngularJ64          ... total J of the LSJ-level
    + parity            ::Basics.Parity       ... total parity J of the LSJ-level
    + energy            ::Float64             ... total energy, compute with the given atomic model.
    + nistEnergy        ::Float64             ... excitation energy, if available in NIST, and zero otherwise.
    + xVector           ::Array{Float64,1}    ... x-vector, generated with the given atomic model. 
"""
struct  XyVector
    L                   ::AngularJ64 
    S                   ::AngularJ64
    J                   ::AngularJ64 
    parity              ::Basics.Parity 
    energy              ::Float64  
    nistEnergy          ::Float64 
    xVector             ::Array{Float64,1} 
end 


# `Base.show(io::IO, v::AtomicFeatures.XyVector)`  ... prepares a proper printout of the variable v::AtomicFeatures.XyVector
function Base.show(io::IO, v::AtomicFeatures.XyVector) 
    println(io, "L:           $(v.L)  ")
    println(io, "S:           $(v.S)  ")
    println(io, "J:           $(v.J)  ")
    println(io, "parity:      $(v.parity)  ")
    println(io, "energy:      $(v.energy)  ")
    println(io, "nistEnergy:  $(v.nistEnergy)  ")
    println(io, "xVector:     $(v.xVector)  ")
end

end # module

