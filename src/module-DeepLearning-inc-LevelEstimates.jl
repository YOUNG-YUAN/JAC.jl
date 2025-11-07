
# Date structures, functions and methods for (level) estimates of levels that are missing in the NIST database 


"""
`struct  DeepLearning.LsjLevel`  
    ... defines a type for collecting information about a single LSJ with two different energies, a calculated or
        predicted energy and the associated NIST, if this assignment is possible.

    + L                 ::AngularJ64          ... total L of the LSJ-level
    + S                 ::AngularJ64          ... total S of the LSJ-level
    + J                 ::AngularJ64          ... total J of the LSJ-level
    + parity            ::Basics.Parity       ... total parity J of the LSJ-level
    + thisEnergy        ::Float64             ... excitation energy, either compute or predicted.
    + nistEnergy        ::Float64             ... excitation energy, either compute or predicted.
"""
struct  LsjLevel
    conf                ::Configuration
    L                   ::AngularJ64 
    S                   ::AngularJ64
    J                   ::AngularJ64  
    parity              ::Basics.Parity
    thisEnergy          ::Float64 
    nistEnergy          ::Float64  
end 


# `Base.show(io::IO, level::DeepLearning.LsjLevel)`  ... prepares a proper printout of the variable level::DeepLearning.LsjLevel
function Base.show(io::IO, level::DeepLearning.LsjLevel) 
    println(io, "L:           $(level.L)  ")
    println(io, "S:           $(level.S)  ")
    println(io, "J:           $(level.J)  ")
    println(io, "parity:      $(level.parity)  ")
    println(io, "thisEnergy:  $(level.thisEnergy)  ")
    println(io, "nistEnergy:  $(level.nistEnergy)  ")
end


"""
`struct  DeepLearning.NistLevel`  
    ... defines a type for collecting information about a single level from the NIST Level database.

    + conf              ::Configuration       ... Configuration of the LSJ-level
    + L                 ::AngularJ64          ... total L of the LSJ-level
    + S                 ::AngularJ64          ... total S of the LSJ-level
    + J                 ::AngularJ64          ... total J of the LSJ-level
    + parity            ::Basics.Parity       ... total parity of the LSJ-level
    + energy            ::Float64             ... total energy, compute with the given atomic model.
"""
struct  NistLevel
    conf                ::Configuration
    L                   ::AngularJ64 
    S                   ::AngularJ64
    J                   ::AngularJ64 
    parity              ::Basics.Parity
    energy              ::Float64  
end 


# `Base.show(io::IO, level::DeepLearning.NistLevel)`  ... prepares a proper printout of the variable level::DeepLearning.NistLevel
function Base.show(io::IO, level::DeepLearning.NistLevel) 
    println(io, "conf:        $(level.conf)  ")
    println(io, "L:           $(level.L)  ")
    println(io, "S:           $(level.S)  ")
    println(io, "J:           $(level.J)  ")
    println(io, "parity:      $(level.parity)  ")
    println(io, "energy:      $(level.energy)  ")
end


"""
`struct  DeepLearning.LevelEstimationRequest  <:  DeepLearning.AbstractNeuralNetworkRequest`  
    ... to define a (deep-learning) request for estimating the level energies for a given set of configurations.

    + configs       ::Array{Configuration,1}      ... List of configuration for which levels are to be estimated.             
    + nistLevels    ::Array{NistLevel,1}          ... Available list of NIST levels used for comparison.            
"""
struct   LevelEstimationRequest  <:  DeepLearning.AbstractNeuralNetworkRequest
    configs         ::Array{Configuration,1}            
    nistLevels      ::Array{NistLevel,1}              
end


"""
`DeepLearning.LevelEstimationRequest()`  ... constructor for an 'default' instance of a DeepLearning.LevelEstimationRequest.
"""
function LevelEstimationRequest()
    LevelEstimationRequest( Configuration[], NistLevel[])
end


# `Base.string(request::LevelEstimationRequest)`  ... provides a String notation for the variable request::LevelEstimationRequest.
function Base.string(request::LevelEstimationRequest)
    sa = "NN request for level estimation:"
    return( sa )
end


# `Base.show(io::IO, request::LevelEstimationRequest)`  ... prepares a proper printout of the request::LevelEstimationRequest.
function Base.show(io::IO, request::LevelEstimationRequest)
    sa = Base.string(request);        print(io, sa, "\n")
    println(io, "configs:             $(request.configs)  ")
    println(io, "nistLevels[1:3]:     $(request.nistLevels[1:3])  ")
end



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


"""
`DeepLearning.checkConfiguration(conf::Configuration, am::AtomicFeatures.AtomicModel)` 
    ... checks that the given configuration can be encoded into/used by the given atomic model;
        a boolian value of true/false is returned. If false is returned, a warning is printed about
        the reasons, usually concerning the maximum n quantum numbers involved.
"""
function checkConfiguration(conf::Configuration, am::AtomicFeatures.AtomicModel)
    shells = Basics.extractShellList([conf]);   val = true

    # No shell should be larger than nMax and lMax
    for  shell in shells
        if  shell.n > am.nMax  ||  shell.n > am.nMax   val = false
            @warn "$conf has shells > $(Shell(am.nMax, am.lMax)) "
            break
        end 
    end

    return( val )
end


"""
`DeepLearning.computeMultiplet(conf::Configuration, am::AtomicFeatures.AtomicModel)` 
    ... computes the level multiplet as associated with the given configuration, (frozen) set of orbitals and asfSettings;
        these settings just control the CI computations for the levels of a single configuration.
"""
function computeMultiplet(conf::Configuration, am::AtomicFeatures.AtomicModel)
    orbitals  = DeepLearning.extractOrbitals(conf::Configuration, am::AtomicFeatures.AtomicModel)
    multiplet = Hamiltonian.performCIwithFrozenOrbitals([conf], orbitals, am.nuclearModel, am.grid, am.asfSettings; 
                                                        printout=true)
    return( multiplet )
end


"""
`DeepLearning.extractLsjLevels(conf::Configuration, multiplet::Multiplet, nistLevels::Array{NistLevel,1})` 
    ... extracts all LsjLevels that are associated with the given configuration. It also assign the correct
        excitation energy from the NIST tables, if available. This procedure is indeed very crucial and intricate
        as it is diffifult to make this assigment accurate. Further work is needed, perhaps, to consider neighbored
        levels. A list of levels::Array{NistLevel,1} is returned, in which the given NIST energies are set
        and 0. otherwise.
        
        The procedure need to perform several steps: 
        (i)    The LSJ decomposition of all levels in the multiplet;
        (ii)   The analysis of the NIST levels which belong to the given configuration;
        (iii)  The assignment of the right energies to each level, as far as available.
"""
function extractLsjLevels(conf::Configuration, multiplet::Multiplet, nistLevels::Array{NistLevel,1})
    lsjLevels = LsjLevel[];   relevantLevels = NistLevel[]
    
    # Perform a jj-LL transformation of all levels in the multiplet
    
    # Extract the correct levels for the given configuration
    for level in nistLevels
        if  conf == level.conf   push!(relevantLevels, level)   end
    end
    
    # Assign the computed LSJ levels and the relevant NIST levels to each other 
    error("No assignment done yet.")
    
    return( lsjLevels )
end


"""
`DeepLearning.displayLsjLevels(conf::Configuration, lsjLevels::Array{LsjLevel,1})` 
    ... displays (and lists) all given LsjLevels in a neat format; nothing is returned.
"""
function displayLsjLevels(conf::Configuration, lsjLevels::Array{LsjLevel,1})
    nx = 105
    println(" ")
    println("  Comparison of calculated/predicted (excitation) energies with NIST energies for $conf:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(14, "level  "       ; na=0);          sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(18, "^(2S+1) L_J^P" ; na=2);          sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(10, "This energy"   ; na=3);              
    sb = sb * TableStrings.center(10, "eV"            ; na=3)
    sa = sa * TableStrings.center(10, "NIST energy"   ; na=3)              
    sb = sb * TableStrings.center(10, "eV"            ; na=2)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    # 
    for  (i, level) in enumerate(lsjLevels)
        sym = LevelSymmetry( level.J, level.parity)
        sm  = round(Int64, Basics.twice(level.S)+1)
        sL  = round(Int64, Basics.twice(level.L)/2.)
        si  = "    " * string(i)
        sy  = string(sym) * "         "
        sa  = "   " * si[end-4:end] * "   " * "^" * string(sm) * " " * string(sL) * "_" * sy[1:8] 
        sa = sa * @sprintf("%.4e", level.thisEnergy)           * "   "
        sa = sa * @sprintf("%.4e", level.nistEnergy)           * "   "
        println( sa )  
    end
    println("  ", TableStrings.hLine(nx), "\n")

    return( nothing )
end


"""
`DeepLearning.extractNistConfigurations(nistLevels::Array{NistLevel,1})` 
    ... extracts all configurations that are defined by the given list of NIST levels.
"""
function extractNistConfigurations(nistLevels::Array{NistLevel,1})
    configs = Configuration[]
    
    for level  in  nistLevels   push!(configs, level.conf)  end 
    configs = unique(configs)

    return( configs )
end


"""
`DeepLearning.extractNistLevels(filesnames::Array{String,1})` 
    ... extracts all available NIST levels from a (number of) ASCII-files, in which these data are encoded.
"""
function extractNistLevels(filesnames::Array{String,1})
    nistLevels = NistLevel[]
    
    for filename  in  filesnames
        raw_lines = readlines(filename)
        for line in raw_lines
            # Extract the configuration, quantum numbers and energies from the line, if appropriate
            push!(nistLevels, NistLevel(conf, L, S, J, parity, energy))
        end
    end

    return( nistLevels )
end


"""
`DeepLearning.extractOrbitals(conf::Configuration, am::AtomicFeatures.AtomicModel)` 
    ... extracts for the given configuration the orbitals from the atomic model.
"""
function extractOrbitals(conf::Configuration, am::AtomicFeatures.AtomicModel)
    q = round(Int64, am.nuclearModel.Z - conf.NoElectrons)

    if      q == 1   orbitals = am.orbitals01
    elseif  q == 2   orbitals = am.orbitals02
    elseif  q == 3   orbitals = am.orbitals03
    else    error("stop a")
    end

    return( orbitals )
end
    

"""
`DeepLearning.generateAtomicModelForLE_Arn4()` 
    ... generates an atomic model for the training, test and feature extraction for the (level) estimation of atomic levels. 
        Here, in particular, we consider the levels of Ar^+, Ar^2+ and Ar^3+ and configurations with maximum nMax = 4 shells.
        All details are hard-coded and follow some 'script-like' style. We expect that a particular generateAtomicModel...() 
        is designed for each neutral network, which we shall train and consider. 
        
        An atomicModel::AtomicFeatures.AtomicModel is returned.
"""
function generateAtomicModelForLE_Arn4()
    nMax = 4;   lMax = 3
    nm          = Nuclear.Model(18.)
    grid        = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0) 
    asfSettings = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
    
    # Generate all shell and subshells
    shells    = Basics.generateShellList(1, nMax, lMax)
    subshells = Basics.generateSubshellList(shells)
    
    # Generate all orbital sets in turn 
    mfSettings  = MeanFieldSettings(AsfSettings.scField)
    exciteE     = Basics.Excitelectrons(1)
    
    configs01   = [ Configuration("[Ne] 3s^2 3p^5")]
    configs01   = Basics.generateConfigurations(exciteE, configs01)
    meanField   = Representation("Ar^+", nm, grid, configs01, MeanFieldBasis(mfSettings) )
    mfrep       = generate(meanField; output=true)
    orbitals01  = mfrep["orbitals"]

    # Generate better orbitals
    orbitals02  = mfrep["orbitals"]
    orbitals03  = mfrep["orbitals"]
    
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    
    atomicModel = AtomicFeatures.AtomicModel(nMax, lMax, grid, asfSettings, shells, subshells, 
                                             orbitals01, orbitals02, orbitals03)

    return( atomicModel )
end


"""
`DeepLearning.generateFeatureVectors(confs::Array{Configuration,1}, am::AtomicFeatures.AtomicModel, probability::Float64,
                                     nistLevels::Array{NistLevel,1})` 
    ... generates a list of -- training and test -- feature vectors which are associated with the given configurations.
        These feature vectors can be written out and used for the training and testing of NN as needed by the 
        applied DeepLearning toolboxes (outside of JAC). The procedure assumes that the atomic model am has
        been fully determined; it simply loops through the (NIST) configurations and extracts all useful
        feature x-vectors (which also contain the NIST energy as y-vector). Two level lists
        (trainXyVectors::Array{AtomicFeatures.XyVector,1}, testXyVectors::Array{AtomicFeatures.XyVector,1})
        are returned.
"""
function generateFeatureVectors(confs::Array{Configuration,1}, am::AtomicFeatures.AtomicModel, probability::Float64,
                                nistLevels::Array{NistLevel,1})
    
    trainXyVectors = AtomicFeatures.XyVector[];   testXyVectors = AtomicFeatures.XyVector[] 
    
    for conf in confs
        DeepLearning.checkConfiguration(conf, am) 
        
        # Compute the level multiplet in the given atomic model and assign the LSJ-levels
        multiplet = DeepLearning.computeMultiplet(conf, am)
        lsjLevels = DeepLearning.extractLsjLevels(conf, multiplet, nistLevels::Array{NistLevel,1})
        DeepLearning.displayLsjLevels(lsjLevels)
        
        # Generate and select the feature x- (and y-) vectors
        xyVectors = DeepLearning.generateFeatureVectors(conf, multiplet:Multiplet, lsjLevels::Array{LsjLevel,1})        
        wxy       = DeepLearning.selectFeatureVectors(xyVectors, probability)
        append!(trainxyVectors, wxy[1]);    append!(testxyVectors, wxy[2])
    end

    return( (trainXyVectors, testXyVectors) )
end


"""
`DeepLearning.generateFeatureVectors(conf::Configuration, multiplet::Multiplet, lsjLevels::Array{LsjLevel,1})` 
    ... generates a list of feature vectors which are associated with the given configuration.
        This generation is based on 'script-like' code to enable the user to play with different assumptions
        and settings about the relevance of different features. Usually , the length of the feature vectors
        increases very rapidly with the number of "physical features" thar are included into the training
        and extraction of data.
        
        A list of xyVectors::AtomicFeatures.XyVector is returned.
"""
function generateFeatureVectors(conf::Configuration, multiplet::Multiplet, lsjLevels::Array{LsjLevel,1})
    xyVectors = AtomicFeatures.XyVector[]

    orbitals  = DeepLearning.extractOrbitals(conf::Configuration, am::AtomicFeatures.AtomicModel)
    
    # Append for each level in multiplet all desired features to the xVector
    for  (i, level)  in  enumerate(multiplet.Levels)
        xVector = Float64[]
        
        # Add shell occupations of the configuration
        append!(xVector,  AtomicFeatures.extractShellOccupations(am.shells, conf) )

        # Add mean occupation of the given LSJ level
        append!(xVector,  AtomicFeatures.extractMeanOccupationNumbers(am.subshells, xvector) )

        # Compute and add LSJ quantum numbers of level
        L = lsjLevels[i].L;    S = lsjLevels[i].S;    J = lsjLevels[i].J
        if  J != level.J  error("stop a")   end
        append!(xVector,  [Basics.twice(L)/2.0, Basics.twice(S)/2.0, Basics.twice(J)/2.0] )

        # Compute and add intermediate coupling quantum numbers of the n-leading CSF of the levels
        # (with zeros, if < n CSF are defined for the given level
        append!(xVector,  AtomicFeatures.extractIntermediateQN(am.subshells, 2, level) )

        # Add <r^k> of the given subshells
        ## append!(xVector,  AtomicFeatures.extractRkExpectation(am.subshells, -2, orbitals) )
        ## append!(xVector,  AtomicFeatures.extractRkExpectation(am.subshells, -1, orbitals) )
        append!(xVector,  AtomicFeatures.extractRkExpectation(am.subshells,  1, orbitals) )
        ## append!(xVector,  AtomicFeatures.extractRkExpectation(am.subshells,  2, orbitals) )

        # Add F^k (a,b) = R^k(a,b,a,b) of the given subshells
        ## append!(xVector,  AtomicFeatures.extractFkIntegrals(am.subshells, -2, orbitals) )
        ## append!(xVector,  AtomicFeatures.extractFkIntegrals(am.subshells, -1, orbitals) )
        append!(xVector,  AtomicFeatures.extractFkIntegrals(am.subshells,  0, orbitals) )
        ## append!(xVector,  AtomicFeatures.extractFkIntegrals(am.subshells,  1, orbitals) )
        append!(xVector,  AtomicFeatures.extractFkIntegrals(am.subshells,  2, orbitals) )
        append!(xVector,  AtomicFeatures.extractFkIntegrals(am.subshells,  4, orbitals) )

        # Add G^k (a,b) = R^k(a,b,b,a) of the given subshells
        ## append!(xVector,  AtomicFeatures.extractGkIntegrals(am.subshells, -2, orbitals) )
        ## append!(xVector,  AtomicFeatures.extractGkIntegrals(am.subshells, -1, orbitals) )
        append!(xVector,  AtomicFeatures.extractGkIntegrals(am.subshells,  0, orbitals) )
        ## append!(xVector,  AtomicFeatures.extractGkIntegrals(am.subshells,  1, orbitals) )
        append!(xVector,  AtomicFeatures.extractGkIntegrals(am.subshells,  2, orbitals) )
        append!(xVector,  AtomicFeatures.extractGkIntegrals(am.subshells,  4, orbitals) )
        
        
        # Append the newxVector to the desired list
        push!(xyVectors, AtomicFeatures.XyVector(L, S, J, level.energy, xVector))
    end
    
    println(">> $(length(xyVectors)) feature vectors have been generate for configuration $conf ")

    return( xyVectors )
end


"""
`DeepLearning.run(request::DeepLearning.LevelEstimationRequest, applic::DeepLearning.Application; output::Bool=true)`  
    ... to run a deep-learning level estimation (request) to estimate level energies for given configurations,
        and based on a given atomic model. The predicted energies are compared explicitly with the data availalbe 
        from the NIST Level database. A dict::Dict{String, Any} is returned if output=true, and nothing otherwise.
"""
function  run(request::DeepLearning.LevelEstimationRequest, applic::DeepLearning.Application; output::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    
    # Run the level estimation in turn for each configuration in the request
    for  conf in request.configs
        xyVectors = DeepLearning.generateFeatureVectors(conf, applic.atomicModel, request.nistLevels)
        lsjLevels = DeepLearning.applyFeatureVectors(xyVectors, application.neuralNetwork)
        # Print the configuration
        DeepLearning.displayLsjLevels(conf, lsjLevels)        
    end
    
    println("\n> Missing NIST level estimation complete ...")
    
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results )
end



"""
`DeepLearning.selectFeatureVectors(xyVectors::Array{AtomicFeatures.XyVector,1}, probability::Float64)`
    ... divides the set of feature vectors into two sets: for training and for testing due to the given 
        probability for selecting test data (probability <= 0.2). Two list of feature vectors 
        (trainXyVectors::Array{AtomicFeatures.XyVector,1}, testXyVectors::Array{AtomicFeatures.XyVector,1})
        are returned.
"""
function selectFeatureVectors(xyVectors::Array{AtomicFeatures.XyVector,1}, probability::Float64)
    trainXyVectors = AtomicFeatures.XyVector[];    testXyVectors = AtomicFeatures.XyVector[]  
    wrn = rand( length(xyVectors) );   nw = 0
    
    for xyVector in xyVectors
        if   xyVector.nistEnergy == 0.   continue    end
        nw = nw + 1
        if   wrn(nw) < probability    push!(testXyVectors, level)
        else                          push!(trainXyVectors, level)
        end 
    end

    return( (trainXyVectors, testXyVectors) )
end


"""
`DeepLearning.writeFeatureVectors(xyVectors::Array{AtomicFeatures.XyVector,1}, filename::String)` 
    ... writes out to a file(name) the xyVectors ins a format suitable for the training and test of neural networks.
        Nothing is returned.
"""
function writeFeatureVectors(xyVectors::Array{AtomicFeatures.XyVector,1}, filename::String)
    # Open a file to append the feature vectors
    ioFeatures = open(filename, "w")
    # write(ioFeatures, "Feature vector generated on $(today()): \n\n")
    
    for  xyVector in  xyVectors
        line = "  " * xyVector.xVector
        line = line * "        " * @sprintf("%.4e", xyVector.nistEnergy)
        write(ioFeatures, line )
    end
    
    close(ioFeatures)

    return(nothing)
end

    
    

