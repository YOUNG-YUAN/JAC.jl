
#==
++  September 2025: The Hfs module now comprises again all basic data structure (HfBasisVector, HfLevel, HfMultiplet)
    which we shall use to generate the representation of all hyperfine levels and to apply these levels in the evaluation
    of autoionization and photoemission amplitudes/rates.
    
    Here; I shall indicate the basic logic and definition of all main procedures that enables you to realize the 
    computation of hyperfine-resolved DR strength, etc. Inside of these procedure, however, I shall only provide some
    pseudo-code of different quality and ask for your help to "fill + test" this code. I shall indicate also 
    a list of "major steps", which can be realized and tested more or less independently.
    
    1) Generate from the final, intermediate and initial (electronic) multiplets the corresponding
       final, intermediate and initial hyperfine multiplets ::Hfs.HfMultiplet.
       While we shall include the full hyperfine interaction into the representation of the initialHfMultiplet,
       we can readily set the nuclear moments mu and Q simply to zero for the finalHfMultiplet and 
       intermediateHfMultiplet. This generation should be done by functions from the Hfs module.
       
    2) We first omit all empirical corrections ... but should keep the relevant code as comments inside,
       so that we can later re-activate this corrections, if appropriate.
       
    3) Determine the HfPassages by a function 
       DielectronicRecombination.determineHfPassages(intermediateHfMultiplet::Hfs.HfMutliplet, 
                                 initialHfMultiplet::Hfs.HfMutliplet, empTreatment, settings)
       The parameter empTreatment is currently not used.
       
    4) Make a new function DielectronicRecombination.displayHyperfinePassages(stdout, hfPassages) 
       that follows similar lines as DielectronicRecombination.displayPassages() but is adapted to
       hfPassages::Vector{HfPassage}.
       
    5) Make a function 
       DielectronicRecombination.computeHyperfineAmplitudes(hfPassage::HfPassage, 
                   finalHfMultiplet::HfMultiplet, nm, grid, nrContinuum, empTreatment, settings)
       The nm variable is redundant but likely useful.
       
    6) Make a function 
       DielectronicRecombination.computeHyperfineResonances(newHfPassages, settings)
       
    7) The two procedures 
       DielectronicRecombination.displayHyperfineResults(stdout, hfResonances,  settings)
       DielectronicRecombination.displayHyperfineRateCoefficients(stdout, hfResonances,  settings)
       need to be adapted; likely, it will be useful to introduce an optional, boolean flag
       hfResolved::Bool=true/false  in order to control whether all hyperfine-resolved data are 
       printed explicitly, or what is typically needed, to be comprised back into electronic
       resonances::Array{DR.Resonance,1} ... and for which you could use the existing procedure.
       
    It looks to me, this is (almost) all what you need to do; please, follow the style of the analogue 
    function --- this simplifies all our lives. Good luck.
    
    I presently include ...Hyperfine... in almost all function names to make this more explixit to you;
    we can later readily return to the previous names by using multiple dispatch. Let's first get the present 
    version running.

==#


#######################################################################################################################
#######################################################################################################################


"""
`struct  DielectronicRecombination.HfPassage`  
    ... defines a type for a dielectronic-recombination hyperfine passage, i.e. a (reduced) hyperfine pathways, 
        that include the definition of channels and their corresponding amplitudes for the individual 
        i --> m resonances, whereas the subsequent radiative stabilization is considered only later.

    + initialHfLevel       ::Hfs.HfLevel        ... initial-(state) hyperfine level
    + intermediateHfLevel  ::Hfs.HfLevel        ... intermediate-(state) hyperfine level
    + electronEnergy       ::Float64            ... energy of the (incoming, captured) electron
    + captureRate          ::Float64            ... rate for the electron capture (Auger rate)
    + photonRate           ::EmProperty         ... rate for the photon emission
    + reducedStrength      ::EmProperty              
        ... reduced (hyperfine) resonance strength Sum_f S(i -> d -> f) * Gamma_d of this hyperfine passage; 
            this reduced strength does not require the knowledge of Gamma_d for the individual hyperfine passage.
    + captureChannels   ::Array{AutoIonization.Channel,1}   
        ... List of |i> -->  |n>   dielectronic (Auger) capture channels which purely refer to the electronic levels.
"""
struct  HfPassage
    initialHfLevel         ::Hfs.HfLevel
    intermediateHfLevel    ::Hfs.HfLevel
    electronEnergy         ::Float64
    captureRate            ::Float64
    photonRate             ::EmProperty
    reducedStrength        ::EmProperty
    captureChannels        ::Array{AutoIonization.Channel,1} 
end 


"""
`DielectronicRecombination.HfPassage()`  
    ... constructor for an 'empty' instance of a dielectronic recombination hyperfine passage between a specified 
        initial and intermediate level.
"""
function HfPassage()
    em = EmProperty(0., 0.)
    HfPassage(Hfs.HfLevel(), Hfs.HfLevel(), 0., 0., em, em, AutoIonization.Channel[])
end


# `Base.show(io::IO, passage::DielectronicRecombination.HfPassage)`  
#   ... prepares a proper printout of the variable passage::DielectronicRecombination.HfPassage.
function Base.show(io::IO, passage::DielectronicRecombination.HfPassage) 
    println(io, "initialHfLevel:             $(passage.initialHfLevel)  ")
    println(io, "intermediateHfLevel:        $(passage.intermediateHfLevel)  ")
    println(io, "electronEnergy:             $(passage.electronEnergy)  ")
    println(io, "captureRate:                $(passage.captureRate)  ")
    println(io, "photonRate:                 $(passage.photonRate)  ")
    println(io, "reducedStrength:            $(passage.reducedStrength)  ")
    println(io, "captureChannels:            $(passage.captureChannels)  ")
end


"""
`struct  DielectronicRecombination.HfResonance`  
    ... defines a type for a dielectronic hyperfine resonance as defined by a given initial and resonance 
        hyprfine levels but by summing over all final (hyperfine) levels.

    + initialHfLevel       ::Hfs.HfLevel       ... initial-(state) hyperfine level
    + intermediateHfLevel  ::Level             ... intermediate-(state) hyperfine level
    + resonanceEnergy      ::Float64           ... energy of the resonance w.r.t. the inital-state
    + resonanceStrength    ::EmProperty        ... strength of this resonance due to the stabilization into any of the allowed final levels.
    + captureRate          ::Float64           ... capture (Auger) rate to form the intermediate (hyperfine) resonance, starting 
                                                   from the initial level.
    + augerRate            ::Float64           ... total (Auger) rate for an electron emission of the intermediate hyperfine resonance
    + photonRate           ::EmProperty        ... total photon rate for a photon emission, i.e. for stabilization.
"""
struct  HfResonance
    initialHfLevel         ::Level
    intermediateHfLevel    ::Level
    resonanceEnergy        ::Float64 
    resonanceStrength      ::EmProperty
    captureRate            ::Float64
    augerRate              ::Float64
    photonRate             ::EmProperty
end 


"""
`DielectronicRecombination.HfResonance()`  
    ... constructor for an 'empty' instance of a dielectronic hyperfine resonance as defined by a given initial 
        and resonance hyperfine level but by summing over all final hyperfine levels.
"""
function HfResonance()
    em = EmProperty(0., 0.)
    HfResonance(Hfs.HfLevel, Hfs.HfLevel, 0., em, 0., 0., em)
end


# `Base.show(io::IO, hfResonance::DielectronicRecombination.HfResonance)`  ... prepares a proper printout of the variable resonance::DielectronicRecombination.HfResonance.
function Base.show(io::IO, resonance::DielectronicRecombination.HfResonance) 
    println(io, "initialHfLevel:             $(resonance.initialLevel)  ")
    println(io, "intermediateHfLevel:        $(resonance.intermediateLevel)  ")
    println(io, "resonanceEnergy:            $(resonance.resonanceEnergy)  ")
    println(io, "resonanceStrength:          $(resonance.resonanceStrength)  ")
    println(io, "captureRate:                $(resonance.captureRate)  ")
    println(io, "augerRate:                  $(resonance.augerRate)  ")
    println(io, "photonRate:                 $(resonance.photonRate)  ")
end


#######################################################################################################################
#######################################################################################################################


"""
`DielectronicRecombination.computeHyperfinePassages(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                                                    nm::Nuclear.Model, grid::Radial.Grid, empTreatment::EmpiricalTreatment, 
                                                    settings::DielectronicRecombination.Settings)`  
    ... to compute the data for all resonances (resonance lines) directly from the given multiplets of the initial-, intermediate- 
        and final states. The computation of hyperfine-resolved passages and resonance strength assumes, however, that all level are 
        hfLevels, i.e. hyperfine-resolved levels. This applies especially for the initial levels, which must include the full 
        hyperfine Hamiltonian and splitting in their representation. For the intermediate and final levels, in contrast, we still
        assume a hyperfine-resolution (i.e. the use hfLevel's) but shall set the nuclear magnetic-dipole and electric-quadrupole
        simply to zero. Hence, the intermediate and final levels are completely degenerate and later be comprised into
        (electronic) resonances and observables. A list of resonances::Array{DielectronicRecombination.HfResonance,1} is 
        returned.
        
        The function is prepared also to (successively) include a set of corrections to the resonance strength to incorporate
        the contributions of shells that were not considered explicitly. However, this branch of the code is not supported in 
        the present version of the code.    
"""
function  computeHyperfinePassages(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                                   nm::Nuclear.Model, grid::Radial.Grid, empTreatment::EmpiricalTreatment,
                                   settings::DielectronicRecombination.Settings)
    # First, we generate from the given multiplets the corresponding hfMultiplet::HfMultiplet's by taking the interaction
    # with the nuclear moments into account.
    nmZero = Nuclear.Model(nm; mu=0., Q=0.)       # This is the same nuclear model but with zero mu = Q = 0 nuclear moments.
    initialHfMultiplet      = Hfs.generateHfMultiplet(initialMultiplet, nm)
    intermediateHfMultiplet = Hfs.generateHfMultiplet(intermediateMultiplet, nmZero)
    finalHfMultiplet        = Hfs.generateHfMultiplet(finalMultiplet, nmZero)
    
    hfPassages = DielectronicRecombination.determineHyperfinePassages(intermediateHfMultiplet, initialHfMultiplet, 
                                                                      empTreatment, settings)
    # Display all selected resonances before the computations start
    if  settings.printBefore    DielectronicRecombination.displayHyperfinePassages(stdout, passages)    end
    # Determine maximum (electron) energy and check for consistency of the grid
    maxEnergy = 0.;   for  passage in hfPassages   maxEnergy = max(maxEnergy, passage.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    
    # Calculate all amplitudes and requested properties; simply copy if the captureChannels have been computed before
    # Here, the selected set of "corrections" can also be considered for each passage.
    #
    if Distributed.nworkers() > 1
        # Distributed loop
        newHfPassages_ = @showprogress desc="Computing Passages ..." pmap(hfPassage -> 
                        DielectronicRecombination.computeHyperfineAmplitudes(hfPassage, finalHfMultiplet, nm, grid, 
                                                                             nrContinuum, empTreatment, settings), hfPassages)
        newHfPassages = convert(Vector{DielectronicRecombination.HfPassage}, newHfPassages_)
    else
        # Multithreading loop
        localHfPassages = [DielectronicRecombination.HfPassage[] for _ in 1:nthreads()]
        @threads for  p in eachindex(hfPassages)
            newHfPassage = DielectronicRecombination.computeHyperfineAmplitudes(hfPassages[p], finalMultiplet, nm, grid, 
                                                                                nrContinuum, empTreatment, settings) 
            push!(localHfPassages[threadid()], newHfPassage)
        end 
        newHfPassages = vcat(localHfPassages...)
    end

    #== Add empirical passages to newPassages, if requested as correction
    if  empTreatment.doEmpiricalCorrections   &&   empTreatment.nUpperEmpirical > 0
        DielectronicRecombination.addEmpiricalPassages!(newPassages, empTreatment)      end  ==#
    # 
    # Calculate all corresponding resonance
    hfResonances = DielectronicRecombination.computeHyperfineResonances(newHfPassages, settings)
    # Print all results to screen
    DielectronicRecombination.displayHyperfineResults(stdout, hfResonances,  settings)
    DielectronicRecombination.displayHyperfineRateCoefficients(stdout, hfResonances,  settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   DielectronicRecombination.displayHyperfineResults(stdout, hfResonances,  settings)
                       DielectronicRecombination.displayHyperfineRateCoefficients(stdout, hfResonances,  settings)    end
                
    return( newHfPassages )
end
    

"""
`DielectronicRecombination.computeHyperfineAmplitudes(hfPassage::DielectronicRecombination.HfPassage, 
                           finalMultiplet::Hfs.HfMultiplet, grid::Radial.Grid, nrContinuum::Int64, 
                           empTreatment::EmpiricalTreatment, settings::DielectronicRecombination.Settings)` 
    ... to compute all amplitudes and properties of the given hyperfine passage; a hfPassage::DielectronicRecombination.HfPassage
        is returned for which the amplitudes and properties have now been evaluated. No parameter nm::Nuclear.Model
        is provided here since the information about the nuclear (isomeric) states need to be taken from the (basis) of the 
        hyperfine levels themselves.
"""
function  computecomputeHyperfineAmplitudes(hfPassage::DielectronicRecombination.HfPassage, 
                                 finalMultiplet::Hfs.HfMultiplet, grid::Radial.Grid, nrContinuum::Int64, 
                                 empTreatment::EmpiricalTreatment, settings::DielectronicRecombination.Settings)
    rateA = 0.
    # ...
   
   
    return( hfPassage )
end


"""
`DielectronicRecombination.computeHyperfineRateCoefficient(hfResonance::DielectronicRecombination.HfResonance, temp::Float64)`  
    ... computes for a delta-like resonance the DR rate coefficient alpha_d (i, Te) from the given resonance strength
        and temperature [K], and for both, Coulomb and Babushkin gauge. All values are directly returned in [cm^3/s].
        An alphaDR::EmProperty is returned. ... We might also first go back to the Resonance data type and to simply 
        use the existing function.
"""
function computeRateHyperfineCoefficient(hfResonance::DielectronicRecombination.HfResonance, temp::Float64)
                
    return( alphaDR )
end


"""
`DielectronicRecombination.computeHyperfineResonances(hfPassages::Array{DielectronicRecombination.HfPassage,1}, 
                                                      settings::DielectronicRecombination.Settings)`  
    ... to compute the data for all hyperfine resonances (hyperfine resonance lines) as defined by the given 
        hyperfine passages and and settings. For hyperfine-resolved spectra, we shall support only the computation of 
        hyperfine passages and extract the (normal) resonances from these hyperfine passage. A list o resonances::Array{DielectronicRecombination.Resonance,1} is returned.
"""
function  computeHyperfineResonances(hfPassages::Array{DielectronicRecombination.HfPassage,1}, 
                                     settings::DielectronicRecombination.Settings)
    hfResonances = DielectronicRecombination.HfResonance[]
    
    return( hfResonances )
end


"""
`DielectronicRecombination.determineHyperfineCaptureChannels(intermediateLevel::Hfs.HfLevel, initialLevel::Hfs.HfLevel, 
                                                             settings::DielectronicRecombination.Settings)` 
    ... to determine a list of AutoIonization.Channel for a (Auger) capture transitions from the initial to an 
        intermediate level, and by taking into account the particular settings of for this computation;  
        an Array{AutoIonization.Channel,1} is returned. The capture/autoionization channels are still purely
        electronic (not hyperfine) channels but need to be properly combined with the hyperfine notation.
"""
function determineHyperfineCaptureChannels(intermediateLevel::Hfs.HfLevel, initialLevel::Hfs.HfLevel, 
                                           settings::DielectronicRecombination.Settings)
    channels = AutoIonization.Channel[];  
    #== Need to be adapted.
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity)
    symn = LevelSymmetry(intermediateLevel.J, intermediateLevel.parity)
    kappaList = AngularMomentum.allowedKappaSymmetries(symi, symn)
    for  kappa in kappaList
        push!( channels, AutoIonization.Channel(kappa, symn, 0., Complex(0.)) )
    end  ==#

    return( channels )  
end


"""
`DielectronicRecombination.determineHyperfinePhotonChannels(finalLevel::Hfs.HfLevel, intermediateLevel::Hfs.HfLevel, 
                                                            settings::DielectronicRecombination.Settings)` 
    ... to determine a list of PhotoEmission.Channel for the photon transitions from the intermediate and to a final level, and by 
        taking into account the particular settings of for this computation;  an Array{PhotoEmission.Channel,1} is returned.
        The photoemission channels are still purely electronic (not hyperfine) channels but need to be properly combined 
        with the hyperfine notation.
"""
function determinePhotonChannels(finalLevel::Hfs.HfLevel, intermediateLevel::Hfs.HfLevel, settings::DielectronicRecombination.Settings)
    channels = PhotoEmission.Channel[];  
    #== Need to be adapted.
    symn = LevelSymmetry(intermediateLevel.J, intermediateLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
    for  mp in settings.multipoles
        if   AngularMomentum.isAllowedMultipole(symn, mp, symf)
            hasMagnetic = false
            for  gauge in settings.gauges
                # Include further restrictions if appropriate
                if     string(mp)[1] == 'E'  &&   gauge == UseCoulomb      push!(channels, PhotoEmission.Channel(mp, Basics.Coulomb,   0.) )
                elseif string(mp)[1] == 'E'  &&   gauge == UseBabushkin    push!(channels, PhotoEmission.Channel(mp, Basics.Babushkin, 0.) )  
                elseif string(mp)[1] == 'M'  &&   !(hasMagnetic)           push!(channels, PhotoEmission.Channel(mp, Basics.Magnetic,  0.) );
                                                    hasMagnetic = true; 
                end 
            end
        end
    end   ==#

    return( channels )  
end


"""
`DielectronicRecombination.determineHyperfinePassages(intermediateMultiplet::Hfs.HfMultiplet, initialMultiplet::Hfs.HfMultiplet, 
                                             empTreatment::EmpiricalTreatment, settings::DielectronicRecombination.Settings)`  
    ... to determine a list of dielectronic-recombination hyperfine resonances between the levels from the given 
        (hyperfine) initial- and intermediate- states, whereas the final states are considered "on-fly"; the particular 
        selections and settings for this computation are taken into account; an Array{DielectronicRecombination.HfPasssage,1} 
        is returned. Apart from the level specification, all physical properties are set to zero during the 
        initialization process. The parameter empTreatment::EmpiricalTreatment is obsolete here but might be "activated"
        in the future.
"""
function  determineHyperfinePassages(intermediateMultiplet::Hfs.HfMultiplet, initialMultiplet::Hfs.HfMultiplet,
                                     empTreatment::EmpiricalTreatment, settings::DielectronicRecombination.Settings)
    hfPassages = DielectronicRecombination.HfPassage[]
    #== electronEnergyShift = Defaults.convertUnits("energy: to atomic", settings.electronEnergyShift)
    @warn("No pathway selection is considered, if settings.calcOnlyPassages=true.")
    #
    for  iLevel  in  initialMultiplet.levels
        for  nLevel  in  intermediateMultiplet.levels
            eEnergy = nLevel.energy - iLevel.energy + electronEnergyShift
            if  eEnergy < 0.  ||   eEnergy <  empTreatment.resonanceEnergyMin  ||   
                                   eEnergy >  empTreatment.resonanceEnergyMax  continue    end
            cChannels = DielectronicRecombination.determineCaptureChannels(nLevel, iLevel, settings) 
            push!( passages, DielectronicRecombination.Passage(iLevel, nLevel, eEnergy, 0., EmProperty(0., 0.), 
                                                                EmProperty(0., 0.), cChannels) )
        end
    end  ==#
    return( hfPassages )
end


"""
`DielectronicRecombination.displayHyperfinePassages(stream::IO, hfPassages::Array{DielectronicRecombination.HfPassage,1})`  
    ... to display a list of (hyperfine) passages and channels that have been selected due to the prior settings. 
        A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayHyperfinePassages(stream::IO, hfPassages::Array{DielectronicRecombination.HfPassage,1})
    nx = 120
    println(stream, " ")
    println(stream, "  Selected dielectronic-recombination hyperfine passages:")
    #== Need to be adapted !
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "     ";   sb = "     "
    sa = sa * TableStrings.center(16, "Levels"; na=4);            sb = sb * TableStrings.center(16, "i  --  m"; na=4);          
    sa = sa * TableStrings.center(16, "J^P symmetries"; na=3);    sb = sb * TableStrings.center(16, "i  --  m"; na=3);
    sa = sa * TableStrings.center(18, "Energies  " * TableStrings.inUnits("energy"); na=5);              
    sb = sb * TableStrings.center(18, "electron     "; na=5)
    sa = sa * TableStrings.flushleft(57, "List of kappas and total symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(57, "partial (total J^P)                  "; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  passage in passages
        sa  = "  ";     isym = LevelSymmetry( passage.initialLevel.J,      passage.initialLevel.parity)
                        msym = LevelSymmetry( passage.intermediateLevel.J, passage.intermediateLevel.parity)
        sa = sa * TableStrings.center(17, TableStrings.levels_if(passage.initialLevel.index, passage.intermediateLevel.index); na=7)
        sa = sa * TableStrings.center(17, TableStrings.symmetries_if(isym, msym);  na=4)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", passage.electronEnergy)) * "        "
        kappaSymmetryList = Tuple{Int64,LevelSymmetry}[]
        for  cChannel in passage.captureChannels
            push!( kappaSymmetryList, (cChannel.kappa, cChannel.symmetry) )
        end
        wa = TableStrings.kappaSymmetryTupelList(85, kappaSymmetryList)
        if  length(wa) > 0    sb = sa * wa[1];    println(stream,  sb )    end  
        for  i = 2:length(wa)
            sb = TableStrings.hBlank( length(sa) ) * wa[i];    println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "\n>> A total of $(length(passages)) dielectronic-recombination passages will be calculated. \n")
    ==#
    #
    return( nothing )
end


"""
`DielectronicRecombination.displayResults(stream::IO, hfResonances::Array{DielectronicRecombination.HfResonance,1},
                                          settings::DielectronicRecombination.Settings)`  
    ... to list all results for the hyperfine resonances. A neat table is printed but nothing is returned otherwise.
"""
function  displayResults(stream::IO, hfResonances::Array{DielectronicRecombination.HfResonance,1},
                                     settings::DielectronicRecombination.Settings)
    nx = 160
    #==   Need to be adapted.
    println(stream, " ")
    println(stream, "  Total Auger rates, radiative rates and resonance strengths:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-m"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--m"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(14, "Energy"   ; na=2);               
    sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(42, "Auger rate     Cou -- rad. rates -- Bab"; na=1);       
    sb = sb * TableStrings.center(16, TableStrings.inUnits("rate"); na=1)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("rate"); na=0)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("rate"); na=6)
    sa = sa * TableStrings.center(30, "Cou -- res. strength -- Bab"; na=3);       
    sb = sb * TableStrings.center(12, TableStrings.inUnits("strength");  na=0)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("strength");  na=2)
    sa = sa * TableStrings.center(18, "Widths Gamma_m"; na=2);       
    sb = sb * TableStrings.center(16, TableStrings.inUnits("energy"); na=6)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  resonance in resonances
        sa  = "";      isym = LevelSymmetry( resonance.initialLevel.J,      resonance.initialLevel.parity)
                       msym = LevelSymmetry( resonance.intermediateLevel.J, resonance.intermediateLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(resonance.initialLevel.index, resonance.intermediateLevel.index); na=4)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, msym);  na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", resonance.resonanceEnergy))          * "      "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", resonance.augerRate))                  * "      "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", resonance.photonRate.Coulomb))         * "  "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("rate: from atomic", resonance.photonRate.Babushkin))       * "        "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("strength: from atomic", resonance.resonanceStrength.Coulomb))    * "  "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("strength: from atomic", resonance.resonanceStrength.Babushkin))  * "     "
        wa = resonance.augerRate + resonance.photonRate.Coulomb 
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", wa))                                 * "   "
        wa = resonance.augerRate + resonance.photonRate.Babushkin 
        sa = sa * "(" * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", wa)) * ")"                     * "   "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    ==#
    #
    return( nothing )
end


"""
`DielectronicRecombination.displayRateCoefficients(stream::IO, hfResonances::Array{DielectronicRecombination.HfResonance,1},
                                                   settings::DielectronicRecombination.Settings)`  
    ... to list, if settings.calcRateAlpha, all rate coefficients for the selected temperatures. Both, the individual as well as
        the total DR plasma rate coefficients are printed in neat tables, though nothing is returned otherwise.
"""
function  displayRateCoefficients(stream::IO, hfResonances::Array{DielectronicRecombination.HfResonance,1},
                                  settings::DielectronicRecombination.Settings)
    ntemps = length(settings.temperatures)
    #==  Need to be adapted.
    if  !settings.calcRateAlpha  ||  ntemps == 0     return(nothing)     end
    #
    nx = 54 + 17 * min(ntemps, 7)
    println(stream, " ")
    println(stream, "  Rate coefficients for delta-like resonances [cm^3/s]:        ... all results in Babushkin gauge")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-m"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--m"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(14, "Energy"   ; na=2);               
    sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=2)
    for  nt = 1:min(ntemps, 7)
        sa = sa * TableStrings.center(14, "T = " * @sprintf("%.2e", settings.temperatures[nt]); na=3);       
        sb = sb * TableStrings.center(14, "[K]"; na=3)
    end
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  resonance in resonances
        sa  = "";      isym = LevelSymmetry( resonance.initialLevel.J,      resonance.initialLevel.parity)
                        msym = LevelSymmetry( resonance.intermediateLevel.J, resonance.intermediateLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(resonance.initialLevel.index, resonance.intermediateLevel.index); na=4)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, msym);  na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", resonance.resonanceEnergy))    * "      "
        for  nt = 1:min(ntemps, 7)
            alphaDR      = DielectronicRecombination.computeRateCoefficient(resonance, settings.temperatures[nt])
            sa = sa * @sprintf("%.4e", alphaDR.Babushkin)  * "       "
        end
        println(stream, sa)
    end
    #
    println(stream, "  ")
    sa = "       alpha^DR (T, i; Coulomb gauge):                      " 
    sb = "       alpha^DR (T, i; Babushkin gauge):                    " 
    for  nt = 1:min(ntemps, 7) 
        alphaDRtotal = EmProperty(0.)
        for  resonance in resonances
            alphaDR      = DielectronicRecombination.computeRateCoefficient(resonance, settings.temperatures[nt])
            alphaDRtotal = alphaDRtotal + alphaDR
        end
        sa = sa * @sprintf("%.4e", alphaDRtotal.Coulomb)    * "       "
        sb = sb * @sprintf("%.4e", alphaDRtotal.Babushkin)  * "       "
    end
    println(stream, sa);    println(stream, sb)
    println(stream, "  ", TableStrings.hLine(nx))
    ==#
    #
    return( nothing )
end

