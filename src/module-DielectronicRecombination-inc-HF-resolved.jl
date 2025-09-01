

"""
`DielectronicRecombination.computeHyperfinePassages(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                                                    nm::Nuclear.Model, grid::Radial.Grid, empTreatment::EmpiricalTreatment, 
                                                    settings::DielectronicRecombination.Settings)`  
    ... to compute the data for all resonances (resonance lines) directly from the given multiplets of the initial-, intermediate- 
        and final states. It also enables one to (successively) include a set of corrections to the resonance strength to incorporate
        the contributions of shells that were not considered explicitly. 
        A list of resonances::Array{DielectronicRecombination.Resonance,1} is returned.
"""
function  computeHyperfinePassages(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                                   nm::Nuclear.Model, grid::Radial.Grid, empTreatment::EmpiricalTreatment,
                                   settings::DielectronicRecombination.Settings)
    
    error("Not yet: This code need to be adapted to hyperfine-resolved DR strength, etc.")
    
    passages = DielectronicRecombination.determinePassages(intermediateMultiplet, initialMultiplet, empTreatment, settings)
    # Display all selected resonances before the computations start
    if  settings.printBefore    DielectronicRecombination.displayPassages(stdout, passages)           end
    # Determine maximum (electron) energy and check for consistency of the grid
    maxEnergy = 0.;   for  passage in passages   maxEnergy = max(maxEnergy, passage.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    #
    # Calculate all amplitudes and requested properties; simply copy if the captureChannels have been computed before
    # Here, the selected set of "corrections" can also be considered for each passage.
    #
    if Distributed.nworkers() > 1
        # Distributed loop
        newPassages_ = @showprogress desc="Computing Passages ..." pmap(passage -> 
                        DielectronicRecombination.computeAmplitudesProperties(passage, finalMultiplet, nm, grid, nrContinuum, empTreatment, settings), passages)
        newPassages = convert(Vector{DielectronicRecombination.Passage}, newPassages_)
    else
    # Multithreading loop
        localPassages = [DielectronicRecombination.Passage[] for _ in 1:nthreads()]
        @threads for  p in eachindex(passages)
            newPassage = DielectronicRecombination.computeAmplitudesProperties(passages[p], finalMultiplet, nm, grid, nrContinuum, 
                                                                            empTreatment, settings) 
            push!(localPassages[threadid()], newPassage)
        end 
        newPassages = vcat(localPassages...)
    end
    #
    # Original implementation
    # newPassages = DielectronicRecombination.Passage[]; 
    # for  passage in passages
    #     newPassage = DielectronicRecombination.computeAmplitudesProperties(passage, finalMultiplet, nm, grid, nrContinuum, 
    #                                                                        empTreatment, settings) 
    #     push!( newPassages, newPassage)
    # end 

    # Add empirical passages to newPassages, if requested as correction
    if  empTreatment.doEmpiricalCorrections   &&   empTreatment.nUpperEmpirical > 0
        DielectronicRecombination.addEmpiricalPassages!(newPassages, empTreatment)      end
    # 
    # Calculate all corresponding resonance
    resonances = DielectronicRecombination.computeResonances(newPassages, settings)
    # Print all results to screen
    DielectronicRecombination.displayResults(stdout, resonances,  settings)
    DielectronicRecombination.displayRateCoefficients(stdout, resonances,  settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   DielectronicRecombination.displayResults(iostream, resonances,  settings)
                       DielectronicRecombination.displayRateCoefficients(iostream, resonances,  settings)    end
                
    return( newPassages )
end
