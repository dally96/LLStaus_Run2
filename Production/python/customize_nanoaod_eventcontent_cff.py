import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.nano_eras_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer

from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.taus_cff import *
from PhysicsTools.NanoAOD.muons_cff import *
from PhysicsTools.NanoAOD.jetsAK4_CHS_cff import *
#from PhysicsTools.NanoAOD.custom_muon_cff import *
#from PhysicsTools.NanoAOD.jets_cff import *


def customize_process_and_associate(process, isMC, disTauTagOutputOpt = 1) :
    # Lost tracks
    process.lostTrackTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("lostTracks"),
        #cut = cms.string(""),
        cut = cms.string("pt > 1"),
        name= cms.string("LostTrack"),
        doc = cms.string("Lost tracks"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            CandVars,
        )
    )

    process.disMuonTable = simpleCandidateFlatTableProducer.clone(
        src = cms.InputTag("slimmedDisplacedMuons"),
        name = cms.string("DisMuon"),
        doc = cms.string("Displaced Muon Collection"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(CandVars,
            ptErr   = Var("bestTrack().ptError()", float, doc = "ptError of the muon track", precision=6),
            tunepRelPt = Var("tunePMuonBestTrack().pt/pt",float,doc="TuneP relative pt, tunePpt/pt",precision=6),
            dz = Var("dB('PVDZ')",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
            dzErr = Var("abs(edB('PVDZ'))",float,doc="dz uncertainty, in cm",precision=6),
            dxybs = Var("dB('BS2D')",float,doc="dxy (with sign) wrt the beam spot, in cm",precision=10),
            dxy = Var("dB('PV2D')",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
            dxyErr = Var("edB('PV2D')",float,doc="dxy uncertainty, in cm",precision=6),
            trkChi2 = Var("? globalTrack().isNonnull() ? globalTrack().normalizedChi2() : ? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().normalizedChi2() : -99",float,doc="Normalized Chi Square from either globalTrack or innerTrack "),
            muonHits = Var("? globalTrack().isNonnull() ? globalTrack().hitPattern().numberOfValidMuonHits() : ?  innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().numberOfValidMuonHits() :-99",float,doc="Number of valid Muon Hits from either globalTrack or innerTrack"),
            pixelHits = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().numberOfValidPixelHits() : -99", float, doc="Numbr of valid pixel hits"),
            validFraction = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().validFraction() : -99", float, doc="Inner Track Valid Fraction"),
            positionChi2 = Var("combinedQuality().chi2LocalPosition", float, doc="chi2 Local Position"),
            trkKink = Var("combinedQuality().trkKink", float, doc="Track Kink"),
            ip3d = Var("abs(dB('PV3D'))",float,doc="3D impact parameter wrt first PV, in cm",precision=10),
            sip3d = Var("abs(dB('PV3D')/edB('PV3D'))",float,doc="3D impact parameter significance wrt first PV",precision=10),
            segmentComp   = Var("segmentCompatibility()", float, doc = "muon segment compatibility", precision=14), # keep higher precision since people have cuts with 3 digits on this
            nStations = Var("numberOfMatchedStations", "uint8", doc = "number of matched stations with default arbitration (segment & track)"),
            nTrackerLayers = Var("?track.isNonnull?innerTrack().hitPattern().trackerLayersWithMeasurement():0", "uint8", doc = "number of layers in the tracker"),
            highPurity = Var("?track.isNonnull?innerTrack().quality('highPurity'):0", bool, doc = "inner track is high purity"),
            jetIdx = Var("?hasUserCand('jet')?userCand('jet').key():-1", "int16", doc="index of the associated jet (-1 if none)"),
            svIdx = Var("?hasUserCand('vertex')?userCand('vertex').key():-1", "int16", doc="index of matching secondary vertex"),
            tkRelIso = Var("isolationR03().sumPt/tunePMuonBestTrack().pt",float,doc="Tracker-based relative isolation dR=0.3 for highPt, trkIso/tunePpt",precision=6),
            #miniPFRelIso_chg = Var("userFloat('miniIsoChg')/pt",float,doc="mini PF relative isolation, charged component"),
            #miniPFRelIso_all = Var("userFloat('miniIsoAll')/pt",float,doc="mini PF relative isolation, total (with scaled rho*EA PU corrections)"),
            pfRelIso03_chg = Var("pfIsolationR03().sumChargedHadronPt/pt",float,doc="PF relative isolation dR=0.3, charged component"),
            pfRelIso03_all = Var("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.3, total (deltaBeta corrections)"),
            pfRelIso04_all = Var("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.4, total (deltaBeta corrections)"),
            #jetRelIso = Var("?userCand('jetForLepJetVar').isNonnull()?(1./userFloat('ptRatio'))-1.:(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt",float,doc="Relative isolation in matched jet (1/ptRatio-1, pfRelIso04_all if no matched jet)",precision=8),
            #jetPtRelv2 = Var("?userCand('jetForLepJetVar').isNonnull()?userFloat('ptRel'):0",float,doc="Relative momentum of the lepton with respect to the closest jet after subtracting the lepton",precision=8),
            tightCharge = Var("?(muonBestTrack().ptError()/muonBestTrack().pt() < 0.2)?2:0", "uint8", doc="Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass)"),
            looseId  = Var("passed('CutBasedIdLoose')",bool, doc="muon is loose muon"),
            isPFcand = Var("isPFMuon",bool,doc="muon is PF candidate"),
            isGlobal = Var("isGlobalMuon",bool,doc="muon is global muon"),
            isTracker = Var("isTrackerMuon",bool,doc="muon is tracker muon"),
            isStandalone = Var("isStandAloneMuon",bool,doc="muon is a standalone muon"),
            mediumId = Var("passed('CutBasedIdMedium')",bool,doc="cut-based ID, medium WP"),
            mediumPromptId = Var("passed('CutBasedIdMediumPrompt')",bool,doc="cut-based ID, medium prompt WP"),
            tightId = Var("passed('CutBasedIdTight')",bool,doc="cut-based ID, tight WP"),
            softId = Var("passed('SoftCutBasedId')",bool,doc="soft cut-based ID"),
            softMvaId = Var("passed('SoftMvaId')",bool,doc="soft MVA ID"),
            softMva = Var("softMvaValue()",float,doc="soft MVA ID score",precision=6),
            highPtId = Var("?passed('CutBasedIdGlobalHighPt')?2:passed('CutBasedIdTrkHighPt')","uint8",doc="high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)"),
            pfIsoId = Var("passed('PFIsoVeryLoose')+passed('PFIsoLoose')+passed('PFIsoMedium')+passed('PFIsoTight')+passed('PFIsoVeryTight')+passed('PFIsoVeryVeryTight')","uint8",doc="PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)"),
            tkIsoId = Var("?passed('TkIsoTight')?2:passed('TkIsoLoose')","uint8",doc="TkIso ID (1=TkIsoLoose, 2=TkIsoTight)"),
            miniIsoId = Var("passed('MiniIsoLoose')+passed('MiniIsoMedium')+passed('MiniIsoTight')+passed('MiniIsoVeryTight')","uint8",doc="MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)"),
            mvaMuID = Var("mvaIDValue()",float,doc="MVA-based ID score ",precision=6),
            #mvaMuID_WP = Var("userFloat('mvaIDMuon_wpMedium') + userFloat('mvaIDMuon_wpTight')","uint8",doc="MVA-based ID selector WPs (1=MVAIDwpMedium,2=MVAIDwpTight)"),
            multiIsoId = Var("?passed('MultiIsoMedium')?2:passed('MultiIsoLoose')","uint8",doc="MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium)"),
            puppiIsoId = Var("passed('PuppiIsoLoose')+passed('PuppiIsoMedium')+passed('PuppiIsoTight')", "uint8", doc="PuppiIsoId from miniAOD selector (1=Loose, 2=Medium, 3=Tight)"),
            triggerIdLoose = Var("passed('TriggerIdLoose')",bool,doc="TriggerIdLoose ID"),
            inTimeMuon = Var("passed('InTimeMuon')",bool,doc="inTimeMuon ID"),
            #jetNDauCharged = Var("?userCand('jetForLepJetVar').isNonnull()?userFloat('jetNDauChargedMVASel'):0", "uint8", doc="number of charged daughters of the closest jet"),
            ),
    )

    process.disMuonsMCMatchForTable = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
        src         = process.disMuonTable.src,                         # final reco collection
        matched     = cms.InputTag("finalGenParticles"),     # final mc-truth particle collection
        mcPdgId     = cms.vint32(13),               # one or more PDG ID (13 = mu); absolute values (see below)
        checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge
        mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
        maxDeltaR   = cms.double(0.3),              # Minimum deltaR for the match
        maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
        resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
        resolveByMatchQuality = cms.bool(True),    # False = just match input in order; True = pick lowest deltaR pair first
    )

    process.disMuonMCTable = cms.EDProducer("CandMCMatchTableProducer",
        src     = process.disMuonTable.src,
        mcMap   = cms.InputTag("disMuonsMCMatchForTable"),
        objName = process.disMuonTable.name,
        objType = cms.string("Muon"), #cms.string("Muon"),
        branchName = cms.string("genPart"),
        docString = cms.string("MC matching to status==1 muons"),
    )
    
    #disMuonTask = cms.Task(disMuonTable)
    #disMuonMCTask = cms.Task(disMuonsMCMatchForTable, disMuonMCTable)

    ### Custom Muon branches (https://github.com/cms-sw/cmssw/blob/9fca027cd9b380c6d9f4244b56f3e2c9ef2899db/PhysicsTools/NanoAOD/python/custom_muon_cff.py)
    #muonWithVariables = "muonWithVariables"
    #setattr(process, muonWithVariables, cms.EDProducer("SimpleTrackFlatTableProducer",
    #                muonSrc=cms.InputTag("slimmedMuons"),
    #                vertexSrc=cms.InputTag("offlineSlimmedPrimaryVertices"),
    #                trkSrc=cms.InputTag("pfTracks"),
    #                )
    #)
    #getattr(process,"muonTask").add(getattr(process,muonWithVariables))
    
    #process.slimmedMuonsUpdated.src = cms.InputTag("muonWithVariables")   
    
    #process.disMuonTable = muonTable.clone()
    #process.disMuonTable.src = cms.InputTag("slimmedDisplacedMuons")
    #process.disMuonTable.name = cms.string("DisMuon")
    #process.disMuonTable.doc = cms.string("Slimmed Displaced Muon Collection") 
    #process.disMuonTable.variables.trkChi2 = Var("? globalTrack().isNonnull() ? globalTrack().normalizedChi2() : ? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().normalizedChi2() : -99",float,doc="Normalized Chi Square from either globalTrack or innerTrack ")
    #process.disMuonTable.variables.muonHits = Var("? globalTrack().isNonnull() ? globalTrack().hitPattern().numberOfValidMuonHits() : ?  innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().numberOfValidMuonHits() :-99",float,doc="Number of valid Muon Hits from either globalTrack or innerTrack")
    #process.disMuonTable.variables.pixelHits = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().numberOfValidPixelHits() : -99", float, doc="Numbr of valid pixel hits")
    #process.disMuonTable.variables.miniPFRelIso_chg = Var("-99", float, doc="Doesn't exist for DisMuons")
    #process.disMuonTable.variables.miniPFRelIso_all = Var("-99", float, doc="Doesn't exist for DisMuons")
    #process.disMuonTable.variables.mvaMuID_WP = Var("-99", float, doc="Doesn't exist for DisMuons")
    #(run2_nanoAOD_106Xv2 | run3_nanoAOD_122).toModify(process.disMuonTable.variables,mvaMuID=None).toModify(process.disMuonTable.variables, mvaMuID = Var("99", float, doc="MVA-based ID score",precision=6))


    #process.disMuonsMCMatchForTable = muonsMCMatchForTable.clone()
    #process.disMuonsMCMatchForTable.src = process.disMuonTable.src
  
    #process.disMuonMCTable = muonMCTable.clone()
    #process.disMuonMCTable.src = process.disMuonTable.src
    #process.disMuonMCTable.objName = process.disMuonTable.name
    #process.disMuonMCTable.mcMap = cms.InputTag("disMuonsMCMatchForTable")
    #
    process.disMuonMCTask = cms.Task(process.disMuonsMCMatchForTable, process.disMuonMCTable)  
    process.disMuonTablesTask = cms.Task(process.disMuonTable)

    myMuonTable = muonTable.clone()
    
    #StandAlone Variables
    myMuonTable.variables.standalonePt = Var("? standAloneMuon().isNonnull() ? standAloneMuon().pt() : -1", float, doc = "pt of the standalone muon", precision=14)
    myMuonTable.variables.standaloneEta = Var("? standAloneMuon().isNonnull() ? standAloneMuon().eta() : -99", float, doc = "eta of the standalone muon", precision=14)
    myMuonTable.variables.standalonePhi = Var("? standAloneMuon().isNonnull() ? standAloneMuon().phi() : -99", float, doc = "phi of the standalone muon", precision=14)
    myMuonTable.variables.standaloneCharge = Var("? standAloneMuon().isNonnull() ? standAloneMuon().charge() : -99", float, doc = "phi of the standalone muon", precision=14)    
    
    # Inner Track Algo variables
    myMuonTable.variables.innerTrackAlgo = Var('? innerTrack().isNonnull() ? innerTrack().algo() : -99', 'int', precision=-1, doc='Track algo enum, check DataFormats/TrackReco/interface/TrackBase.h for details.')
    myMuonTable.variables.innerTrackOriginalAlgo = Var('? innerTrack().isNonnull() ? innerTrack().originalAlgo() : -99', 'int', precision=-1, doc='Track original algo enum')

    #Spark Tool Iso 03 variables
    myMuonTable.variables.pfAbsIso03_neu = Var("pfIsolationR03().sumNeutralHadronEt",float,doc="PF absolute isolation dR=0.3, neutral component")
    myMuonTable.variables.pfAbsIso03_pho = Var("pfIsolationR03().sumPhotonEt",float,doc="PF absolute isolation dR=0.3, photon component")
    myMuonTable.variables.pfAbsIso03_sumPU = Var("pfIsolationR03().sumPUPt",float,doc="PF absolute isolation dR=0.3, pu component (no deltaBeta corrections)")
    
    # Spark Tool Iso 04 variables
    myMuonTable.variables.pfAbsIso04_chg = Var("pfIsolationR04().sumChargedHadronPt",float,doc="PF absolute isolation dR=0.4, charged component")
    myMuonTable.variables.pfAbsIso04_neu = Var("pfIsolationR04().sumNeutralHadronEt",float,doc="PF absolute isolation dR=0.4, neutral component")
    myMuonTable.variables.pfAbsIso04_pho = Var("pfIsolationR04().sumPhotonEt",float,doc="PF absolute isolation dR=0.4, photon component")
    myMuonTable.variables.pfAbsIso04_sumPU = Var("pfIsolationR04().sumPUPt",float,doc="PF absolute isolation dR=0.4, pu component (no deltaBeta corrections)")

    #Mini PF Isolation
    myMuonTable.variables.miniPFAbsIso_chg = Var("userFloat('miniIsoChg')",float,doc="mini PF absolute isolation, charged component")
    myMuonTable.variables.miniPFAbsIso_all = Var("userFloat('miniIsoAll')",float,doc="mini PF absolute isolation, total (with scaled rho*EA PU corrections)")
    myMuonTable.variables.miniPFAbsIso_neu = Var("miniPFIsolation().neutralHadronIso()",float,doc="mini PF absolute isolation, neutral component")
    myMuonTable.variables.miniPFAbsIso_pho = Var("miniPFIsolation().photonIso()", float, doc="mini PF absolute isolation, photon component")
    
    # Absolute Isolations for variables already present in Standard NanoAOD as Relative Isolation
    myMuonTable.variables.tkAbsIso = Var("isolationR03().sumPt",float,doc="Tracker-based absolute isolation dR=0.3 for highPt, trkIso",precision=6)
    myMuonTable.variables.pfAbsIso03_chg = Var("pfIsolationR03().sumChargedHadronPt",float,doc="PF absolute isolation dR=0.3, charged component")
    myMuonTable.variables.pfAbsIso03_all = Var("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))",float,doc="PF absolute isolation dR=0.3, total (deltaBeta corrections)")
    myMuonTable.variables.pfAbsIso04_all = Var("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))",float,doc="PF absolute isolation dR=0.4, total (deltaBeta corrections)")
    myMuonTable.variables.jetAbsIso = Var("?userCand('jetForLepJetVar').isNonnull()?(1./userFloat('ptRatio'))-1.:(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))",float,doc="Absolute isolation in matched jet (1/ptRatio-1, pfRelIso04_all if no matched jet)",precision=8)
    #myMuonTable.variables.relTrkiso4 = Var("userFloat('relTrkiso4')",float,doc="Realtive Tracker Iso with cone size 0.4")

    # Muon Quality Variables
    myMuonTable.variables.expectedMatchedStations = Var("expectedNnumberOfMatchedStations()",int,doc="Expected Number of Matched stations")
    myMuonTable.variables.RPCLayers = Var("numberOfMatchedRPCLayers()",int,doc="Number of RPC Layers")
    myMuonTable.variables.stationMask = Var("stationMask()","uint8",doc="Number of masked station")
    myMuonTable.variables.nShowers = Var("numberOfShowers()",int,doc="Number of Showers")
    myMuonTable.variables.muonHits = Var("? globalTrack().isNonnull() ? globalTrack().hitPattern().numberOfValidMuonHits() : ?  innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().numberOfValidMuonHits() :-99",float,doc="Number of valid Muon Hits from either globalTrack or innerTrack")
    ## For completeness I save here also the muonHits for the outer tracker also
    myMuonTable.variables.outerTrackMuonHits = Var("? outerTrack().isNonnull() ? outerTrack().hitPattern().numberOfValidMuonHits() : -99", float, doc = "Number of valid Muon Hits from OuterTrack")
    ##
    myMuonTable.variables.pixelLayers = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().pixelLayersWithMeasurement() : -99", float,doc="Number of Pixel Layers") # No of tracker layers are already saved in the standard NanoAODs
    myMuonTable.variables.validFraction = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().validFraction() : -99", float, doc="Inner Track Valid Fraction")
    myMuonTable.variables.pixelHits = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().numberOfValidPixelHits() : -99", float, doc="Numbr of valid pixel hits")
    myMuonTable.variables.muonStations = Var("? outerTrack().isNonnull() && outerTrack().isAvailable() ? outerTrack().hitPattern().muonStationsWithValidHits() : -99", float, doc="No of valid hits in muon stations")
    myMuonTable.variables.DTHits = Var("? outerTrack().isNonnull() && outerTrack().isAvailable() ? outerTrack().hitPattern().numberOfValidMuonDTHits() : -99", float, doc="No of valid hits in DT")
    myMuonTable.variables.CSCHits = Var("? outerTrack().isNonnull() && outerTrack().isAvailable() ? outerTrack().hitPattern().numberOfValidMuonCSCHits() : -99", float, doc="No of valid hits in CSC")
    myMuonTable.variables.RPCHits = Var("? outerTrack().isNonnull() && outerTrack().isAvailable() ? outerTrack().hitPattern().numberOfValidMuonRPCHits() : -99", float, doc="No of valid hits in RPC")
    
    
    
    # Chi2 related to different tracks
    myMuonTable.variables.trkChi2 = Var("? globalTrack().isNonnull() ? globalTrack().normalizedChi2() : ? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().normalizedChi2() : -99",float,doc="Normalized Chi Square from either globalTrack or innerTrack ")
    myMuonTable.variables.trkChi2_outerTrack = Var("? outerTrack().isNonnull() && outerTrack().isAvailable() ? outerTrack().normalizedChi2() : -99",float,doc="Normalized Chi Square from outerTrack ")
    myMuonTable.variables.trkChi2_innerTrack = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().normalizedChi2() : -99",float,doc="Normalized Chi Square from outerTrack ")
   

    #pt, ptErr, eta, phi, charge for different tracks
    ## ptErr in standard NanoAOD are saved from bestTrack()
    ## For Spark tool it is needed from innerTrack. For completeness outerTrack
    ## variables are also saved
    myMuonTable.variables.innerTrack_ptErr = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().ptError()/innerTrack().pt() : -99", float, doc="InnerTrack Pt Error")
    myMuonTable.variables.outerTrack_ptErr = Var("? outerTrack().isNonnull() && outerTrack().isAvailable() ? outerTrack().ptError()/outerTrack().pt() : -99", float, doc="OuterTrack Pt Error")
    myMuonTable.variables.outerTrack_pt = Var("? outerTrack().isNonnull() && outerTrack().isAvailable() ? outerTrack().pt(): -99", float, doc="OuterTrack Pt")
    myMuonTable.variables.outerTrack_eta = Var("? outerTrack().isNonnull() && outerTrack().isAvailable() ? outerTrack().eta(): -99", float, doc="OuterTrack Eta")
    myMuonTable.variables.outerTrack_phi = Var("? outerTrack().isNonnull() && outerTrack().isAvailable() ? outerTrack().phi(): -99", float, doc="OuterTrack Phi")
    myMuonTable.variables.outerTrack_charge = Var("? outerTrack().isNonnull() && outerTrack().isAvailable() ? outerTrack().charge(): -99", float, doc="OuterTrack charge")
    myMuonTable.variables.innerTrack_charge = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().charge(): -99", float, doc="OuterTrack charge")


    # TuneP related variables
    myMuonTable.variables.tuneP_pt = Var("? tunePMuonBestTrack().isNonnull() ? tunePMuonBestTrack().pt() : -99", float, doc = "pT from tunePMuonBestTrack")
    myMuonTable.variables.tuneP_pterr = Var("? tunePMuonBestTrack().isNonnull() ? tunePMuonBestTrack().ptError() : -99", float, doc = "pTerr from tunePMuonBestTrack")
    myMuonTable.variables.tuneP_muonHits = Var("? tunePMuonBestTrack().isNonnull() ? tunePMuonBestTrack().hitPattern().numberOfValidMuonHits() : -99", int, doc="No of valid muon hists from tunePMuonBestTrack")
   

    #CombinedQuality Variables
    myMuonTable.variables.positionChi2 = Var("combinedQuality().chi2LocalPosition", float, doc="chi2 Local Position")
    myMuonTable.variables.momentumChi2 = Var("combinedQuality().chi2LocalMomentum", float, doc="chi2 Local Momentum")
    myMuonTable.variables.trkKink = Var("combinedQuality().trkKink", float, doc="Track Kink")
    myMuonTable.variables.glbKink = Var("combinedQuality().glbKink", float, doc="Glb Kink")
    myMuonTable.variables.glbTrackProbability = Var("combinedQuality().glbTrackProbability", float, doc="Glb Track Probability")
    myMuonTable.variables.trkRelChi2 = Var("combinedQuality().trkRelChi2",float,doc="Track Rel Chi2")
    
    #timAtIpInOutErr
    myMuonTable.variables.timAtIpInOutErr = Var("time().timeAtIpInOutErr",float,doc="timAtIpInOutErr")
    
    #isArbitratedTracker
    #myMuonTable.variables.isArbitratedTracker = Var("userInt('isArbitratedTracker')", bool, doc = "s Arbitrated Tracker")

    #ExtraidX
    myMuonTable.variables.standaloneExtraIdx = Var('? standAloneMuon().isNonnull() ? standAloneMuon().extra().key() : -99', 'int', precision=-1, doc='Index of the StandAloneTrack TrackExtra in the original collection')
    myMuonTable.variables.innerTrackExtraIdx = Var('? innerTrack().isNonnull() ? innerTrack().extra().key() : -99', 'int', precision=-1, doc='Index of the innerTrack TrackExtra in the original collection')

    #Jet Related Variables
    myMuonTable.variables.jetPtRatio = Var("?userCand('jetForLepJetVar').isNonnull()?min(userFloat('ptRatio'),1.5):1.0/(1.0+(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt)", float, doc="ptRatio using the LepAware JEC approach, for muon MVA")
    myMuonTable.variables.jetDF = Var("?userCand('jetForLepJetVar').isNonnull()?max(userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:probbb')+userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:probb')+userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:problepb'),0.0):0.0",float,doc="b-tagging discriminator of the jet matched to the lepton, for muon MVA")
    myMuonTable.variables.jetCSVv2 = Var("?userCand('jetForLepJetVar').isNonnull()?max(userCand('jetForLepJetVar').bDiscriminator('pfCombinedSecondaryVertexV2BJetTags'),0.0):0.0",float,doc="CSVv2 b-tagging discriminator of the jet matched to the lepton, for muon MVA")
    
    #Dxy Dz variables as of Spark tool
    #myMuonTable.variables.innerTrackDxy = Var("? userInt('isGoodVertex') ? userFloat('innerTrackDxy') : -99.9",float,doc = "dxy from Primary Vertex calculated with Inner Track")
    #myMuonTable.variables.innerTrackDz = Var("? userInt('isGoodVertex') ? userFloat('innerTrackDz') : -99.9",float,doc= "dz from Primary Vertex calculated with Inner Track")

    #nSegments
    #myMuonTable.variables.nsegments = Var("userInt('nsegments')", int, doc = "nsegments as of Spark-tool")

    #Sim Variables
    myMuonTable.variables.simType = Var("? simType() ? simType() : -99",int,doc="simType")
    myMuonTable.variables.simExtType = Var("? simExtType() ? simExtType() : -99",int,doc="simExtType")
    myMuonTable.variables.simFlavour = Var("? simFlavour() ? simFlavour() : -99",int,doc="simFlavour")
    myMuonTable.variables.simHeaviestMotherFlavour = Var(" ? simHeaviestMotherFlavour() ? simHeaviestMotherFlavour() : -99",int,doc="simHeaviestMotherFlavour")
    myMuonTable.variables.simPdgId = Var("? simPdgId() ? simPdgId() : -99",int,doc="simPdgId")
    myMuonTable.variables.simMotherPdgId = Var("? simMotherPdgId() ? simMotherPdgId() : -99",int,doc="simMotherPdgId")
    myMuonTable.variables.simBX = Var("? simBX() ? simBX() : -99",int,doc="simBX")
    myMuonTable.variables.simProdRho = Var("? simProdRho() ? simProdRho(): -99",float,doc="simProdRho")
    myMuonTable.variables.simProdZ = Var("? simProdZ() ? simProdZ(): -99",float,doc="simProdZ")
    myMuonTable.variables.simPt = Var("? simPt() ? simPt(): -99",float,doc="simPt")
    myMuonTable.variables.simEta = Var("? simEta() ? simEta(): -99",float,doc="simEta")
    myMuonTable.variables.simPhi = Var("? simPhi() ? simPhi(): -99",float,doc='simPhi')

    process.globalReplace("muonTable", myMuonTable)


    # PF candidates
    process.isFromTauForPfCand = cms.EDProducer("IsFromPatTauMapProducer",
        packedPFCandidates = cms.InputTag("packedPFCandidates"),
        #patTaus = cms.InputTag("slimmedTaus"),
        patTaus = cms.InputTag("linkedObjects", "taus"),
        #patTaus = cms.InputTag("selectedPatTaus"),
    )
    
    trk_cond = "hasTrackDetails"
    
    # https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#Packed_ParticleFlow_Candidates
    # https://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_12_4_0/doc/html/d8/d79/classpat_1_1PackedCandidate.html
    # lostInnerHits: https://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_12_4_0/doc/html/d8/d79/classpat_1_1PackedCandidate.html#ab9ef9a12f92e02fa61653ba77ee34274
    # fromPV: https://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_12_4_0/doc/html/d8/d79/classpat_1_1PackedCandidate.html#a1e86b4e893738b7cbae410b7f106f339
    process.pfCandTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("packedPFCandidates"),
        #cut = cms.string(""),
        cut = cms.string("pt > 1"),
        name= cms.string("PFCandidate"),
        doc = cms.string("PF candidates"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            CandVars,
            fromPV                  = Var("fromPV"                              , int       , doc = "isolated track comes from PV"),
            lostInnerHits           = Var("lostInnerHits"                       , int       , doc = "Lost inner hits"),
            hasTrackDetails         = Var("hasTrackDetails"                     , bool      , doc = "True if a bestTrack can be extracted from this Candidate"),
            phiAtVtx                = Var("phiAtVtx"                            , float     , doc = "Phi of the candidate's track at the vertex; this is identical to phi() for the vast majority of the particles, but the two might differ for some of them if the calorimeters had contributed significantly in defining the 4-vector of the particle"),
            dxy                     = Var(f"?{trk_cond}?dxy:-999"                , float     , doc = "dxy w.r.t. associated PV"),
            dxyError                = Var(f"?{trk_cond}?dxyError:-999"           , float     , doc = "Error on dxy"),
            dz                      = Var(f"?{trk_cond}?dzAssociatedPV:-999"     , float     , doc = "dz w.r.t. associated PV"),
            dzError                 = Var(f"?{trk_cond}?dzError:-999"            , float     , doc = "Error on dz"),
            vx                      = Var("vx"                                  , float     , doc = "Vertex x"),
            vy                      = Var("vx"                                  , float     , doc = "Vertex y"),
            vz                      = Var("vz"                                  , float     , doc = "Vertex z"),
        ),
        externalVariables = cms.PSet(
            isTauIdxSignalCand     = ExtVar("isFromTauForPfCand:isTauIdxSignalCand"       , int, doc = "Index of the tau if it belongs to pat::Tau::signalCands(); else -1"),
            isTauIdxIsoCand        = ExtVar("isFromTauForPfCand:isTauIdxIsoCand"          , int, doc = "Index of the tau if it belongs to pat::Tau::isolationCands(); else -1"),
            isTauIdxLeadChHadCand  = ExtVar("isFromTauForPfCand:isTauIdxLeadChHadCand"    , int, doc = "Index of the tau if it is pat::Tau::leadChargedHadrCand(); else -1"),
        )
    )
    
    
    # Unfiltered taus
    process.finalTaus.cut = cms.string("pt > 18")
    process.tauTable.doc = cms.string("slimmedTaus after basic selection (" + process.finalTaus.cut.value()+")")
    
    
    # CaloJets
    process.caloJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("slimmedCaloJets"),
        #cut = cms.string(""),
        cut = cms.string("pt > 15"),
        name= cms.string("CaloJet"),
        doc = cms.string("AK4 calo jets"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            P4Vars,
            emEnergyFraction                = Var("emEnergyFraction"                , float),
            #emEnergyInEB                    = Var("emEnergyInEB"                    , float),
            #emEnergyInEE                    = Var("emEnergyInEE"                    , float),
            #emEnergyInHF                    = Var("emEnergyInHF"                    , float),
            energyFractionHadronic          = Var("energyFractionHadronic"          , float),
            #hadEnergyInHB                   = Var("hadEnergyInHB"                   , float),
            #hadEnergyInHE                   = Var("hadEnergyInHE"                   , float),
            #hadEnergyInHF                   = Var("hadEnergyInHF"                   , float),
            #hadEnergyInHO                   = Var("hadEnergyInHO"                   , float),
            #maxEInEmTowers                  = Var("maxEInEmTowers"                  , float),
            #maxEInHadTowers                 = Var("maxEInHadTowers"                 , float),
            #towersArea                      = Var("towersArea"                      , float),
            detectorP4pt                    = Var("detectorP4.Pt"                   , float),
            detectorP4eta                   = Var("detectorP4.Eta"                  , float),
            detectorP4phi                   = Var("detectorP4.Phi"                  , float),
            detectorP4mass                  = Var("detectorP4.M"                    , float),
            detectorP4energy                = Var("detectorP4.E"                    , float),
        ),
    )
    
    if isMC:     
        # GenParticles
        myGenParticleTable = genParticleTable.clone()
        myGenParticleTable.variables.vertexX        = Var("vertex.X"      , float)
        myGenParticleTable.variables.vertexY        = Var("vertex.Y"      , float)
        myGenParticleTable.variables.vertexZ        = Var("vertex.Z"      , float)
        myGenParticleTable.variables.vertexRho      = Var("vertex.Rho"    , float)
        myGenParticleTable.variables.vertexR        = Var("vertex.R"      , float)
    
        process.globalReplace("genParticleTable", myGenParticleTable)
    
    
    ## GenVisTau
    #myGenVisTauTable = genVisTauTable.clone()
    #myGenVisTauTable.variables.vertexX        = Var("vertex.X"      , float)
    #myGenVisTauTable.variables.vertexY        = Var("vertex.Y"      , float)
    #myGenVisTauTable.variables.vertexZ        = Var("vertex.Z"      , float)
    #myGenVisTauTable.variables.vertexRho      = Var("vertex.Rho"    , float)
    #myGenVisTauTable.variables.vertexR        = Var("vertex.R"      , float)
    #
    #process.globalReplace("genVisTauTable", myGenVisTauTable)
    
    if (disTauTagOutputOpt > 0) :
        print("DisTauTagOutputOpt is ", disTauTagOutputOpt)
        process.disTauTag = cms.EDProducer(
            "DisTauTag",
            graphPath = cms.string("data/particlenet_v1_a27159734e304ea4b7f9e0042baa9e22.pb"),
            #jets = cms.InputTag("finalJets"),
            jets = process.jetTable.src,
            pfCandidates = cms.InputTag('packedPFCandidates'),
            save_inputs  = cms.bool(False)
        )
        
        d_disTauTagVars = {
            "disTauTag_score0":     ExtVar("disTauTag:score0"       , float, doc = "Score 0"),
            "disTauTag_score1":     ExtVar("disTauTag:score1"       , float, doc = "Score 1"),
        }
    
    ##process.jetTable.externalVariables = process.jetTable.externalVariables.clone(
    ##    #disTauTag_score0         = ExtVar("disTauTag:score0"       , float, doc = "Score 0"),
    ##    #disTauTag_score1         = ExtVar("disTauTag:score1"       , float, doc = "Score 1"),
    ##    **d_disTauTagVars
    ##)
    
    
    # Create the task
    if (disTauTagOutputOpt == 0) :
         
        process.custom_nanoaod_task = cms.Task(
            process.lostTrackTable,
            process.isFromTauForPfCand,
            process.disMuonTablesTask,
            process.disMuonMCTask,
            process.pfCandTable,
            process.caloJetTable,
        )
    
    elif (disTauTagOutputOpt == 1) :
        
        process.jetTable.externalVariables = process.jetTable.externalVariables.clone(**d_disTauTagVars)
        
        process.custom_nanoaod_task = cms.Task(
            process.lostTrackTable,
            process.isFromTauForPfCand,
            process.pfCandTable,
            process.caloJetTable,
            process.disMuonTablesTask,
            process.disMuonMCTask,
            process.disTauTag,
        )
    
    elif (disTauTagOutputOpt == 2) :
        
        process.jetTable.variables = cms.PSet()
        process.jetTable.externalVariables = cms.PSet(**d_disTauTagVars)
        
        process.custom_nanoaod_task = cms.Task(process.disTauTag)
    
    
    
    # Associate the task to the associate
    process.schedule.associate(process.custom_nanoaod_task)
