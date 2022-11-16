# coding: utf-8

"""
Test script using FWLite to interactively explore the content of MiniAOD files.
"""


import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch()


# the LFN to open
# lfn = "/store/mc/RunIISummer20UL17MiniAODv2/GluGluToBulkGravitonToHHTo2B2Tau_M-1000_TuneCP5_PSWeights_narrow_13TeV-madgraph-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v2/30000/679F9DD7-2735-2C4C-8672-5A790CA5812F.root"
lfn = "/store/mc/RunIISummer20UL17MiniAODv2/DYJetsToLL_LHEFilterPtZ-50To100_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v2/2430000/0D401064-953F-0446-9BEA-BBD1E6BDD2E1.root "

# PFN using the IT redirector
pfn = "root://xrootd-cms.infn.it/" + lfn

# handles of objects to open
handles = {
    "jets": {
        "type": "std::vector<pat::Jet>",
        "label": "slimmedJets",
    },
    "fat_jets": {
        "type": "std::vector<pat::Jet>",
        "label": "slimmedJetsAK8",
    },
    "gen_jets": {
        "type": "std::vector<reco::GenJet>",
        "label": "slimmedGenJets",
    },
    "lhe_particles": {
        "type": "LHEEventProduct",
        "label": "externalLHEProducer",
    },
    "gen_particles": {
        "type": "std::vector<reco::GenParticle>",
        "label": "prunedGenParticles",
    },
}


def fwlite_loop(path, handle_data=None, start=0, end=-1, object_type="Event"):
    """
    Opens one or more ROOT files defined by *path* and yields the FWLite event. When *handle_data*
    is not *None*, it is supposed to be a dictionary ``key -> {"type": ..., "label": ...}``. In that
    case, the handle products are yielded as well in a dictionary, mapped to the key, as
    ``(event, objects dict)``.
    """
    ROOT.gSystem.Load("libFWCoreFWLite.so")
    ROOT.gSystem.Load("libDataFormatsFWLite.so")
    ROOT.FWLiteEnabler.enable()
    from DataFormats.FWLite import Events, Runs, Handle  # noqa

    paths = path if isinstance(path, (list, tuple)) else [path]

    handles = {}
    if handle_data:
        for key, data in handle_data.items():
            handles[key] = Handle(data["type"])

    objects = locals()[object_type + "s"](paths)
    if start > 0:
        objects.to(start)

    for i, obj in enumerate(objects):
        if end >= 0 and (start + i) >= end:
            break

        if handle_data:
            products = {}
            for key, data in handle_data.items():
                obj.getByLabel(data["label"], handles[key])
                products[key] = handles[key].product()
            yield obj, products
        else:
            yield obj


# loop through events and start ipython shells for each event (stop with ctrl+c/z)
for i, (event, prod) in enumerate(fwlite_loop(pfn, handles, end=-1)):
    print("\nevent {}".format(i))

    print("Zpt test")
    print("LHE particles")
    lhe = prod["lhe_particles"]
    lhe_particles = lhe.hepeup().PUP
    lhe_ids = lhe.hepeup().IDUP
    lhe_statuses = lhe.hepeup().ISTUP
    lhe_z = ROOT.TLorentzVector()
    for i in range(lhe_particles.size()):
        # debug log
        if abs(lhe_ids[i]) == 23:
            print("LHE Z!", lhe_ids[i], "status", lhe_statuses[i], "vector", lhe_particles[i][3], lhe_particles[i][0], lhe_particles[i][1], lhe_particles[i][2])
        # filter as in https://github.com/cms-sw/cmssw/blob/master/GeneratorInterface/GenFilters/plugins/LHEVpTFilter.cc
        # must be outgoing particle
        if lhe_statuses[i] != 1:
            continue
        # must be lepton
        if not (11 <= abs(lhe_ids[i]) <= 16):
            continue
        print("id", lhe_ids[i], "status", lhe_statuses[i], "vector", lhe_particles[i][3], lhe_particles[i][0], lhe_particles[i][1], lhe_particles[i][2])
        lhe_z += ROOT.TLorentzVector(lhe_particles[i][0], lhe_particles[i][1], lhe_particles[i][2], lhe_particles[i][3])
    print("LHE Z sum", lhe_z.Energy(), lhe_z.Pt(), lhe_z.M())

    print("\nGen particles")
    gen = prod["gen_particles"]
    gen_z = ROOT.TLorentzVector()
    for i in range(gen.size()):
        p = gen[i]
        # must be from hard process
        if not p.statusFlags().fromHardProcess():
            continue
        # must be lepton
        if not (11 <= abs(p.pdgId()) <= 16):
            continue
        print("id", p.pdgId(), "status", p.status(), "vector", p.energy(), p.px(), p.py(), p.pz())
        gen_z += ROOT.TLorentzVector(p.px(), p.py(), p.pz(), p.energy())
    print("Gen Z sum", gen_z.Energy(), gen_z.Pt(), gen_z.M())

    print("\nZ Gen0")
    n_gen_z = 0
    for i in range(gen.size()):
        p = gen[i]
        # must be Z boson
        if abs(p.pdgId()) != 23:
            continue
        print("Z", "status", p.status(), "vector", p.energy(), p.pt(), p.mass())
        n_gen_z += 1
    print("Zpt test end")

    from IPython import embed; embed()
