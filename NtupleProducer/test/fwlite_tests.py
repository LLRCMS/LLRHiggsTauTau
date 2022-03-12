# coding: utf-8

"""
Test script using FWLite to interactively explore the content of MiniAOD files.
"""


import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch()


# the LFN to open
lfn = "/store/mc/RunIISummer20UL17MiniAODv2/GluGluToBulkGravitonToHHTo2B2Tau_M-1000_TuneCP5_PSWeights_narrow_13TeV-madgraph-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v2/30000/679F9DD7-2735-2C4C-8672-5A790CA5812F.root"

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
for i, (event, prod) in enumerate(fwlite_loop(pfn, handles, end=10)):
    print("event {}".format(i))
    from IPython import embed; embed()
