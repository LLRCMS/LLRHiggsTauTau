TAG_GENERALE="prod24Lu2015"

python submitOnTier3_LLRHTauTau.py -s Data_DoubleMuon.py -w True -m -1 -t "Data_DoubleMuon_$TAG_GENERALE" -n 40 -i False
python submitOnTier3_LLRHTauTau.py -s Data_DoubleEG.py -w True -m -1 -t "Data_DoubleEG_$TAG_GENERALE" -n 40 -i False
python submitOnTier3_LLRHTauTau.py -s Data_SingleElectron.py -w True -m -1 -t "Data_SingleElectron_$TAG_GENERALE" -n 40 -i False
python submitOnTier3_LLRHTauTau.py -s Data_SingleMu.py -w True -m -1 -t "Data_SinglMu_$TAG_GENERALE" -n 40 -i False
python submitOnTier3_LLRHTauTau.py -s Data_SingleMuon.py -w True -m -1 -t "Data_SingleMuon_$TAG_GENERALE" -n 40 -i False
python submitOnTier3_LLRHTauTau.py -s Data_Tau.py -w True -m -1 -t "Data_Tau_$TAG_GENERALE" -n 40 -i False