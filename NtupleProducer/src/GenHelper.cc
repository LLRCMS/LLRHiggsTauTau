/* 
**
** Helpers for gen info (implementation)
** 
** 
** \date:    13 May 2015
** \author:  L. Cadamuro (LLR)
*/

#include <iostream>
#include "LLRHiggsTauTau/NtupleProducer/interface/GenHelper.h"
#include <DataFormats/Candidate/interface/Candidate.h>

bool genhelper::IsLastCopy (const reco::GenParticle& part)
{
    bool isLast = true;
    int thisPdgId = part.pdgId();

    if (abs(thisPdgId) == 25 || abs(thisPdgId) == 23 || abs(thisPdgId) == 15) // H, Z, tau must decay
        if (part.numberOfDaughters() == 0) return false; // can happen to have a fake "clone" that does not decay --> reject (it is not a real "last")
    
    // other particles, or H/Z/tau with sons
    for (unsigned int iDau = 0; iDau < part.numberOfDaughters(); iDau++)
    {
        const reco::Candidate * Dau = part.daughter(iDau);
        if (Dau->pdgId() == thisPdgId && Dau->numberOfDaughters() > 0) // sometimes a "fake" clone is produced but not decayed
        {
            isLast = false;
            break;
        }
    }
    return isLast;
}

bool genhelper::IsFirstCopy (const reco::GenParticle& part)
{
    bool isFirst = true;
    int thisPdgId = part.pdgId();
    for (unsigned int iMo = 0; iMo < part.numberOfMothers(); iMo++)
    {
        const reco::Candidate * Mo = part.mother(iMo);
        if (Mo->pdgId() == thisPdgId)
        {
            isFirst = false;
            break;
        }
    }
    return isFirst;
}


int genhelper::GetTauDecay (const reco::Candidate* part)
{
    if (abs(part->pdgId()) != 15) return -1; // only on taus
    int decay = -1;
    int nele = 0;
    int nmu = 0;
    for (unsigned int iDau = 0; iDau < part->numberOfDaughters(); iDau++)
    {
        const reco::Candidate * Dau = part->daughter(iDau);
        int dauId = abs(Dau->pdgId());
        if (dauId == 11) nele++;
        if (dauId == 13) nmu++;
    }
    
    if (nmu == 1 && nele == 0) decay = 0;
    if (nmu == 0 && nele == 1) decay = 1;
    if (nmu == 0 && nele == 0) decay = 2;

    return decay; // -1 if strange things happen
}

int genhelper::GetTauDecay (const reco::GenParticle& part)
{
    const reco::Candidate* p = &part;
    return genhelper::GetTauDecay(p);
}

const reco::Candidate* genhelper::GetLastCopy (const reco::Candidate* part)
{
    int cloneInd = -1;
    int id = part->pdgId();
    for (unsigned int iDau = 0; iDau < part->numberOfDaughters(); iDau++)
    {
        const reco::Candidate * Dau = part->daughter( iDau );
        if (id == Dau->pdgId())
        {
            cloneInd = iDau;
            break;
        }
    }
    
    if (cloneInd == -1) return part;
    else return (GetLastCopy (part->daughter(cloneInd)));
    
}

genhelper::HZDecay genhelper::GetHZDecay (const reco::Candidate* part)
{
    int ntau = 0;
    int nele = 0;
    int nmu = 0;
    
    for (unsigned int iDau = 0; iDau < part->numberOfDaughters(); iDau++)
    {
        const reco::Candidate * Dau = part->daughter( iDau );
        if (abs(Dau->pdgId()) == 11 ) nele++;
        if (abs(Dau->pdgId()) == 13 ) nmu++;
        if (abs(Dau->pdgId()) == 15 )
        {
            ntau++;
            int decay = genhelper::GetTauDecay (genhelper::GetLastCopy(Dau));
            if (decay == 0) nmu++;
            if (decay == 1) nele++;
        }
    }
    
    // determine decay mode
    if (ntau == 0)
    {
        if      (nele == 0 && nmu == 2) return genhelper::HZDecay::MuMuPrompt;
        else if (nele == 2 && nmu == 0) return genhelper::HZDecay::EEPrompt;
        else return genhelper::HZDecay::Other;
    }
    
    else if (ntau == 2)
    {
        if (nmu == 0 && nele == 0) return genhelper::HZDecay::HadHad;    
        if (nmu == 0 && nele == 1) return genhelper::HZDecay::EHad;    
        if (nmu == 0 && nele == 2) return genhelper::HZDecay::EE;    
        if (nmu == 1 && nele == 0) return genhelper::HZDecay::MuHad;    
        if (nmu == 1 && nele == 1) return genhelper::HZDecay::EMu;    
        if (nmu == 2 && nele == 0) return genhelper::HZDecay::MuMu;    
    }
    
    return genhelper::HZDecay::Other;
    
}


reco::GenParticle genhelper::GetTauHad (const reco::Candidate* part)
{
    if (abs(part->pdgId()) != 15)
    {
        reco::GenParticle fakeTauH = reco::GenParticle (0, reco::Candidate::LorentzVector(0.,0.,0.,0.), reco::Candidate::Point (0.,0.,0.), -999999, 0, true);
        std::cout << "Warning: building had tau from a particle with pdgId != 15 --> dummy entry returned" << std::endl;
        return fakeTauH;        
    }
    
    reco::Candidate::LorentzVector p4Had (0,0,0,0);
    for (unsigned int iDau = 0; iDau < part->numberOfDaughters(); iDau++)
    {
        const reco::Candidate * Dau = part->daughter( iDau );
        int dauId = abs(Dau->pdgId());
        if (dauId != 12 && dauId != 14 && dauId != 16) // no neutrinos
            p4Had += Dau->p4();
    }
    
    int sign = part->pdgId() / abs(part->pdgId());
    reco::GenParticle TauH = reco::GenParticle (part->charge(), p4Had, part->vertex(), sign*66615, part->status(), true);
    return TauH;
}

const reco::Candidate* genhelper::IsFromID (const reco::Candidate* part, int targetPDGId)
{
    if (abs(part->pdgId()) == targetPDGId) return part;

    for (unsigned int i = 0; i < part->numberOfMothers(); i++)
    {
        const reco::Candidate* matchMoth = genhelper::IsFromID(part->mother(i), targetPDGId);
        if ( matchMoth != NULL) return matchMoth;
    }
    
    // nothing found, mothers finished, exiting...
    return NULL;
    
}

int genhelper::GetIndexInOutput (const reco::Candidate* part, std::vector<const reco::Candidate *> cands)
{
    int index = -1;
    std::vector<const reco::Candidate *>::const_iterator found = find(cands.begin(), cands.end(), part);
    if(found != cands.end()) index = found - cands.begin();
    return index;

}