/*
** author: L.Cadamuro (LLR)
** date:   11/06/2015
** 
** Class to manage multiple related histos.
** All histos are referenced with a std::map and retrieved using a string tag
** An unique class string tag is also prepended to make all histos unique even when more
** HistoManager objects are created
**
** For more flexibility the map stores pointers to TObjects even if in general only TH1
** will be stored.
*/

#ifndef HISTOMANAGER_H
#define HISTOMANAGER_H

#include <iostream>
#include <map>

#include "TObject.h"
#include "TH1D.h"
#include "TFile.h"

using namespace std;
typedef std::map<std::string, TObject*>::iterator it_type;

class HistoManager
{
    public:
        HistoManager (const char* tag);
        ~HistoManager();
        int AddElement (TObject* ptr, const char* objTag); // 0: added, -1: not added
        TObject* GetElement (const char* objTag);
        void AddNewHisto (const char* name, const char* title, int nbinsx, double xlow, double xup); // creates a new histo
        TH1D* GetHisto(const char* name);
        void SaveAllToFile (TFile* fOut);
        
    private:
        std::string _tag;
        std::map <std::string, TObject*> _map;
        string MakeStoredName(const char * objTag); // adds internal tag to name for the search
        
};

HistoManager::HistoManager (const char* tag) : _tag (tag) {}

HistoManager::~HistoManager()
{
    // delete all allocated objs
    // for this, always store pointers to studd allocated on the heap!!
    for (it_type it = _map.begin(); it != _map.end(); ++it)
    {
        delete (it->second);
    }
}

int HistoManager::AddElement (TObject* ptr, const char* objTag)
{
    string name (MakeStoredName(objTag));
    
    // not optimal, but I hope to have not so many histos!
    if (_map.find (name) != _map.end())
    {
        cout << "TObject with name " << objTag << " already exists, no object added" << endl;
        return -1;
    }
    
    else
    {
        _map[name] = ptr;
        return 0;
    }
}

TObject* HistoManager::GetElement (const char * objTag)
{
    it_type it = _map.find (MakeStoredName(objTag));
    if (it != _map.end())
        return it->second;
    else
    {
        cout << "No TObject found with name " << objTag << endl;
        return NULL;
    }
}

void HistoManager::AddNewHisto (const char* name, const char* title, int nbinsx, double xlow, double xup)
{
    string fullName = MakeStoredName (name);
    TH1D* h = new TH1D (fullName.c_str(), title, nbinsx, xlow, xup);
    int i = AddElement (h, name);
    if (i == -1) delete h; // if not added
}

TH1D* HistoManager::GetHisto(const char* name)
{
    return ((TH1D*) GetElement(name));
}

string HistoManager::MakeStoredName(const char * objTag)
{
    string name (_tag);
    name += "_";
    name += objTag;
    return name;

}

void HistoManager::SaveAllToFile (TFile* fOut)
{
    fOut->cd();
    for (it_type it = _map.begin(); it != _map.end(); ++it)
    {
        it->second->Write();
    }
}
#endif
