//C++
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
//ROOT
#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <cmath>
//#include <Tobject.h>
//#include <TBranch.h>
//#include<Isotropy.C>
#if !defined(__CINT__) || defined(__MAKECINT__)
//WCSim
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
//BONSAI
#include "WCSimBonsai.hh"
#endif


float DWALL(float x, float y, float r)
{

    float posFromCent = sqrt(pow(x, 2) + pow(y, 2));
    float dwall = r - posFromCent;
    return dwall;
}

int main() {

    bool neut = true;
    std::string filename = "FullNeutronSimulation.root";

    //Add libraries here

    WCSimBonsai* bonsai = new WCSimBonsai();

    TFile *file = new TFile(filename.c_str(),"read");
    TTree *tree = (TTree*)file->Get("wcsimT");
    int nevent = tree->GetEntries();
    WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
    TBranch *branch = tree->GetBranch("wcsimrootevent");
    branch->SetAddress(&wcsimrootsuperevent);
    tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

    TTree *geotree = (TTree*)file->Get("wcsimGeoT");
    WCSimRootGeom *geo = 0;
    geotree->SetBranchAddress("wcsimrootgeom", &geo);
    geotree->GetEntry(0);
    bonsai->Init(geo);

    WCSimRootTrigger* wcsimrootevent;

    for (int i =1000; i<1002;i++)
    {


        TString outName = "Bonsai_" + std::to_string(i) + filename;
        TFile *output = new TFile(outName,"RECREATE");


        int numCherDigi = 0;
        float dWall = 0;
        float recoX = 0;
        float recoY = 0;
        float recoZ = 0;
        float trueX = 0;
        float trueY = 0;
        float trueZ = 0;
        float recoE = 0;
        float trueE = 0;
        float recoR = 0;
        float recoT = 0;
        float recoTheta = 0;
        float recoPhi = 0;
        float recoAlpha = 0;
        float recoCone = 0;
        float recoGoodness = 0;
        float recoEllipticity = 0;
        float recoLikelihood = 0;
        float trueT = 0;
        int ncaps = 0;
        int ngamma = 0;
        int nucleus = 0;



        TTree *deets = new TTree("Data","Data");
        deets->Branch("numCherDigi",&numCherDigi, "numCherDigi/I");
        deets->Branch("dWall",&dWall, "dWall/F");
        deets->Branch("recoX",&recoX, "recoX/F");
        deets->Branch("recoY",&recoY, "recoY/F");
        deets->Branch("recoZ",&recoZ, "recoZ/F");
        deets->Branch("trueX",&trueX, "trueX/F");
        deets->Branch("trueY",&trueY, "trueY/F");
        deets->Branch("trueZ",&trueZ, "trueZ/F");
        deets->Branch("trueE",&trueE, "trueE/F");
        deets->Branch("recoE",&recoE, "recoE/F");
        deets->Branch("trueT",&trueT, "trueT/F");
        deets->Branch("recoT",&recoT, "recoT/F");
        deets->Branch("recoR",&recoR, "recoR/F");
        deets->Branch("recoTheta",&recoTheta, "recoTheta/F");
        deets->Branch("recoPhi",&recoPhi, "recoPhi/F");
        deets->Branch("recoAlpha",&recoAlpha, "recoAlpha/F");
        deets->Branch("recoCone",&recoCone, "recoCone/F");
        deets->Branch("recoGoodness",&recoGoodness, "recoGoodness/F");
        deets->Branch("recoEllipticity",&recoEllipticity, "recoEllipticity/F");
        deets->Branch("recoLikelihood",&recoLikelihood, "recoLikelihood/F");


        if (neut)
        {
            deets->Branch("nucleus",&nucleus, "nucleus/I");
            deets->Branch("ngamma",&ngamma, "ngamma/I");

        }

        int startJ = (i-1000)*1000;
        for (int j = startJ;j<2000;j++)
        {

            if (j % 1000 == 0 && j != startJ) {break;} // Check to make sure the event
            std::cout <<"Event : " << j<<std::endl;

            float bsT[500],bsQ[500];
            float bsvertex[4],bsresult[6];
            float bsgood[500];
            int bsCAB[500];
            int bsnhit[1];
            int bsnsel[2];

            tree->GetEntry(j);
//                      wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);

            for (int index = 0; index < wcsimrootsuperevent->GetNumberOfEvents(); index++)
            {
                wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);

                        numCherDigi = wcsimrootevent->GetNcherenkovdigihits();

                if (neut)
                {
                    ncaps = wcsimrootevent->GetNcaptures();

                    if (ncaps != 0)
                    {
                        TObject *captureN = (wcsimrootevent->GetCaptures())->At(0);
                        WCSimRootCapture *capInfo = dynamic_cast<WCSimRootCapture*>(captureN);

                        ngamma = capInfo->GetNGamma();
                        trueT = capInfo->GetCaptureT();
                        trueX = capInfo->GetCaptureVtx(0);
                        trueY = capInfo->GetCaptureVtx(1);
                        trueZ = capInfo->GetCaptureVtx(2);
                        nucleus = capInfo->GetCaptureNucleus();


                    }

                }
                else
                {
                    trueX = wcsimrootevent->GetVtx(0);
                    trueY = wcsimrootevent->GetVtx(1);
                    trueZ = wcsimrootevent->GetVtx(2);
                    trueT = 0;
                }



                bsnhit[0] = numCherDigi;
//                std::cout<<"Num Hits: "<<numCherDigi<<std::endl;
                for (int hitNum = 0; hitNum < numCherDigi; hitNum++ )
                {

                    TObject *element = (wcsimrootevent->GetCherenkovDigiHits())->At(hitNum);

                    WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit =
                            dynamic_cast<WCSimRootCherenkovDigiHit*>(element);
//                        std::cout<< "HitNum" << hitNum<<"Hit Time: "<< wcsimrootcherenkovdigihit->GetT()<<std::endl;
                    bsT[hitNum]=wcsimrootcherenkovdigihit->GetT();
                    bsQ[hitNum]=wcsimrootcherenkovdigihit->GetQ();
                    bsCAB[hitNum]=wcsimrootcherenkovdigihit->GetTubeId();

                }//End loop over digitized hits

                bool bonsaiStatus = true;
                if (bonsai->BonsaiFit( bsvertex, bsresult, bsgood, bsnsel, bsnhit, bsCAB, bsT, bsQ) ==0)
                {
                    bonsaiStatus = false;
                }

                bool savedata = false;

                if (neut && ncaps != 0 ) savedata = true;
                if (!neut) savedata = true;

                if (savedata && bonsaiStatus)
                {
                    recoX = bsvertex[0] ;
                    recoY = bsvertex[1] ;
                    recoZ = bsvertex[2] ;
                    recoT = bsvertex[3] ;
                    recoR = sqrt(pow(recoX-trueX, 2) + pow(recoY-trueY, 2) + pow(recoZ-trueZ, 2));
                    dWall = DWALL(recoX,recoY,400.);
                    recoTheta = bsresult[0];
                    recoPhi = bsresult[1];
                    recoAlpha = bsresult[2];
                    recoCone = bsresult[3];
                    recoEllipticity = bsresult[4];
                    recoLikelihood = bsresult[5];
                    recoGoodness = bsgood[2];
//                    std::cout<<"x: "<<trueX <<" y: "<<trueY <<" z: "<<trueZ <<" rx: "<<recoX <<" ry: "<<recoY <<" rz: "<<recoZ<<std::endl;
                    deets->Fill();
                }


            }//End loop over subevents

        }//loop stops every 1000 events

        output->cd();
        deets->Write();
        output->Close();

    }// Ends the loop of events

    return 0;
}