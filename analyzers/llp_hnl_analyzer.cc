#include "llp_hnl_analyzer.h"
#include "RazorHelper.h"
#include "HNLMuonSystemTree.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
#include "DBSCAN.h"
//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"


using namespace std::chrono;
using namespace std;

struct greater_than_pt
{
  inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){return p1.Pt() > p2.Pt();}
};

struct leptons
{
  TLorentzVector lepton;
  int pdgId;
  float dZ;
  // bool passLooseId;
  // bool passMediumId;
  bool passId;
  bool passVetoId;
  bool passLooseIso;
  bool passTightIso;
  bool passVTightIso;
  bool passVVTightIso;

};


struct jets
{
  TLorentzVector jet;
  float time;
  bool passId;
  // bool passLooseId;
  // bool passMediumId;
  // bool passTightId;
  bool isCSVL;
  int ecalNRechits;
  float ecalRechitE;
  float jetChargedEMEnergyFraction;
  float jetNeutralEMEnergyFraction;
  float jetChargedHadronEnergyFraction;
  float jetNeutralHadronEnergyFraction;
  bool jetPassMuFrac;
  float jetPtJESUp;
  float jetPtJESDown;
  float jetEJESUp;
  float jetEJESDown;
  float JecUnc;

  float electronEnergyFraction;
  float neutralEmEnergyFraction;
  float chargedHadronEnergyFraction;
  float neutralHadronEnergyFraction;
  float muonEnergyFraction;

};

//lepton highest pt comparator
struct largest_pt
{
  inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
} my_largest_pt;

//jet highest pt comparator
struct largest_pt_jet
{
  inline bool operator() (const jets& p1, const jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
} my_largest_pt_jet;

void llp_hnl_analyzer::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  //options format: MH/MX/ctau/condor: 1000/300/0/1
  // mh can be 3-4 digits, mx is always 3 digits, ctau is one digit(number of zeros), last digit is condor option
  // mh can be 3-4 digits, mx is always 3 digits, ctau is 2 digit(number of zeros), last digit is condor option
  //
  //
  // int mx = int(options/1000)%1000;
  // int mh = options/1000000;
  // int ctau = pow(10, int(options/10)%10) * int(int(options/100)%10);
  //
  // cout<<"mh "<<mh<<", mx "<<mx<<", ctau "<<ctau<<endl;

  bool signalScan = int(options/10) == 1;
  int option = options%10;
  // if (options % 1){
  //   option = 1; // used when running condor
  // }
  // else{
  //   option = 0;// used when running locally
  // }

  if( isData )
    {
      std::cout << "[INFO]: running on data with option: " << option << std::endl;
    }
  else
    {
      std::cout << "[INFO]: running on MC with option: " << option << std::endl;
    }
  if( signalScan )
    {
      std::cout << "[INFO]: running with Signal scan" << std::endl;
    }
  else
    {
      std::cout << "[INFO]: running without Signal scan " << option << std::endl;
    }



  const float ELE_MASS = 0.000511;
  const float MU_MASS  = 0.105658;
  const float Z_MASS   = 91.2;

  if (analysisTag == ""){
    analysisTag = "Razor2016_80X";

  }


  //-----------------------------------------------
  //Set up Output File
  //-----------------------------------------------
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "HeavyNeutralLepton_Tree.root";
  TFile *outFile;
  //if (isData || !signalScan) outFile = new TFile(outfilename.c_str(), "RECREATE");
  outFile = new TFile(outfilename.c_str(), "RECREATE");


  HNLMuonSystemTree *MuonSystem = new HNLMuonSystemTree;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();

  // for signals, need one output file for each signal point
  map<pair<int,int>, TFile*> Files2D;
  map<pair<int,int>, TTree*>Trees2D;
  map<pair<int,int>, TH1F*> NEvents2D;
  map<pair<int,int>, TH1F*> accep2D;
  map<pair<int,int>, TH1F*> accep_met2D;
  map<pair<int,int>, TH1F*> Total2D;


  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *Total = new TH1F("Total", "Total", 1, 1, 2);

  TH1F *accep = new TH1F("accep", "acceptance", 1, 1, 2);
  TH1F *accep_met = new TH1F("accep_met", "acceptance_met", 1, 1, 2);

  TH1F *Nmet200 = new TH1F("Nmet200", "Nmet200", 1, 1, 2);
  TH1F *NmetFilter = new TH1F("NmetFilter", "NmetFilter", 1, 1, 2);
  TH1F *Nlep0 = new TH1F("Nlep0", "Nlep0", 1, 1, 2);
  TH1F *Njet1 = new TH1F("Njet1", "Njet1", 1, 1, 2);
  TH1F *NcosmicVeto = new TH1F("NcosmicVeto", "NcosmicVeto", 1, 1, 2);

  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
  if(cmsswPath != NULL and option == 1) pathname = "JEC/"; //run on condor if option == 1

  //--------------------------------
  //Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  helper = new RazorHelper(analysisTag, isData, false);

  std::vector<FactorizedJetCorrector*> JetCorrector = helper->getJetCorrector();
  std::vector<std::pair<int,int> > JetCorrectorIOV = helper->getJetCorrectionsIOV();

  //*************************************************************************
  //Loop over Input File Events
  //*************************************************************************
  if (fChain == 0) return;
  cout << "Total Events: " << fChain->GetEntries() << "\n";
  Long64_t nbytes = 0, nb = 0;
  clock_t start, end;
  start = clock();

  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {
    //begin event
    if(jentry % 10000 == 0) {
      end = clock();
      double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
      cout << "Processing entry " << jentry+1 << endl;
      cout << "Time taken by program is : " << time_taken << endl;
      start = clock();
    }
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    // cout << "Load event" << endl;
    //GetEntry(ientry);
    nb = fChain->GetEntry(jentry); nbytes += nb;

    //std::cout << "deb0 " << jentry << std::endl;
    MuonSystem->InitVariables();
    //std::cout << "deb1 " << jentry << std::endl;

    // cout<<*lheComments<<endl;
    if (!isData && signalScan) {

      string mh_substring = lheComments->substr(lheComments->find("MH-")+3);
      int mh = stoi(mh_substring.substr(0,mh_substring.find('_')));
      string mx_substring = lheComments->substr(lheComments->find("MS-")+3);
      int mx = stoi(mx_substring.substr(0,mx_substring.find('_')));
      string ctau_substring = lheComments->substr(lheComments->find("ctauS-")+6);
      int ctau = stoi(ctau_substring.substr(0,ctau_substring.find('_')));
      MuonSystem->mH = mh;
      MuonSystem->mX = mx;
      MuonSystem->ctau = ctau;

      // if (mh2 != mh || mx2!=mx || ctau2!=ctau) continue;
      // cout<<*lheComments<<endl;

      pair<int,int> signalPair = make_pair(mx, ctau);

      if (Files2D.count(signalPair) == 0){ //create file and tree
        //format file name
        string thisFileName = outfilename;
        thisFileName.erase(thisFileName.end()-5, thisFileName.end());
        thisFileName += "_" + to_string(mx) + "_" + to_string(ctau) + ".root";

        Files2D[signalPair] = new TFile(thisFileName.c_str(), "recreate");
        Trees2D[signalPair] =  MuonSystem->tree_->CloneTree(0);
        NEvents2D[signalPair] = new TH1F(Form("NEvents%d%d", mx, ctau), "NEvents", 1,0.5,1.5);
        Total2D[signalPair] = new TH1F(Form("Total%d%d", mx, ctau), "Total", 1,0.5,1.5);
        accep2D[signalPair] = new TH1F(Form("accep2D%d%d", mx, ctau), "acceptance", 1,0.5,1.5);
        accep_met2D[signalPair] = new TH1F(Form("accep_met2D%d%d", mx, ctau), "acceptance_met", 1,0.5,1.5);
      }
      //Fill NEvents hist
      // cout<<"here"<<endl;
      NEvents2D[signalPair]->Fill(1.0, genWeight);
    }

    // cout << "Get event level info" << endl;

    //event info
    if (isData) {
      NEvents->Fill(1);
      MuonSystem->weight = 1;
    }
    else {
      // cout<<*lheComments<<endl;
      MuonSystem->weight = genWeight;
      NEvents->Fill(1, genWeight);
      // NEvents_genweight->Fill(1);
    }
    MuonSystem->runNum = runNum;
    MuonSystem->lumiSec = lumiNum;
    MuonSystem->evtNum = eventNum;

    //cout << "Get Gen level info" << endl;

    bool wzFlag = false;
    if (!isData) {
      for (int i=0; i < nGenParticle; i++) {

        if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == 24) {
          // choosing only the W->munu events
          wzFlag = true;
          MuonSystem->gLepId = gParticleId[i];
          MuonSystem->gLepPt = gParticlePt[i];
          MuonSystem->gLepEta = gParticleEta[i];
          MuonSystem->gLepE = gParticleE[i];
          MuonSystem->gLepPhi = gParticlePhi[i];
        }
        else if (abs(gParticleId[i]) == 15 && gParticleStatus[i] == 2 && abs(gParticleMotherId[i]) == 24){
          wzFlag = true;
          MuonSystem->gLepId = gParticleId[i];
          MuonSystem->gLepPt = gParticlePt[i];
          MuonSystem->gLepEta = gParticleEta[i];
          MuonSystem->gLepE = gParticleE[i];
          MuonSystem->gLepPhi = gParticlePhi[i];
        }
        if (abs(gParticleId[i]) == 24 && gParticleStatus[i]==62) {
          MuonSystem->gWPt = gParticlePt[i];
        }
        if (abs(gParticleId[i])== 25 || abs(gParticleId[i] == 35)) {
          MuonSystem->gHiggsPt = gParticlePt[i];
          MuonSystem->gHiggsEta = gParticleEta[i];
          MuonSystem->gHiggsPhi = gParticlePhi[i];
          MuonSystem->gHiggsE = gParticleE[i];
        }
        if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == 23) {
          //  Z->mumu/Z->ee
          MuonSystem->ZCategory  = 0;
        }
        else if (abs(gParticleId[i]) == 15 && gParticleStatus[i] == 2 && abs(gParticleMotherId[i]) == 23){
          //  Z->tautau
          MuonSystem->ZCategory  = 0;
        }
        else if ((abs(gParticleId[i]) == 12 || abs(gParticleId[i]) == 14 || abs(gParticleId[i]) == 16) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == 23){
          //Z->nunu
          MuonSystem->ZCategory  = 1;
        }
        else if ((abs(gParticleId[i]) < 6) && gParticleStatus[i] == 23 && abs(gParticleMotherId[i]) == 23){
          //Z->qq
          MuonSystem->ZCategory  = 2;
        }
         //MuonSystem->gParticleStatus[MuonSystem->nGenParticle] = gParticleStatus[i];
         //MuonSystem->gParticleId[MuonSystem->nGenParticle]  = gParticleId[i];
         //MuonSystem->gParticleMotherId[MuonSystem->nGenParticle]  = gParticleMotherId[i];
         //MuonSystem->gParticleMotherIndex[MuonSystem->nGenParticle]  = gParticleMotherIndex[i];
         //MuonSystem->gParticlePt[MuonSystem->nGenParticle]  = gParticlePt[i];
         //MuonSystem->gParticleEta[MuonSystem->nGenParticle]  = gParticleEta[i];
         //MuonSystem->gParticlePhi[MuonSystem->nGenParticle]  = gParticlePhi[i];
         //MuonSystem->gParticleE[MuonSystem->nGenParticle]  = gParticleE[i];
         //MuonSystem->gParticleProdVertexX[MuonSystem->nGenParticle]  = gParticleProdVertexX[i];
         //MuonSystem->gParticleProdVertexY[MuonSystem->nGenParticle]  = gParticleProdVertexY[i];
         //MuonSystem->gParticleProdVertexZ[MuonSystem->nGenParticle]  = gParticleProdVertexZ[i];
         //MuonSystem->nGenParticle++;
         //cout<<"genparticles: "<<MuonSystem->nGenParticle<<endl;
      }
      // MuonSystem->nGenJets = nGenJets;
      // for(int i=0; i < nGenJets; i++) {
      //   MuonSystem->genJetE[i]  = genJetE[i];
      //   MuonSystem->genJetPt[i]  = genJetPt[i];
      //   MuonSystem->genJetEta[i]  = genJetEta[i];
      //   MuonSystem->genJetPhi[i]  = genJetPhi[i];
      // }



      //cout << "Get Gen jet level info" << endl;
      MuonSystem->genMetPtTrue = genMetPtTrue;
      MuonSystem->genMetPhiTrue = genMetPhiTrue;
      MuonSystem->genMetPtCalo = genMetPtCalo;
      MuonSystem->genMetPhiCalo = genMetPhiCalo;

      for(int i = 0; i < 1;i++) {
        MuonSystem->gLLP_eta[i] = gLLP_eta[i];
        MuonSystem->gLLP_phi[i] = gLLP_phi[i];
        MuonSystem->gLLP_decay_vertex_r[i] = sqrt(gLLP_decay_vertex_x[i]*gLLP_decay_vertex_x[i]+gLLP_decay_vertex_y[i]*gLLP_decay_vertex_y[i]);
        MuonSystem->gLLP_decay_vertex_x[i] = gLLP_decay_vertex_x[i];
        MuonSystem->gLLP_decay_vertex_y[i] = gLLP_decay_vertex_y[i];
        MuonSystem->gLLP_decay_vertex_z[i] = gLLP_decay_vertex_z[i];
        float beta = gLLP_beta[i];
        float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[i], 2) + pow(MuonSystem->gLLP_decay_vertex_z[i],2));
        float gamma = 1.0/sqrt(1-beta*beta);
        MuonSystem->gLLP_ctau[i] = gLLP_decay_vertex/(beta * gamma);
        MuonSystem->gLLP_beta[i] = gLLP_beta[i];
        MuonSystem->gLLP_e[i] = gLLP_e[i];
        MuonSystem->gLLP_pt[i] = gLLP_pt[i];
        MuonSystem->gLLP_lepdPhi[i] = deltaPhi(MuonSystem->gLepPhi,MuonSystem->gLLP_phi[i]);

        // endcap
        if (abs(MuonSystem->gLLP_eta[i]) < 2.4 && abs(MuonSystem->gLLP_eta[i]) > 0.9
            && abs(MuonSystem->gLLP_decay_vertex_z[i])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[i])>568
            && MuonSystem->gLLP_decay_vertex_r[i] < 695.5) {
          MuonSystem->gLLP_csc[i] = true;
        }
        if (abs(MuonSystem->gLLP_decay_vertex_z[i])< 661.0
          && MuonSystem->gLLP_decay_vertex_r[i] < 738.0
           && MuonSystem->gLLP_decay_vertex_r[i] > 380.0) MuonSystem->gLLP_dt[i] = true;

      }//end of nLLP loop
      // Fill daughter info
      for(int i = 0; i < 3;i++)
      {
       MuonSystem->gLLP_daughter_id[i] = gLLP_daughter_id[i];
       MuonSystem->gLLP_daughter_pt[i] = gLLP_daughter_pt[i];
       MuonSystem->gLLP_daughter_eta[i] = gLLP_daughter_eta[i];
       MuonSystem->gLLP_daughter_phi[i] = gLLP_daughter_phi[i];
       MuonSystem->gLLP_daughter_e[i] = gLLP_daughter_e[i];
       MuonSystem->gLLP_daughter_mass[i] = gLLP_daughter_mass[i];


      }




      for (int i=0; i < nBunchXing; i++) {
        if (BunchXing[i] == 0) {
          MuonSystem->npu = nPUmean[i];
        }
      }
      MuonSystem->pileupWeight = helper->getPileupWeight(MuonSystem->npu);
      MuonSystem->pileupWeightUp = helper->getPileupWeightUp(MuonSystem->npu) / MuonSystem->pileupWeight;
      MuonSystem->pileupWeightDown = helper->getPileupWeightDown(MuonSystem->npu) / MuonSystem->pileupWeight;
    }//end MConly if loop


    //cout << "Get PU level info" << endl;

    //get NPU
    MuonSystem->npv = nPV;
    MuonSystem->rho = fixedGridRhoFastjetAll;
    MuonSystem->met = metType1Pt;
    MuonSystem->metPhi = metType1Phi;
    MuonSystem->metJESUp = MuonSystem->met;
    MuonSystem->metJESDown = MuonSystem->met;

    //cout << "Extract MET" << endl;
    MuonSystem->metSF = helper->getMetTriggerSF(MuonSystem->met);

    if(signalScan && !isData)Total2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->pileupWeight);
    // comment this out because we want LLPs either in barrel or in endcap
    // need branch gLLP_dt
    // if(!isData)
    //   {
    //     if (MuonSystem->gLLP_csc[0] == false && MuonSystem->gLLP_csc[1] == false)continue;
    //   }

    if(signalScan && !isData)accep2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->pileupWeight);
    else if (!isData) accep->Fill(1.0, genWeight*MuonSystem->pileupWeight);

    //cout << "Corrected MET" << endl;
    std::pair<double,double> corrected_met;
    if (analysisTag=="Razor2016_07Aug2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2016, !isData, nPV);
    else if (analysisTag=="Razor2017_17Nov2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2017, !isData, nPV);
    else if (analysisTag=="Razor2018_17SeptEarlyReReco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2018, !isData, nPV);

    MuonSystem->metXYCorr = corrected_met.first;
    MuonSystem->metPhiXYCorr = corrected_met.second;

    MuonSystem->metEENoise = MuonSystem->metXYCorr;
    MuonSystem->metHEM = MuonSystem->metXYCorr;
    MuonSystem->metPhiEENoise = MuonSystem->metPhiXYCorr;
    MuonSystem->metPhiHEM = MuonSystem->metPhiXYCorr;

    //Triggers (trigger_names_llp_v3.dat!)
    for(int i = 0; i < NTriggersMAX; i++){
      MuonSystem->HLTDecision[i] = HLTDecision[i];
    }
    MuonSystem->SingleMuonTrigger = (HLTDecision[135] or
                                     HLTDecision[136] or
                                     HLTDecision[196]);;
    MuonSystem->SingleEleTrigger = (HLTDecision[87] or
                                    HLTDecision[521] or
                                    HLTDecision[516] or
                                    HLTDecision[867] or
                                    HLTDecision[868]);
    MuonSystem->SingleLepTrigger = (MuonSystem->SingleMuonTrigger or
                                    MuonSystem->SingleEleTrigger);

    //cout << "Check if trigger by single lepton " << MuonSystem->SingleLepTrigger << endl;

    // flags
    MuonSystem->Flag_HBHENoiseFilter = Flag_HBHENoiseFilter;
    MuonSystem->Flag_HBHEIsoNoiseFilter = Flag_HBHEIsoNoiseFilter;
    MuonSystem->Flag_BadPFMuonFilter = Flag_BadPFMuonFilter;
    // MuonSystem->Flag_CSCTightHaloFilter = Flag_CSCTightHaloFilter;
    MuonSystem->Flag_globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016Filter;
    MuonSystem->Flag_goodVertices = Flag_goodVertices;
    MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;
    MuonSystem->Flag_BadChargedCandidateFilter = Flag_BadChargedCandidateFilter;
    MuonSystem->Flag_eeBadScFilter = Flag_eeBadScFilter;
    // MuonSystem->Flag_all = Flag_goodVertices && Flag_HBHEIsoNoiseFilter && Flag_BadPFMuonFilter && Flag_CSCTightHaloFilter && Flag_goodVertices && Flag_ecalBadCalibFilter;

    MuonSystem->Flag2_HBHENoiseFilter = Flag2_HBHENoiseFilter;
    MuonSystem->Flag2_HBHEIsoNoiseFilter = Flag2_HBHEIsoNoiseFilter;
    MuonSystem->Flag2_BadPFMuonFilter = Flag2_BadPFMuonFilter;
    MuonSystem->Flag2_globalSuperTightHalo2016Filter = Flag2_globalSuperTightHalo2016Filter;
    MuonSystem->Flag2_globalTightHalo2016Filter = Flag2_globalTightHalo2016Filter;
    MuonSystem->Flag2_BadChargedCandidateFilter =Flag2_BadChargedCandidateFilter;
    // Flag2_goodVertices = Flag2_goodVertices;
    MuonSystem->Flag2_EcalDeadCellTriggerPrimitiveFilter = Flag2_EcalDeadCellTriggerPrimitiveFilter;
    MuonSystem->Flag2_ecalBadCalibFilter = Flag2_ecalBadCalibFilter;
    MuonSystem->Flag2_eeBadScFilter = Flag2_eeBadScFilter;
    MuonSystem->Flag2_all = Flag2_HBHENoiseFilter && Flag2_HBHEIsoNoiseFilter && Flag2_BadPFMuonFilter && Flag2_globalSuperTightHalo2016Filter && Flag2_EcalDeadCellTriggerPrimitiveFilter;
    if (isData) MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_eeBadScFilter;

    if (analysisTag!="Razor2016_07Aug2017Rereco")
      {
        MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_ecalBadCalibFilter;
      }

    // if (!MuonSystem->Flag2_all) continue;
    //*************************************************************************
    //Start Object Selection
    //*************************************************************************

    std::vector<leptons> Leptons;
    //-------------------------------
    //Muons
    //-------------------------------
    for( int i = 0; i < nMuons; i++ ) {

      if(!isMuonPOGTightMuon(i)) continue;
      if(muonPt[i] < 25) continue;
      if(fabs(muonEta[i]) > 2.4) continue;

      //remove overlaps
      bool overlap = false;
      for(auto& lep : Leptons)
        {
          if (RazorAnalyzer::deltaR(muonEta[i],muonPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
      if(overlap) continue;

      leptons tmpMuon;
      tmpMuon.lepton.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i], MU_MASS);
      tmpMuon.pdgId = 13 * -1 * muonCharge[i];
      tmpMuon.dZ = muon_dZ[i];
      tmpMuon.passId = isMuonPOGTightMuon(i);
      float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];

      tmpMuon.passLooseIso = muonIso<0.25;
      tmpMuon.passTightIso = muonIso<0.15;
      tmpMuon.passVTightIso = muonIso<0.10;
      tmpMuon.passVVTightIso = muonIso<0.05;

      tmpMuon.passVetoId = false;
      Leptons.push_back(tmpMuon);
    }

    //-------------------------------
    //Electrons
    //-------------------------------
    for( int i = 0; i < nElectrons; i++ ) {

      if (!isEGammaPOGTightElectron(i, true, false, true, "vid")) continue;

      if(elePt[i] < 35) continue;
      if(fabs(eleEta[i]) > 2.4) continue;

      //remove overlaps
      bool overlap = false;
      for(auto& lep : Leptons)
        {
          if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
      if(overlap) continue;
      leptons tmpElectron;
      tmpElectron.lepton.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i], ELE_MASS);
      tmpElectron.pdgId = 11 * -1 * eleCharge[i];
      tmpElectron.dZ = ele_dZ[i];
      tmpElectron.passId = isEGammaPOGTightElectron(i, true, true, true, "Summer16");
      tmpElectron.passVetoId = isEGammaPOGVetoElectron(i, true, true, true, "Summer16");
      Leptons.push_back(tmpElectron);
    }

    sort(Leptons.begin(), Leptons.end(), my_largest_pt);

    //cout<< "Count leptons " << endl;

    // set the number of leptons to 0 initially
    MuonSystem->nLeptons = 0;

    for (const auto& tmp : Leptons ) {
      MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
      MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
      MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
      MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
      MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
      MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
      MuonSystem->lepPassId[MuonSystem->nLeptons] = tmp.passId;
      MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
      MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
      MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
      MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;
      MuonSystem->nLeptons++;
    }
    for(int i = 0; i < nMuons; i++) {
      if (fabs(muonEta[i]>3.0)) continue;
      float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
        MuonSystem->muonPt[i]   = muonPt[i];
        MuonSystem->muonPt[i]   = muonE[i];
        MuonSystem->muonEta[i]  = muonEta[i];
        MuonSystem->muonPhi[i]  = muonPhi[i];
        MuonSystem->muonIso[i]  = muonIso;
        MuonSystem->muonIsGlobal[i]  = muon_isGlobal[i];
        MuonSystem->muonTightId[i]  = isMuonPOGTightMuon(i);
        MuonSystem->muonLooseId[i]  = isMuonPOGLooseMuon(i);
    }

    MuonSystem->nMuons = nMuons;
    MuonSystem->category = MuonSystem->nLeptons;


    TLorentzVector met;
    met.SetPtEtaPhiE(metType1Pt,0,metType1Phi,metType1Pt);
    if ( Leptons.size() > 0 ) {
      TLorentzVector visible = Leptons[0].lepton;
      MuonSystem->MT = GetMT(visible,met);
    }


    //-----------------------------------------------
    //Select Jets
    //-----------------------------------------------


    //cout<< "Select jets " << endl;

    std::vector<jets> Jets;
    float MetXCorr_JESUp = 0.;
    float MetYCorr_JESUp = 0.;
    float MetXCorr_JESDown = 0.;
    float MetYCorr_JESDown = 0.;
    float MetXCorr_HEM = 0.;
    float MetYCorr_HEM = 0.;
    float MetXCorr_EENoise = 0.;
    float MetYCorr_EENoise = 0.;

    for(int i = 0; i < nJets; i++) {

      //------------------------------------------------------------
      //exclude selected muons and electrons from the jet collection
      //------------------------------------------------------------
      double deltaR = -1;
      for(auto& lep : Leptons){
        double thisDR = RazorAnalyzer::deltaR(jetEta[i],jetPhi[i],lep.lepton.Eta(),lep.lepton.Phi());
        if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
      }
      if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

      //------------------------------------------------------------
      //Apply Jet Energy and Resolution Corrections
      //------------------------------------------------------------
      // double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
      // fixedGridRhoFastjetAll, jetJetArea[i] , JetCorrector);
      // cout<<"before JEC"<<endl;
      // double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
      //  				   fixedGridRhoFastjetAll, jetJetArea[i],
      //  				   runNum,
      //  				   JetCorrectorIOV,JetCorrector);
      double JEC = 1.0;
      double jetCorrPt = jetPt[i]*JEC;
      double jetCorrE = jetE[i]*JEC;
      TLorentzVector thisJet = makeTLorentzVector( jetCorrPt, jetEta[i], jetPhi[i], jetCorrE );



      if (thisJet.Eta()>-3.0 && thisJet.Eta()<-1.3 && thisJet.Phi() >-1.57 && thisJet.Phi() <-0.87 && analysisTag == "Razor2018_17SeptEarlyReReco")
        {
          MetXCorr_HEM += thisJet.Px();
          MetYCorr_HEM += thisJet.Py();
        }
      if (fabs(thisJet.Eta())> 2.65 && fabs(thisJet.Eta())<3.139 && thisJet.Pt() < 50  && analysisTag == "Razor2017_17Nov2017Rereco")
        {
          MetXCorr_EENoise += thisJet.Px();
          MetYCorr_EENoise += thisJet.Py();
          if (eventNum==123969624)
            {
              cout<<i<<", "<<thisJet.Eta()<<","<<jetEta[i]<<endl;
              cout<<MetXCorr_EENoise<<", "<<MetYCorr_EENoise<<", "<<thisJet.Px()<<", "<<thisJet.Py()<<","<<thisJet.Pt()<<","<<jetPt[i]<<endl;
            }
        }
      if (fabs(thisJet.Eta())> 2.25 && fabs(thisJet.Eta())<3.0 && thisJet.Pt() > 100 && (analysisTag == "Razor2016_07Aug2017Rereco" || analysisTag == "Razor2017_17Nov2017Rereco"))
        {
          MuonSystem->EE_prefiring = false;
        }


      if (fabs(thisJet.Eta()) >= 3.0)continue;

      jets tmpJet;
      tmpJet.jet    = thisJet;
      tmpJet.passId = jetPassIDTight[i];

//      tmpJet.jetPassMuFrac = jetPassMuFrac[i];
//      tmpJet.jetNeutralHadronEnergyFraction = jetNeutralHadronEnergyFraction[i];
//      tmpJet.jetNeutralEMEnergyFraction = jetNeutralEMEnergyFraction[i];
//      tmpJet.jetChargedEMEnergyFraction = jetChargedEMEnergyFraction[i];
//      tmpJet.jetChargedHadronEnergyFraction = jetChargedHadronEnergyFraction[i];


      // calculate jet energy scale uncertainty
      double unc = helper->getJecUnc( jetCorrPt, jetEta[i], runNum ); //use run=999 as default
      tmpJet.jetPtJESUp = jetCorrPt*(1+unc);
      tmpJet.jetPtJESDown = jetCorrPt*(1-unc);
      tmpJet.jetEJESUp = jetCorrE*(1+unc);
      tmpJet.jetEJESDown = jetCorrE*(1-unc);
      tmpJet.JecUnc = unc;
      TLorentzVector thisJetJESUp = makeTLorentzVector(tmpJet.jetPtJESUp, jetEta[i], jetPhi[i], tmpJet.jetEJESUp);
      TLorentzVector thisJetJESDown = makeTLorentzVector(tmpJet.jetPtJESDown, jetEta[i], jetPhi[i], tmpJet.jetEJESDown);
      if (tmpJet.jetPtJESUp > 10) {
        MetXCorr_JESUp += -1 * (thisJetJESUp.Px() - thisJet.Px());
        MetYCorr_JESUp += -1 * (thisJetJESUp.Py() - thisJet.Py());
      }
      if (tmpJet.jetPtJESDown > 10) {
        MetXCorr_JESDown += -1 * (thisJetJESDown.Px() - thisJet.Px());
        MetYCorr_JESDown += -1 * (thisJetJESDown.Py() - thisJet.Py());
      }
      // if ( !jetPassIDLoose[i] ) continue;
      if(!isPFTightJet(i, true,analysisTag))continue;
      if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
      Jets.push_back(tmpJet);
    }

    //-----------------------------
    //Require at least 2 jets
    //-----------------------------

    //cout<< "Count jets " << endl;

    MuonSystem->nJets = 0;

    // if( Jets.size() < 2 ) continue;
    // if (triggered) trig_lepId_dijet->Fill(1);
    sort(Jets.begin(), Jets.end(), my_largest_pt_jet);
    if (Jets.size()>0) {
      MuonSystem->jetMet_dPhi = RazorAnalyzer::deltaPhi(jetPhi[0],metType1Phi);
    }
    else{
      MuonSystem->jetMet_dPhi = -999.;
    }
    double jetMet_dPhiMin_temp = 999.;
    double jetMet_dPhiMin4_temp = 999.;
    for ( auto &tmp : Jets ) {
      // if(tmp.jet.Pt()<50)continue;
      if(abs(tmp.jet.Eta())>2.4)continue;
      MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
      MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
      MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
      MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
      MuonSystem->jetTime[MuonSystem->nJets] = tmp.time;
      MuonSystem->jetPtJESUp[MuonSystem->nJets] = tmp.jetPtJESUp;
      MuonSystem->jetPtJESDown[MuonSystem->nJets] = tmp.jetPtJESDown;
      MuonSystem->jetEJESUp[MuonSystem->nJets] = tmp.jetEJESUp;
      MuonSystem->jetEJESDown[MuonSystem->nJets] = tmp.jetEJESDown;
      MuonSystem->JecUnc[MuonSystem->nJets] = tmp.JecUnc;
      if (jetMet_dPhiMin4_temp > abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi)) && MuonSystem->nJets < 4) {
        jetMet_dPhiMin4_temp = abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi));

      }
      if (jetMet_dPhiMin_temp > abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi))) {
        if (tmp.jet.Pt()>30 && abs(tmp.jet.Eta())<2.4) {
          jetMet_dPhiMin_temp = abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi));
        }
      }
      // if (tmp.jet.Pt()>20 )MuonSystem->HT = MuonSystem->HT + tmp.jet.Pt();
      MuonSystem->jetTightPassId[MuonSystem->nJets] = tmp.passId;
      // MuonSystem->jetChargedEMEnergyFraction[MuonSystem->nJets] = tmp.jetChargedEMEnergyFraction;
      // MuonSystem->jetNeutralEMEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralEMEnergyFraction;
      // MuonSystem->jetChargedHadronEnergyFraction[MuonSystem->nJets] = tmp.jetChargedHadronEnergyFraction;
      // MuonSystem->jetNeutralHadronEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralHadronEnergyFraction;
      // MuonSystem->jetPassMuFrac[MuonSystem->nJets] = tmp.jetPassMuFrac;

      // float min_deltaR = 15.;
      // int index = 999;
      // for(int i=0; i < nGenJets; i++) {
      //
      //   double current_delta_r = RazorAnalyzer::deltaPhi(genJetPhi[i],jetPhi[MuonSystem->nJets]);
      //   if (current_delta_r < min_deltaR) {
      //     min_deltaR = current_delta_r;
      //     index = i;
      //   }
      // }
      // if (min_deltaR < 0.4) {
      //   MuonSystem->jet_match_genJet_minDeltaR[MuonSystem->nJets] = min_deltaR;
      //   MuonSystem->jet_match_genJet_index[MuonSystem->nJets] = index;
      //   MuonSystem->jet_match_genJet_pt[MuonSystem->nJets] = jetPt[index];
      // }
      MuonSystem->nJets++;
    }


    //cout<< "Select MET " << endl;

    MuonSystem-> jetMet_dPhiMin = jetMet_dPhiMin_temp;
    MuonSystem-> jetMet_dPhiMin4 = jetMet_dPhiMin4_temp;
    // TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metType1Pt, 0, metType1Phi, 0);
    TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(MuonSystem->metXYCorr, 0, MuonSystem->metPhiXYCorr, 0);

    //JES up
    float PFMetXJESUp   = PFMET.Px() + MetXCorr_JESUp;
    float PFMetYJESUp   = PFMET.Py() + MetYCorr_JESUp;
    MuonSystem->metJESUp    = sqrt( pow(PFMetXJESUp,2) + pow(PFMetYJESUp,2) );
    MuonSystem->metPhiJESUp    = atan(PFMetYJESUp/PFMetXJESUp);
    if  (PFMetXJESUp < 0.0) MuonSystem->metPhiJESUp = RazorAnalyzer::deltaPhi(TMath::Pi() + MuonSystem->metPhiJESUp,0.0);
    MuonSystem->metJESUpSF = helper->getMetTriggerSF(MuonSystem->metJESUp)/MuonSystem->metSF;

    //JES down
    float PFMetXJESDown   = PFMET.Px() + MetXCorr_JESDown;
    float PFMetYJESDown   = PFMET.Py() + MetYCorr_JESDown;
    MuonSystem->metJESDown    =  sqrt( pow(PFMetXJESDown,2) + pow(PFMetYJESDown,2) );
    MuonSystem->metPhiJESDown    = atan(PFMetYJESDown/PFMetXJESDown);
    if  (PFMetXJESUp < 0.0) MuonSystem->metPhiJESDown = RazorAnalyzer::deltaPhi(TMath::Pi() + MuonSystem->metPhiJESDown,0.0);
    MuonSystem->metJESDownSF = helper->getMetTriggerSF(MuonSystem->metJESDown)/MuonSystem->metSF;

    //HEM
    float PFMetXHEM   = PFMET.Px() + MetXCorr_HEM;
    float PFMetYHEM   = PFMET.Py() + MetYCorr_HEM;
    MuonSystem->metHEM    = sqrt( pow(PFMetXHEM,2) + pow(PFMetYHEM,2) );
    MuonSystem->metPhiHEM    = atan(PFMetYHEM/PFMetXHEM);
    if  (PFMetXHEM < 0.0) MuonSystem->metPhiHEM = RazorAnalyzer::deltaPhi(TMath::Pi() + MuonSystem->metPhiHEM,0.0);


    //EENoise
    float PFMetXEENoise   = PFMET.Px() + MetXCorr_EENoise;
    float PFMetYEENoise   = PFMET.Py() + MetYCorr_EENoise;
    MuonSystem->metEENoise    = sqrt( pow(PFMetXEENoise,2) + pow(PFMetYEENoise,2) );
    MuonSystem->metPhiEENoise    = atan(PFMetYEENoise/PFMetXEENoise);
    if  (PFMetXEENoise < 0.0) MuonSystem->metPhiEENoise = RazorAnalyzer::deltaPhi(TMath::Pi() + MuonSystem->metPhiEENoise,0.0);

    // if (!MuonSystem->SingleLepTrigger) continue;
    // if (MuonSystem->metEENoise < 200) continue;

    if(signalScan && !isData)accep_met2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->pileupWeight*MuonSystem->metSF);
    else if(!isData) accep_met->Fill(1.0, genWeight*MuonSystem->pileupWeight*MuonSystem->metSF);
    else Nmet200->Fill(1.0);

    //cout<< "Rechit clustering " << endl;

    // DT rechit clustering
    vector<Point> dt_points;
    vector<int> dtRechitsClusterId;
    dt_points.clear();
    MuonSystem->nDTRechits  = 0;

    for (int i = 0; i < nDtRechits; i++) {
      Point p;
      p.phi = dtRechitPhi[i];
      p.eta = dtRechitEta[i];
      p.x = dtRechitX[i];
      p.y = dtRechitY[i];
      p.z = dtRechitZ[i];
      // the timing is actually not used in the clustering, the distance metric only uses eta and phi
      // p.t = dtRechitTime[i];
      p.station = dtRechitStation[i];
      // p.wheel = dtRechitWheel[i];
      // p.superlayer = 0;
      p.clusterID = UNCLASSIFIED;
      dt_points.push_back(p);
      dtRechitsClusterId.push_back(-1);

        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -2) MuonSystem->nDTRechitsChamberMinus12++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -1) MuonSystem->nDTRechitsChamberMinus11++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 0) MuonSystem->nDTRechitsChamber10++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 1) MuonSystem->nDTRechitsChamberPlus11++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 2) MuonSystem->nDTRechitsChamberPlus12++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -2) MuonSystem->nDTRechitsChamberMinus22++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -1) MuonSystem->nDTRechitsChamberMinus21++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 0) MuonSystem->nDTRechitsChamber20++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 1) MuonSystem->nDTRechitsChamberPlus21++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 2) MuonSystem->nDTRechitsChamberPlus22++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -2) MuonSystem->nDTRechitsChamberMinus32++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -1) MuonSystem->nDTRechitsChamberMinus31++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 0) MuonSystem->nDTRechitsChamber30++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 1) MuonSystem->nDTRechitsChamberPlus31++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 2) MuonSystem->nDTRechitsChamberPlus32++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -2) MuonSystem->nDTRechitsChamberMinus42++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -1) MuonSystem->nDTRechitsChamberMinus41++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 0) MuonSystem->nDTRechitsChamber40++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 1) MuonSystem->nDTRechitsChamberPlus41++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 2) MuonSystem->nDTRechitsChamberPlus42++;
        if (dtRechitStation[i] == 1) MuonSystem->nDTRechitsStation1++;
        if (dtRechitStation[i] == 2) MuonSystem->nDTRechitsStation2++;
        if (dtRechitStation[i] == 3) MuonSystem->nDTRechitsStation3++;
        if (dtRechitStation[i] == 4) MuonSystem->nDTRechitsStation4++;

        if (dtRechitWheel[i] == -2) MuonSystem->nDTRechitsWheelMinus2++;
        if (dtRechitWheel[i] == -1) MuonSystem->nDTRechitsWheelMinus1++;
        if (dtRechitWheel[i] == 0) MuonSystem->nDTRechitsWheel0++;
        if (dtRechitWheel[i] == 1) MuonSystem->nDTRechitsWheelPlus1++;
        if (dtRechitWheel[i] == 2) MuonSystem->nDTRechitsWheelPlus2++;

      MuonSystem->nDTRechits++;
    }

      if ( MuonSystem->nDTRechitsChamberMinus12 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberMinus11 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamber10 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus11 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus12 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberMinus22 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberMinus21 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamber20 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus21 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus22 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberMinus32 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberMinus31 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamber30 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus31 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus32 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberMinus42 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberMinus41 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamber40 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus41 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus42 > 50) MuonSystem->nDtRings++;

      if (MuonSystem->nDTRechitsStation1 > 25) MuonSystem->nDtStations25++;
      if (MuonSystem->nDTRechitsStation2 > 25) MuonSystem->nDtStations25++;
      if (MuonSystem->nDTRechitsStation3 > 25) MuonSystem->nDtStations25++;
      if (MuonSystem->nDTRechitsStation4 > 25) MuonSystem->nDtStations25++;

      if (MuonSystem->nDTRechitsWheelMinus2 > 25) MuonSystem->nDtWheels25++;
      if (MuonSystem->nDTRechitsWheelMinus1 > 25) MuonSystem->nDtWheels25++;
      if (MuonSystem->nDTRechitsWheel0 > 25) MuonSystem->nDtWheels25++;
      if (MuonSystem->nDTRechitsWheelPlus1 > 25) MuonSystem->nDtWheels25++;
      if (MuonSystem->nDTRechitsWheelPlus2 > 25) MuonSystem->nDtWheels25++;

    //Do DT DBSCAN Clustering
    int min_dt_point = 50;  //minimum number of segments to call it a cluster
    float dt_epsilon = 0.2; //cluster radius parameter
    DBSCAN dtds(min_dt_point, dt_epsilon, dt_points);

    dtds.run();
    dtds.result();
    dtds.clusterMoments();
    dtds.sort_clusters();
    dtds.merge_clusters();
    dtds.result();
    dtds.clusterMoments();
    dtds.sort_clusters();


    //cout << "Do CSC rechit clustering" << endl;

    // CSC rechit clustering
    vector<Point> points;
    vector<int> cscRechitsClusterId;
    points.clear();
    cscRechitsClusterId.clear();
    MuonSystem->nCscRechits  = 0;
    float cscPosRechitsAvgT =0.0;
    float cscNegRechitsAvgT =0.0;

    for (int i = 0; i < ncscRechits; i++) {

      //pick out the right bits for chamber
      Point p;
      p.phi = cscRechitsPhi[i];
      p.eta = cscRechitsEta[i];
      p.x = cscRechitsX[i];
      p.y = cscRechitsY[i];
      p.z = cscRechitsZ[i];
      p.t = cscRechitsTpeak[i];
      p.twire = cscRechitsTwire[i];
      p.station = cscRechitsStation[i];
      p.chamber = cscRechitsChamber[i];
      p.clusterID = UNCLASSIFIED;
      points.push_back(p);

      cscRechitsClusterId.push_back(-1);
      // these branches were actually not used in the analysis as cuts, but they were used to select for a control region that had a large number of early hits, which should have a lot of cosmic muon showers
      if (cscRechitsY[i]>=0.0)
        {
          MuonSystem->nCscPositiveYRechits++;
          cscPosRechitsAvgT +=  cscRechitsTpeak[i];
        }
      else
        {
          MuonSystem->nCscNegativeYRechits++;
          cscNegRechitsAvgT +=  cscRechitsTpeak[i];
        }
      if (cscRechitsTpeak[i]<-12.5)MuonSystem->nEarlyCscRechits++;
      if (cscRechitsTpeak[i]>12.5)MuonSystem->nLateCscRechits++;
      if (cscRechitsTpeak[i]<-25)MuonSystem->nEarly2CscRechits++;
      if (cscRechitsTpeak[i]>25)MuonSystem->nLate2CscRechits++;
      if (cscRechitsChamber[i] == 11) MuonSystem->nCscRechitsChamberPlus11++;
      if (cscRechitsChamber[i] == 12) MuonSystem->nCscRechitsChamberPlus12++;
      if (cscRechitsChamber[i] == 13) MuonSystem->nCscRechitsChamberPlus13++;
      if (cscRechitsChamber[i] == 21) MuonSystem->nCscRechitsChamberPlus21++;
      if (cscRechitsChamber[i] == 22) MuonSystem->nCscRechitsChamberPlus22++;
      if (cscRechitsChamber[i] == 31) MuonSystem->nCscRechitsChamberPlus31++;
      if (cscRechitsChamber[i] == 32) MuonSystem->nCscRechitsChamberPlus32++;
      if (cscRechitsChamber[i] == 41) MuonSystem->nCscRechitsChamberPlus41++;
      if (cscRechitsChamber[i] == 42) MuonSystem->nCscRechitsChamberPlus42++;
      if (cscRechitsChamber[i] == -11) MuonSystem->nCscRechitsChamberMinus11++;
      if (cscRechitsChamber[i] == -12) MuonSystem->nCscRechitsChamberMinus12++;
      if (cscRechitsChamber[i] == -13) MuonSystem->nCscRechitsChamberMinus13++;
      if (cscRechitsChamber[i] == -21) MuonSystem->nCscRechitsChamberMinus21++;
      if (cscRechitsChamber[i] == -22) MuonSystem->nCscRechitsChamberMinus22++;
      if (cscRechitsChamber[i] == -31) MuonSystem->nCscRechitsChamberMinus31++;
      if (cscRechitsChamber[i] == -32) MuonSystem->nCscRechitsChamberMinus32++;
      if (cscRechitsChamber[i] == -41) MuonSystem->nCscRechitsChamberMinus41++;
      if (cscRechitsChamber[i] == -42) MuonSystem->nCscRechitsChamberMinus42++;
      MuonSystem->nCscRechits++;
    }

    //TODO: Check what to fill with  cscPosTpeak
    MuonSystem->cscPosTpeak = float(cscPosRechitsAvgT/ncscRechits);
    MuonSystem->cscNegTpeak = float(cscNegRechitsAvgT/ncscRechits);
    MuonSystem->nCscRings = 0;
    if ( MuonSystem->nCscRechitsChamberPlus11 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus12 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus13 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus21 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus22 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus31 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus32 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus41 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus42 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus11 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus12 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus13 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus21 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus22 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus31 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus32 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus41 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus42 > 50) MuonSystem->nCscRings++;

    //MuonSystem->tree_->Fill();
    //cout << "Do CSC DBScan" << endl;
    //Do DBSCAN Clustering
    int min_point = 50;  //minimum number of segments to call it a cluster
    float epsilon = 0.2; //cluster radius parameter
    DBSCAN ds(min_point, epsilon, points);
    ds.run();
    ds.result();
    ds.clusterMoments();
    ds.sort_clusters();

    ds.merge_clusters();
    ds.result();
    ds.clusterMoments();
    ds.sort_clusters();

    //cout << "CSC cluster analysis" << endl;

    MuonSystem->nCscRechitClusters3 = 0;
    for ( auto &tmp : ds.CscCluster ) {

      MuonSystem->cscRechitCluster3X[MuonSystem->nCscRechitClusters3] =tmp.x;
      MuonSystem->cscRechitCluster3Y[MuonSystem->nCscRechitClusters3] =tmp.y;
      MuonSystem->cscRechitCluster3Z[MuonSystem->nCscRechitClusters3] =tmp.z;
      MuonSystem->cscRechitCluster3Size[MuonSystem->nCscRechitClusters3] =tmp.nCscSegments;
      MuonSystem->cscRechitCluster3Time[MuonSystem->nCscRechitClusters3] = tmp.t;
      MuonSystem->cscRechitCluster3TimeTotal[MuonSystem->nCscRechitClusters3] = tmp.tTotal;
      MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3] =tmp.eta;
      MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3] = tmp.phi;
      MuonSystem->cscRechitCluster3MajorAxis[MuonSystem->nCscRechitClusters3] =tmp.MajorAxis;
      MuonSystem->cscRechitCluster3MinorAxis[MuonSystem->nCscRechitClusters3] =tmp.MinorAxis;
      MuonSystem->cscRechitCluster3XSpread[MuonSystem->nCscRechitClusters3] =tmp.XSpread;
      MuonSystem->cscRechitCluster3YSpread[MuonSystem->nCscRechitClusters3] =tmp.YSpread;
      MuonSystem->cscRechitCluster3ZSpread[MuonSystem->nCscRechitClusters3] =tmp.ZSpread;
      MuonSystem->cscRechitCluster3EtaPhiSpread[MuonSystem->nCscRechitClusters3] =tmp.EtaPhiSpread;
      MuonSystem->cscRechitCluster3XYSpread[MuonSystem->nCscRechitClusters3] =tmp.XYSpread;
      MuonSystem->cscRechitCluster3RSpread[MuonSystem->nCscRechitClusters3] =tmp.RSpread;


      MuonSystem->cscRechitCluster3EtaSpread[MuonSystem->nCscRechitClusters3] =tmp.EtaSpread;
      MuonSystem->cscRechitCluster3PhiSpread[MuonSystem->nCscRechitClusters3] = tmp.PhiSpread;
      MuonSystem->cscRechitCluster3TimeSpread[MuonSystem->nCscRechitClusters3] = tmp.TSpread;

      MuonSystem->cscRechitCluster3MaxChamber[MuonSystem->nCscRechitClusters3] = tmp.maxChamber;
      MuonSystem->cscRechitCluster3NChamber[MuonSystem->nCscRechitClusters3] = tmp.nChamber;
      MuonSystem->cscRechitCluster3MaxStation[MuonSystem->nCscRechitClusters3] = tmp.maxStation;
      MuonSystem->cscRechitCluster3NStation10[MuonSystem->nCscRechitClusters3] = tmp.nStation10;
      MuonSystem->cscRechitCluster3AvgStation10[MuonSystem->nCscRechitClusters3] = tmp.avgStation10;


      MuonSystem->cscRechitCluster3Me11Ratio[MuonSystem->nCscRechitClusters3] = tmp.Me11Ratio;
      MuonSystem->cscRechitCluster3Me12Ratio[MuonSystem->nCscRechitClusters3] = tmp.Me12Ratio;

      //Jet veto/ muon veto
      MuonSystem->cscRechitCluster3JetVetoPt[MuonSystem->nCscRechitClusters3]     = -999.0;
      MuonSystem->cscRechitCluster3JetVetoE[MuonSystem->nCscRechitClusters3]      = -999.0;
      MuonSystem->cscRechitCluster3MuonVetoPt[MuonSystem->nCscRechitClusters3]    = -999.0;
      MuonSystem->cscRechitCluster3MuonVetoE[MuonSystem->nCscRechitClusters3]     = -999.0;
      MuonSystem->cscRechitCluster3GenMuonVetoPt[MuonSystem->nCscRechitClusters3] = -999.0;
      MuonSystem->cscRechitCluster3GenMuonVetoE[MuonSystem->nCscRechitClusters3]  = -999.0;
      MuonSystem->cscRechitCluster3IsoMuonVetoPt[MuonSystem->nCscRechitClusters3] = -999.0;

      // jet veto
      //cout << "Check jet veto" <<endl;
      for(int i = 0; i < nJets; i++) {
        if (fabs(jetEta[i]>3.0)) continue;
        if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4
            && jetPt[i] > MuonSystem->cscRechitCluster3JetVetoPt[MuonSystem->nCscRechitClusters3] ) {
          MuonSystem->cscRechitCluster3JetVetoPt[MuonSystem->nCscRechitClusters3]  = jetPt[i];
          MuonSystem->cscRechitCluster3JetVetoEta[MuonSystem->nCscRechitClusters3]  = jetEta[i];
          MuonSystem->cscRechitCluster3JetVetoPhi[MuonSystem->nCscRechitClusters3]  = jetPhi[i];

        }
        if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3], MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4
            && jetE[i] > MuonSystem->cscRechitCluster3JetVetoE[MuonSystem->nCscRechitClusters3] ) {
          MuonSystem->cscRechitCluster3JetVetoE[MuonSystem->nCscRechitClusters3]  = jetE[i];
        }
        }
      float min_deltaR = 15.;
      int index = 999;


      //cout << "Check muon veto" <<endl;
      for(int i = 0; i < nMuons; i++) {
        if (fabs(muonEta[i]>3.0)) continue;
        float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
        if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.8
            && muonPt[i] > MuonSystem->cscRechitCluster3MuonVetoPt[MuonSystem->nCscRechitClusters3] ) {
          MuonSystem->cscRechitCluster3MuonVetoPt[MuonSystem->nCscRechitClusters3]  = muonPt[i];
          MuonSystem->cscRechitCluster3MuonVetoE[MuonSystem->nCscRechitClusters3]  = muonE[i];
          MuonSystem->cscRechitCluster3MuonVetoPhi[MuonSystem->nCscRechitClusters3]  = muonPhi[i];
          MuonSystem->cscRechitCluster3MuonVetoEta[MuonSystem->nCscRechitClusters3]  = muonEta[i];
          MuonSystem->cscRechitCluster3MuonVetoLooseIso[MuonSystem->nCscRechitClusters3]  = muonIso<0.25;
          MuonSystem->cscRechitCluster3MuonVetoTightIso[MuonSystem->nCscRechitClusters3]  = muonIso<0.15;
          MuonSystem->cscRechitCluster3MuonVetoVTightIso[MuonSystem->nCscRechitClusters3]  = muonIso<0.10;
          MuonSystem->cscRechitCluster3MuonVetoVVTightIso[MuonSystem->nCscRechitClusters3]  = muonIso<0.05;
          MuonSystem->cscRechitCluster3MuonVetoTightId[MuonSystem->nCscRechitClusters3]  = isMuonPOGTightMuon(i);
          MuonSystem->cscRechitCluster3MuonVetoLooseId[MuonSystem->nCscRechitClusters3]  = isMuonPOGLooseMuon(i);
        }
      }

      //cout << "match to gen-level muon" <<endl;

      // match to gen-level muon
      if(!isData) {
        for(int i = 0; i < nGenParticle; i++) {
          if (abs(gParticleId[i])!=13) continue;
          if (abs(gParticleMotherId[i])>24 || abs(gParticleMotherId[i])<23) continue;
          if (abs(gParticleStatus[i])!=1)continue;
          // if (fabs(muonEta[i]>3.0)) continue;
          if (RazorAnalyzer::deltaR(gParticleEta[i], gParticlePhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 ) {
            MuonSystem->cscRechitCluster3GenMuonVetoPt[MuonSystem->nCscRechitClusters3]  = gParticlePt[i];
            MuonSystem->cscRechitCluster3GenMuonVetoE[MuonSystem->nCscRechitClusters3]  = gParticleE[i];
          }
        }
        min_deltaR = 15.;
        index = 999;
        for(int j = 0; j < nGenParticle; j++) {

          double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3], MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3], gParticleEta[j], gParticlePhi[j]);
          if (current_delta_r < min_deltaR) {
            min_deltaR = current_delta_r;
            index = j;
          }
        }
        if (min_deltaR < 0.4) {

          MuonSystem->cscRechitCluster3_match_gParticle[MuonSystem->nCscRechitClusters3] = true;
          MuonSystem->cscRechitCluster3_match_gParticle_minDeltaR[MuonSystem->nCscRechitClusters3] = min_deltaR;
          MuonSystem->cscRechitCluster3_match_gParticle_index[MuonSystem->nCscRechitClusters3] = index;
          MuonSystem->cscRechitCluster3_match_gParticle_id[MuonSystem->nCscRechitClusters3] = gParticleId[index];

          MuonSystem->cscRechitCluster3_match_gParticle_eta[MuonSystem->nCscRechitClusters3] = gParticleEta[index];
          MuonSystem->cscRechitCluster3_match_gParticle_phi[MuonSystem->nCscRechitClusters3] = gParticlePhi[index];
          MuonSystem->cscRechitCluster3_match_gParticle_E[MuonSystem->nCscRechitClusters3] = gParticleE[index];
          MuonSystem->cscRechitCluster3_match_gParticle_pt[MuonSystem->nCscRechitClusters3] = gParticlePt[index];
          MuonSystem->cscRechitCluster3_match_gParticle_MotherId[MuonSystem->nCscRechitClusters3] = gParticleMotherId[index];
        }
        // match to gen level LLP
        min_deltaR = 15.;
        index = 999;
        for(int j = 0; j < 2;j++) {

          double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3], MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3], gLLP_eta[j], gLLP_phi[j]);
          if (current_delta_r < min_deltaR) {
            min_deltaR = current_delta_r;
            index = j;
          }
        }
        if (min_deltaR < 0.4) {
          MuonSystem->cscRechitCluster3_match_gLLP[MuonSystem->nCscRechitClusters3] = true;
          MuonSystem->cscRechitCluster3_match_gLLP_minDeltaR[MuonSystem->nCscRechitClusters3] = min_deltaR;
          MuonSystem->cscRechitCluster3_match_gLLP_eta[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_eta[index];
          MuonSystem->cscRechitCluster3_match_gLLP_phi[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_phi[index];
          MuonSystem->cscRechitCluster3_match_gLLP_decay_r[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_decay_vertex_r[index];
          MuonSystem->cscRechitCluster3_match_gLLP_decay_x[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_decay_vertex_x[index];
          MuonSystem->cscRechitCluster3_match_gLLP_decay_y[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_decay_vertex_y[index];
          MuonSystem->cscRechitCluster3_match_gLLP_decay_z[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_decay_vertex_z[index];
          MuonSystem->cscRechitCluster3_match_gLLP_ctau[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_ctau[index];
          MuonSystem->cscRechitCluster3_match_gLLP_beta[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_beta[index];
          MuonSystem->cscRechitCluster3_match_gLLP_csc[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_csc[index];
        }
      }

      //cout << "match to DT hits and segments" << endl;

      //match to MB1 DT hits
      for (int i = 0; i < nDtRechits; i++) {
        if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 )
          {
            MuonSystem->cscRechitCluster3_match_dtRechits_0p4[MuonSystem->nCscRechitClusters3] ++;
            if (dtRechitStation[i] == 1) MuonSystem->cscRechitCluster3_match_MB1_0p4[MuonSystem->nCscRechitClusters3] ++;
          }

        if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.6 )
          {
            MuonSystem->cscRechitCluster3_match_dtRechits_0p6[MuonSystem->nCscRechitClusters3] ++;
            if (dtRechitStation[i] == 1) MuonSystem->cscRechitCluster3_match_MB1_0p6[MuonSystem->nCscRechitClusters3] ++;
          }

      }
      //match to MB1 DT segments
      for (int i = 0; i < nDtSeg; i++) {
        if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 )
          {
            MuonSystem->cscRechitCluster3_match_dtSeg_0p4[MuonSystem->nCscRechitClusters3] ++;
            if (dtSegStation[i] == 1) MuonSystem->cscRechitCluster3_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters3] ++;
          }

        if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.6 )
          {
            MuonSystem->cscRechitCluster3_match_dtSeg_0p6[MuonSystem->nCscRechitClusters3] ++;
            if (dtSegStation[i] == 1) MuonSystem->cscRechitCluster3_match_MB1Seg_0p6[MuonSystem->nCscRechitClusters3] ++;
          }

      }
      //match to RPC hits in RE1/2
      for (int i = 0; i < nRpc; i++) {
        float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
        if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 )
          {
            if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730)
              {
                MuonSystem->cscRechitCluster3_match_RE12_0p4[MuonSystem->nCscRechitClusters3] ++;
              }
            if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
              {
                MuonSystem->cscRechitCluster3_match_RB1_0p4[MuonSystem->nCscRechitClusters3] ++;
              }

          }
        if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.6 )
          {
            if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730)
              {
                MuonSystem->cscRechitCluster3_match_RE12_0p6[MuonSystem->nCscRechitClusters3] ++;
              }
            if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
              {
                MuonSystem->cscRechitCluster3_match_RB1_0p6[MuonSystem->nCscRechitClusters3] ++;
              }
          }
      }

      MuonSystem->cscRechitCluster3Met_dPhi[MuonSystem->nCscRechitClusters3] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3],MuonSystem->metPhi);
      MuonSystem->cscRechitCluster3MetXYCorr_dPhi[MuonSystem->nCscRechitClusters3] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3],MuonSystem->metPhiXYCorr);
      MuonSystem->cscRechitCluster3MetEENoise_dPhi[MuonSystem->nCscRechitClusters3] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3],MuonSystem->metPhiEENoise);
      MuonSystem->cscRechitCluster3MetHEM_dPhi[MuonSystem->nCscRechitClusters3] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3],MuonSystem->metPhiHEM);

      MuonSystem->nCscRechitClusters3++;
    } // end loop on CSC clusters



    // RPC hits
    MuonSystem->nRpc = nRpc;
    for(int i = 0; i < nRpc; i++) {
        MuonSystem -> rpcX[i]       =rpcX[i]       ;
        MuonSystem -> rpcY[i]       =rpcY[i]       ;
        MuonSystem -> rpcZ[i]       =rpcZ[i]       ;
        MuonSystem -> rpcEta[i]     =rpcEta[i]     ;
        MuonSystem -> rpcPhi[i]     =rpcPhi[i]     ;
        MuonSystem -> rpcBx[i]      =rpcBx[i]      ;
        MuonSystem -> rpcStation[i] =rpcStation[i] ;
        MuonSystem -> rpcRing[i]    =rpcRing[i]    ;
        MuonSystem -> rpcSector[i]  =rpcSector[i]  ;
        MuonSystem -> rpcLayer[i]   =rpcLayer[i]   ;
    }

    // DT cluster



    MuonSystem->nDtRechitClusters = 0;

    for ( auto &tmp : dtds.CscCluster ) {

        MuonSystem->dtRechitClusterX[MuonSystem->nDtRechitClusters] =tmp.x;
        MuonSystem->dtRechitClusterY[MuonSystem->nDtRechitClusters] =tmp.y;
        MuonSystem->dtRechitClusterZ[MuonSystem->nDtRechitClusters] =tmp.z;
        if (abs(tmp.z) < 126.8) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 0;
        else if (tmp.z > 126.8 && tmp.z < 395.4) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 1;
        else if (tmp.z < -126.8 && tmp.z > -395.4)MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -1;
        else if (tmp.z<0) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -2;
        else MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 2;


        MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters] =tmp.eta;
        MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters] =tmp.phi;
        MuonSystem->dtRechitClusterSize[MuonSystem->nDtRechitClusters] = tmp.nCscSegments;

        MuonSystem->dtRechitClusterNSegmentStation1[MuonSystem->nDtRechitClusters] = tmp.nDtSegmentStation1;
        MuonSystem->dtRechitClusterNSegmentStation2[MuonSystem->nDtRechitClusters] = tmp.nDtSegmentStation2;
        MuonSystem->dtRechitClusterNSegmentStation3[MuonSystem->nDtRechitClusters] = tmp.nDtSegmentStation3;
        MuonSystem->dtRechitClusterNSegmentStation4[MuonSystem->nDtRechitClusters] = tmp.nDtSegmentStation4;

        MuonSystem->dtRechitClusterMaxChamber[MuonSystem->nDtRechitClusters] = tmp.maxChamber;
        MuonSystem->dtRechitClusterMaxChamberRatio[MuonSystem->nDtRechitClusters] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
        MuonSystem->dtRechitClusterNChamber[MuonSystem->nDtRechitClusters] = tmp.nChamber;
        MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters] = tmp.maxStation;
        MuonSystem->dtRechitClusterMaxStationRatio[MuonSystem->nDtRechitClusters] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;
        MuonSystem->dtRechitClusterNStation[MuonSystem->nDtRechitClusters] = tmp.nStation;
        MuonSystem->dtRechitClusterNStation10[MuonSystem->nDtRechitClusters] = tmp.nStation10;
        MuonSystem->dtRechitClusterAvgStation10[MuonSystem->nDtRechitClusters] = tmp.avgStation10;

        //Jet veto/ muon veto
        MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] = -999.9;
        MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters]  = -999.9;
        MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters]= -999.9;
        MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = -999.9;


        // jet veto
        for(int i = 0; i < nJets; i++)
        {
          if (fabs(jetEta[i]>3.0)) continue;
          if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] ) {
            MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters]  = jetPt[i];
            MuonSystem->dtRechitClusterJetVetoEta[MuonSystem->nDtRechitClusters]  = jetEta[i];
            MuonSystem->dtRechitClusterJetVetoPhi[MuonSystem->nDtRechitClusters]  = jetPhi[i];
          }
          if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetE[i] > MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] ) {
            MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters]  = jetE[i];
          }
        }


        for(int i = 0; i < nMuons; i++)
        {
          if (fabs(muonEta[i]>3.0)) continue;
          float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
          if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.8 && muonPt[i] > MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] ) {
            MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters]  = muonPt[i];
            MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters]  = muonE[i];
            MuonSystem->dtRechitClusterMuonVetoPhi[MuonSystem->nDtRechitClusters]  = muonPhi[i];
            MuonSystem->dtRechitClusterMuonVetoEta[MuonSystem->nDtRechitClusters]  = muonEta[i];
            MuonSystem->dtRechitClusterMuonVetoGlobal[MuonSystem->nDtRechitClusters]  = muon_isGlobal[i];
            MuonSystem->dtRechitClusterMuonVetoLooseId[MuonSystem->nDtRechitClusters]  = isMuonPOGLooseMuon(i);
          }
        }
        // match to gen-level muon
        float min_deltaR = 15.;
        int index = 999;
        if(!isData)
        {
          for(int i = 0; i < nGenParticle; i++)
          {
            if (abs(gParticleStatus[i])>2)continue;
            double current_delta_r = RazorAnalyzer::deltaR(tmp.eta, tmp.phi, gParticleEta[i], gParticlePhi[i]);
            if (current_delta_r < min_deltaR)
            {
              min_deltaR = current_delta_r;
              index = i;
            }
          }

          MuonSystem->dtRechitCluster_match_gParticle_deltaR[MuonSystem->nDtRechitClusters]  = min_deltaR;
          MuonSystem->dtRechitCluster_match_gParticle_Id[MuonSystem->nDtRechitClusters]  = gParticleId[index];
          MuonSystem->dtRechitCluster_match_gParticle_Pt[MuonSystem->nDtRechitClusters]  = gParticlePt[index];
          MuonSystem->dtRechitCluster_match_gParticle_Eta[MuonSystem->nDtRechitClusters]  = gParticleEta[index];
          MuonSystem->dtRechitCluster_match_gParticle_Phi[MuonSystem->nDtRechitClusters]  = gParticlePhi[index];
          MuonSystem->dtRechitCluster_match_gParticle_E[MuonSystem->nDtRechitClusters]  = gParticleE[index];
          MuonSystem->dtRechitCluster_match_gParticle_Status[MuonSystem->nDtRechitClusters]  = gParticleStatus[index];
          MuonSystem->dtRechitCluster_match_gParticle_MotherId[MuonSystem->nDtRechitClusters]  = gParticleMotherId[index];

        }



        min_deltaR = 15.;
        index = 999;
        for(int j = 0; j < 2;j++)
        {
          double current_delta_r = RazorAnalyzer::deltaR(tmp.eta, tmp.phi, gLLP_eta[j], gLLP_phi[j]);
          if (current_delta_r < min_deltaR)
          {
            min_deltaR = current_delta_r;
            index = j;
          }
        }
        if (min_deltaR < 0.4)MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = true;
        else MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = false;

         MuonSystem->dtRechitCluster_match_gLLP_minDeltaR[MuonSystem->nDtRechitClusters] = min_deltaR;
         MuonSystem->dtRechitCluster_match_gLLP_eta[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_eta[index];
         MuonSystem->dtRechitCluster_match_gLLP_phi[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_phi[index];
         MuonSystem->dtRechitCluster_match_gLLP_decay_x[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_x[index];
         MuonSystem->dtRechitCluster_match_gLLP_decay_y[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_y[index];
         MuonSystem->dtRechitCluster_match_gLLP_decay_z[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_z[index];
         MuonSystem->dtRechitCluster_match_gLLP_ctau[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_ctau[index];
         MuonSystem->dtRechitCluster_match_gLLP_beta[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_beta[index];
         MuonSystem->dtRechitCluster_match_gLLP_csc[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_csc[index];
         MuonSystem->dtRechitCluster_match_gLLP_dt[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_dt[index];
         MuonSystem->dtRechitCluster_match_gLLP_e[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_e[index];


        //match to MB1 DT hits

        for (int i = 0; i < nDtRechits; i++) {
          if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.5 )
          {
            if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1hits_0p5[MuonSystem->nDtRechitClusters] ++;
          }
          if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
          {
            if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1hits_0p4[MuonSystem->nDtRechitClusters] ++;
          }
          if(abs(dtRechitWheel[i]-MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])==1 && dtRechitStation[i] == 1)
          {
            if (abs(RazorAnalyzer::deltaPhi(dtRechitPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < TMath::Pi()/4.0 )
            {
              if (dtRechitWheel[i]-MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] == 1)
              {
                MuonSystem->dtRechitCluster_match_MB1hits_cosmics_plus[MuonSystem->nDtRechitClusters] ++;
              }
              else
              {
                MuonSystem->dtRechitCluster_match_MB1hits_cosmics_minus[MuonSystem->nDtRechitClusters] ++;
              }
            }
          }
        }


        std::vector<int> dtRechitCluster_match_rpcBx;

        //match to RPC hits with dPhi<0.5 and same wheel in DT
        for (int i = 0; i < nRpc; i++) {
          float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
          if (rpcRegion[i]!=0) continue;
          if (abs(RazorAnalyzer::deltaPhi(rpcPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5 )
          {
            if (rpcRing[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])
            {
              dtRechitCluster_match_rpcBx.push_back(rpcBx[i]);
              MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters]++;
              if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_dPhi0p5[MuonSystem->nDtRechitClusters] ++;

            }
          }
          if(RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
          {
            if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_0p4[MuonSystem->nDtRechitClusters] ++;
          }


        }

        int max_occurence = 0;
        int max_bx = -999;
        for (unsigned int l = 0; l < dtRechitCluster_match_rpcBx.size(); l++)
        {
          int counter = 0;
          for(unsigned int j = 0; j < dtRechitCluster_match_rpcBx.size(); j ++)
          {
            if (dtRechitCluster_match_rpcBx[j] == dtRechitCluster_match_rpcBx[l]) counter++;
          }
          if (counter>max_occurence)
          {
            max_occurence = counter;
            max_bx = dtRechitCluster_match_rpcBx[l];
          }
          // cout<<dtRechitCluster_match_rpcBx[l]<<endl;
        }
        // cout<<max_occurence<<","<<max_bx<<endl;
        // cout<<"###################################################"<<endl;
          MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5[MuonSystem->nDtRechitClusters] = max_bx;

        MuonSystem->dtRechitClusterMetEENoise_dPhi[MuonSystem->nDtRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->metPhiEENoise);
        if (MuonSystem->nLeptons > 0)
        {
          MuonSystem->dtRechitClusterLep_dPhi[MuonSystem->nDtRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->lepPhi[0]);

        }





        MuonSystem->nDtRechitClusters++;
      }//end loop for dt clusters
      // NCLUSTER
      if(isData && MuonSystem->nDtRechitClusters == 0 && MuonSystem->nCscRechitClusters3 == 0)continue;
      if(MuonSystem->nDtRechitClusters == 0 && MuonSystem->nCscRechitClusters3 == 0)continue;


    if(!isData && signalScan) {
      pair<int,int> smsPair = make_pair(MuonSystem->mX, MuonSystem->ctau);
      Trees2D[smsPair]->Fill();
    }
    else {
      if (MuonSystem->tree_) {
        //cout << "Fill tree" << endl;
        MuonSystem->tree_->Fill();
      }
      else {
        cout << "No MuonSystem tree" << endl;
      }
    }
  } // end event loop

  // writing output
  if(!isData && signalScan) {
    for(auto &filePtr : Files2D) {
      cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
      filePtr.second->cd();
      Trees2D[filePtr.first]->Write();
      NEvents2D[filePtr.first]->Write("NEvents");
      Total2D[filePtr.first]->Write("Total");
      accep2D[filePtr.first]->Write("acceptance");
      accep_met2D[filePtr.first]->Write("acceptance_met");
      filePtr.second->Close();
    }
  }
  else if (!isData) {
    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->cd();
    MuonSystem->tree_->Write();
    NEvents->Write();
    accep->Write("acceptance");
    accep_met->Write("acceptance_met");
    outFile->Close();
  }
  else {
    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->cd();
    MuonSystem->tree_->Write();
    Nmet200->Write();
    NmetFilter->Write();
    Nlep0->Write();
    Njet1->Write();
    NcosmicVeto->Write();
    NEvents->Write();
    // outFile->Write();
    outFile->Close();
  }
}
