// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey for the benefit of the
// BioDynaMo collaboration. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef ACTIN_H_
#define ACTIN_H_

#include "biodynamo.h"
#include "neuroscience/neuroscience.h"

namespace bdm {

struct SimParam : public ParamGroup {
  BDM_PARAM_GROUP_HEADER(SimParam, 1);
  
  double mempressure = 10;
  int ACTearly = 10000;
  double BeKaATP  = 11.6;
  double BeKaADPP = 3.4;
  double BeKaADP  = 2.9;
  double BeKdATP  = 1.4;
  double BeKdADPP = 0.2;
  double BeKdADP  = 5.4;
  double PeKaATP  = 1.3;
  double PeKaADPP = 0.11;
  double PeKaADP  = 0.09;
  double PeKdATP  = 0.9;
  double PeKdADPP = 0.02;
  double PeKdADP  = 0.25;

  double BeKaATD = (BeKaATP + (BeKaADPP + BeKaADP)/2)/2;
  double BeKdATD = (BeKdATP + (BeKdADPP + BeKdADP)/2)/2;
  double PeKaATD = (PeKaATP + (PeKaADPP + PeKaADP)/2)/2;
  double PeKdATD = (PeKdATP + (PeKdADPP + PeKdADP)/2)/2;

  double KaATD = round(((BeKaATD + PeKaATD)/2)*100)/100;
  double KdATD = round(((BeKdATD + PeKdATD)/2)*100)/100;
  double TKa = KaATD;
  double TKd = KdATD;

  double SPYneckXY = 300;
  double SPYheadZN = 1000;
  double SPYheadZS = 700;
  double SPYheadX  = 500;
  double SPYneckRad = SPYneckXY/2;
  double SPYheadRad = SPYheadX/2;
  double Vneck = Math::kPi * pow(SPYneckRad,2) * SPYheadZS;
  double Vhead = Math::kPi * pow(SPYheadRad,2) * (SPYheadZN-SPYheadZS);
  double SpyV = (Vneck+Vhead) * 1e-24;
  
  // Actin
  double MOL = 6.022e23; 
  double uM_Act = 300;        
  double GActinN0 = uM_Act / 1e6 * MOL * SpyV;
  double Gactin_uM = GActinN0 / SpyV *(1/MOL)*1e6;
  double FActinN     = 1;
  double GActinN     = ceil(GActinN0) - FActinN;
  double dt = 0.06;
  double fKa = TKa * Gactin_uM * dt;
  double fKd = TKd * dt;
  int NFact = 1;
  
  // Arp
  double ArpBR = 1;
  double Arp_Sc  = ArpBR * dt;
  int nStartFils = 1;
  double FArpN  = nStartFils;
  double GArpN0 = 9 / 1e6 * 6e23 * SpyV;
  double GArpN = ceil(GArpN0);
  double Arp_uM = GArpN / SpyV *(1/6e23)*1e6;

  double ArpAdd = 5;
  double ARPmax = 1000;

  // Coffilin
  double uM0_Cofilin = 1;  
  double N0_Cofilin = uM0_Cofilin / 1e6 * MOL * SpyV;    
  double Cofilin_uM = N0_Cofilin / SpyV *(1/MOL)*1e6;   

  double Ka_Cofilin     = 1e-3;      
  double Ka_Cofilin_LTP = 1e-3;     
  double Ka_Cofilin_LTD = 0.01;      

  double Factin_uM = FActinN / SpyV *(1/MOL)*1e6;    
  double KaCof = (dt * Ka_Cofilin * Factin_uM * Cofilin_uM);

  // Thymosin
  double THYM_ACT_uM  = 125;
  double THYM_uM      = 250;

  double THYM_ACT_N   = round(THYM_ACT_uM / 1e6 * MOL * SpyV);
  double THYM_N       = round(THYM_uM / 1e6 * MOL * SpyV);

  double Ka_Thymosin_LTP = 0.01;
  double Kd_Thymosin_LTP = 0.1;
  double Ka_Thymosin  = 0.08;
  double Kd_Thymosin  = 0.04;

  double tKa = Ka_Thymosin * THYM_uM * Gactin_uM * dt;
  double tKd = Kd_Thymosin * THYM_ACT_uM * dt;

  double THYM_D = 1e-3;

  // Constants
  int TOTAL_ACTIN     = GActinN + FActinN + THYM_ACT_N;
  int TOTAL_THYMOSIN  = THYM_N + THYM_ACT_N;
  int TOTAL_ARP       = GArpN + FArpN;
  int TOTAL_COFILIN   = N0_Cofilin;
  int timeStep = 1;
  int ThymTimeStart = 12000;
  int LTPon = 40000;
  int LTPoff = 43000;
};

struct Polymerization : public Behavior {
  BDM_BEHAVIOR_HEADER(Polymerization, Behavior, 1);
  Polymerization() { AlwaysCopyToNew(); }
  virtual ~Polymerization() {}
  
  double fKa;

  void fKaUpdate(){
    auto* sim = Simulation::GetActive();
    auto* P = const_cast<Param*>(sim->GetParam())->Get<SimParam>();

    P->Factin_uM = P->FActinN / P->SpyV * (1/P->MOL) * 1e6;
	  P->Gactin_uM = P->GActinN / P->SpyV * (1/P->MOL) * 1e6;
	  fKa = P->TKa * P->Gactin_uM * P->dt;
  }

  bool TipOk(NeuriteElement* filament) {
    auto P = const_cast<Param*>(Simulation::GetActive()->GetParam())->Get<SimParam>();

    auto Tip_xyz = filament->DistalEnd(); 

    auto XYtipLoc = sqrt(pow(Tip_xyz[0],2) + pow(Tip_xyz[1],2));

    auto ZtipInHead = Tip_xyz[2] >= P->SPYheadZS;
    auto ZtipInNeck = Tip_xyz[2] < P->SPYheadZS;

    auto XYneckOut = XYtipLoc > P->SPYneckRad;
    auto XYheadOut = XYtipLoc > P->SPYheadRad;
    auto ZtopOut = Tip_xyz[2] > P->SPYheadZN;
    auto ZbotOut = Tip_xyz[2] < 0;
    
    auto TipOut = ((XYneckOut && ZtipInNeck)+(XYheadOut && ZtipInHead)+ZtopOut+ZbotOut) > 0;
	   
    return !TipOut;
  }

  double LO(NeuriteElement* filament) {
    auto P = const_cast<Param*>(Simulation::GetActive()->GetParam())->Get<SimParam>();

    auto Tip_xyz = filament->DistalEnd(); 

    auto XYtipLoc = sqrt(pow(Tip_xyz[0],2) + pow(Tip_xyz[1],2));

    auto ZtipInHead = Tip_xyz[2] >= P->SPYheadZS;
    auto ZtipInNeck = Tip_xyz[2] < P->SPYheadZS;

    auto LngXYneckOut = (XYtipLoc - P->SPYneckRad) * ZtipInNeck;  
	  auto LngXYheadOut = (XYtipLoc - P->SPYheadRad) * ZtipInHead;
	  auto LngZtopOut = (Tip_xyz[2] - P->SPYheadZN) * ZtipInHead;
    
    LngXYneckOut = LngXYneckOut<0 ? 0 : LngXYneckOut;
    LngXYheadOut = LngXYheadOut<0 ? 0 : LngXYheadOut;
    LngZtopOut = LngZtopOut<0 ? 0 : LngZtopOut;
    
    auto LngOut = LngXYneckOut + LngXYheadOut + LngZtopOut;
    
    return LngOut / 100 * P->mempressure;
  }

  void Run(Agent* agent) override {
    //std::cout << "Poly" << std::endl;
    auto* sim = Simulation::GetActive();
    auto random = sim->GetRandom()->Uniform(0,1);
    auto* P = const_cast<Param*>(sim->GetParam())->Get<SimParam>();
    
    fKaUpdate();
    auto* filament = bdm_static_cast<NeuriteElement*>(agent);
    auto ZbotOut = filament->DistalEnd()[2] < 0;

    if (filament->GetDaughterLeft()!= nullptr) return;

    if (sim->GetScheduler()->GetSimulatedSteps() < P->ACTearly){
      //std::cout << "early" << std::endl;
      if (fKa > random && TipOk(filament) && (P->GActinN > 1)){
        if (fKa > P->FArpN){
          filament->ElongateTerminalEnd(2*2.71, filament->GetSpringAxis());
          P->FActinN += 2;
          P->GActinN -= 2;
        } else {
          filament->ElongateTerminalEnd(2.71, filament->GetSpringAxis());
          P->FActinN += 1;
          P->GActinN -= 1;
        }
      }
    } else {
      if ((fKa - LO(filament) > random) && !ZbotOut && (P->GActinN > 1) && TipOk(filament)) {
        if (fKa > P->FArpN){
          filament->ElongateTerminalEnd(2*2.71, filament->GetSpringAxis());
          P->FActinN += 2;
          P->GActinN -= 2;
        } else {
          filament->ElongateTerminalEnd(2.71, filament->GetSpringAxis());
          P->FActinN += 1;
          P->GActinN -= 1;
        }
      } 
    }
  }

 private:
  bool init_ = false;
  DiffusionGrid* dg_guide_ = nullptr;
};

struct Branching : public Behavior {
  BDM_BEHAVIOR_HEADER(Branching, Behavior, 1);
  Branching() { AlwaysCopyToNew(); }
  virtual ~Branching() {}
  double ArpBR;

  void updateArpBR(NeuriteElement* filament){
    auto* sim = Simulation::GetActive();
    auto* P = const_cast<Param*>(sim->GetParam())->Get<SimParam>();

    P->Arp_uM = P->GArpN / P->SpyV * (1/P->MOL) * 1e6;
    int FactinUnits = round(filament->GetLength()/2.71);
    
    double Fact_Branch_uM = FactinUnits / P->SpyV * (1/P->MOL) * 1e6;
    ArpBR = P->Arp_Sc * (P->Arp_uM / 1000 * Fact_Branch_uM);
  }
  
  Double3 GetNewDirection(NeuriteElement* filament){
    auto old_direction = filament->GetSpringAxis();
    auto r = filament->GetLength();
    auto theta = acos(old_direction[2]/r) + Math::ToRadian(70);
    auto phi = atan2(old_direction[1],old_direction[0]);

    auto x = r * cos(phi) * sin(theta); 
    auto y = r * sin(phi) * sin(theta);
    auto z = r * cos(theta);

    Double3 new_direction = {x,y,z};
    
    auto random_rotation_angle = round(Simulation::GetActive()->GetRandom()->Uniform(360));

    return Math::RotAroundAxis(new_direction,Math::ToRadian(random_rotation_angle),filament->GetSpringAxis());
  }

  void Run(Agent* agent) override {
    //std::cout << "Branch" << std::endl;
    auto* sim = Simulation::GetActive();
    auto random = sim->GetRandom();
    auto* P = const_cast<Param*>(sim->GetParam())->Get<SimParam>();
    auto* filament = bdm_static_cast<NeuriteElement*>(agent);
    if (filament->GetLength() < 2.71) return;
    updateArpBR(filament);
    auto branch_direction = GetNewDirection(filament);
    
    if (ArpBR > random->Uniform(0,1) && P->GActinN > P->ArpAdd) {
      /*std::cout << filament->GetLength() << std::endl;
      std::cout << filament->GetActualLength() << std::endl;*/
      auto event = bdm::neuroscience::NeuriteBranchingEvent(random->Uniform(0,1),P->ArpAdd * 2.71, filament->GetDiameter(), branch_direction);
      filament->CreateNewAgents(event, {filament,filament});
      auto branch = bdm_static_cast<NeuriteElement*>(event.new_agents[1]);
      P->FArpN+=1;
      P->GArpN-=1;
      P->GActinN-=P->ArpAdd;
    }
  }

 private:
  bool init_ = false;
  DiffusionGrid* dg_guide_ = nullptr;
};

struct Depolymerization : public Behavior {
  BDM_BEHAVIOR_HEADER(Depolymerization, Behavior, 1);
  Depolymerization() { AlwaysCopyToNew(); }
  virtual ~Depolymerization() {}

  void Run(Agent* agent) override {
    //std::cout << "Depoly" << std::endl;
    auto* sim = Simulation::GetActive();
    auto random = sim->GetRandom();
    auto* P = const_cast<Param*>(sim->GetParam())->Get<SimParam>();
    auto* filament = bdm_static_cast<NeuriteElement*>(agent);
    
    if (filament->GetDaughterLeft() != nullptr) return;

    if (P->fKd > random->Uniform(0,1)) {
      if (filament->GetLength() < 5.42) {
        auto mother = dynamic_cast<NeuriteElement*>(filament->GetMother().Get());
        auto ptr = filament->GetAgentPtr<NeuriteElement>();
        auto left_daughter = mother->GetDaughterLeft();
        if (ptr == left_daughter){
          mother->RemoveDaughter(left_daughter);
          filament->RemoveFromSimulation();
          P->GArpN+=1;
          P->FArpN-=1;
        }
      }else {
        filament->RetractTerminalEnd(2.71);
      }
      P->FActinN-=1;
      P->GActinN+=1;
    }
  }

 private:
  bool init_ = false;
  DiffusionGrid* dg_guide_ = nullptr;
};

struct Severing : public Behavior {
  BDM_BEHAVIOR_HEADER(Severing, Behavior, 1);
  Severing() { AlwaysCopyToNew(); }
  virtual ~Severing() {}

  double KaCof;

  void updateKaCof(){
    auto* sim = Simulation::GetActive();
    auto* P = const_cast<Param*>(sim->GetParam())->Get<SimParam>();
  
    int Nfillaments = P->FArpN;
    KaCof = (P->dt/Nfillaments * P->Ka_Cofilin * P->Factin_uM * P->Cofilin_uM);
  }

  void Run(Agent* agent) override {
    //std::cout << "Severing" << std::endl;
    auto* sim = Simulation::GetActive();
    auto random = sim->GetRandom();
    auto* P = const_cast<Param*>(sim->GetParam())->Get<SimParam>();
    auto* filament = bdm_static_cast<NeuriteElement*>(agent);

    updateKaCof();

    if (filament->GetDaughterLeft()!= nullptr) return;

    if (KaCof > random->Uniform(0,1) && filament->GetLength() > 10 * 2.71) {
      double n_severed_actin = round(random->Uniform(1,filament->GetLength()/2.71));
      filament->RetractTerminalEnd(n_severed_actin*2.71);
      P->FActinN-=n_severed_actin;
      P->GActinN+=n_severed_actin;
    };
  }

 private:
  bool init_ = false;
  DiffusionGrid* dg_guide_ = nullptr;
};

struct Thymosin : public Behavior {
  BDM_BEHAVIOR_HEADER(Thymosin, Behavior, 1);
  Thymosin() { AlwaysCopyToNew(); }
  virtual ~Thymosin() {}
  double tKa;
  double tKd;

  void updateThymosinRates(){
    auto* sim = Simulation::GetActive();
    auto* P = const_cast<Param*>(sim->GetParam())->Get<SimParam>();
    int Nfillaments = P->FArpN;

    P->THYM_uM     = P->THYM_N     / P->SpyV * (1/P->MOL) * 1e6;
    P->THYM_ACT_uM = P->THYM_ACT_N / P->SpyV * (1/P->MOL) * 1e6;

    tKa = P->Ka_Thymosin * P->THYM_uM * P->Gactin_uM * P->dt/Nfillaments;
    tKd = P->Kd_Thymosin * P->THYM_ACT_uM * P->dt/Nfillaments;
  }

  void Run(Agent* agent) override {
    auto* sim = Simulation::GetActive();
    auto random = sim->GetRandom();
    auto* P = const_cast<Param*>(sim->GetParam())->Get<SimParam>();
    auto* filament = bdm_static_cast<NeuriteElement*>(agent);
    int time_step = sim->GetScheduler()->GetSimulatedSteps();
    
    updateThymosinRates();

    if (time_step > P->ThymTimeStart){
      if ((tKa > random->Uniform(0,1)) && (P->GActinN>1) && (P->THYM_N>1)){
          P->THYM_N = P->THYM_N - 1;
          P->THYM_ACT_N = P->THYM_ACT_N + 1;
          P->GActinN = P->GActinN - 1;
      }
      if ((tKd > random->Uniform(0,1)) && (P->THYM_ACT_N>1)){
          P->THYM_N = P->THYM_N + 1;
          P->THYM_ACT_N = P->THYM_ACT_N - 1;
          P->GActinN = P->GActinN + 1;
      }
    }
    
}

 private:
  bool init_ = false;
  DiffusionGrid* dg_guide_ = nullptr;
};

struct Time : public Behavior {
  BDM_BEHAVIOR_HEADER(Time, Behavior, 1);
  Time() { NeverCopyToNew(); }
  virtual ~Time() {}
  double PctTeq;
  double PctTAeq;
 
  void Run(Agent* agent) override {
    auto* sim = Simulation::GetActive();
    auto* P = const_cast<Param*>(sim->GetParam())->Get<SimParam>();
    int time_step = sim->GetScheduler()->GetSimulatedSteps();


    if (time_step == P->LTPon){
      int TAeq = P->THYM_ACT_N;
      int Teq = P->THYM_N;
      double PctTeq = Teq / (TAeq + Teq);
      double PctTAeq = TAeq / (TAeq + Teq);
      P->Ka_Thymosin = P->Ka_Thymosin_LTP;
      P->Kd_Thymosin = P->Kd_Thymosin_LTP;
    }
    
    if (time_step > P->LTPon && time_step < P->LTPoff){
      int NrepT = round(P->THYM_N * P->THYM_D);
      int NrepTA = round(P->THYM_ACT_N * P->THYM_D);
      int NrepTAT = NrepT + NrepTA;

      P->THYM_N += -NrepT + round(PctTeq*NrepTAT);
      P->THYM_ACT_N += -NrepTA + round(PctTAeq*NrepTAT);
      
      P->TOTAL_ACTIN = P->GActinN + P->FActinN + P->THYM_ACT_N;
    }

    if (time_step == P->LTPoff){
      P->Ka_Thymosin = 0.08;
      P->Kd_Thymosin = 0.04;
    }

    if ( time_step % 1000 == 0){
      std::cout << "Time step: "<< time_step
      << "  GActin: "<< P->GActinN
      << "  FActin: "<< P->FActinN
      << "  GArpN: "<< P->GArpN
      << "  FArpN: "<< P->FArpN
      << "  TA: "<< P->THYM_ACT_N
      << "  T: "<< P->THYM_N << std::endl;
    }
  }

 private:
  bool init_ = false;
  DiffusionGrid* dg_guide_ = nullptr;
};


inline int Simulate(int argc, const char** argv) {
  neuroscience::InitModule();
  Param::RegisterParamGroup(new SimParam());
  
  Simulation simulation(argc, argv);

  auto* rm = simulation.GetResourceManager();
  auto* soma = new neuroscience::NeuronSoma();
  soma->SetDiameter(10);
  soma->AddBehavior(new Time());
  auto*neurite = soma->ExtendNewNeurite({0,0,1});
  neurite->AddBehavior(new Polymerization());
  neurite->AddBehavior(new Branching());
  neurite->AddBehavior(new Depolymerization());
  neurite->AddBehavior(new Severing());
  neurite->AddBehavior(new Thymosin());
  rm->AddAgent(soma);
  auto* P = const_cast<Param*>(simulation.GetParam())->Get<SimParam>();
  
  // Run simulation for one timestep
  simulation.GetScheduler()->Simulate(80000);
  
  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

}  // namespace bdm

#endif  // ACTIN_H_
