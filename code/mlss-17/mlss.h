#ifndef snap_mlss_h
#define snap_mlss_h

#include "Snap.h"

typedef TNodeEDatNet<TFlt, TFlt> TFltNEDNet;
typedef TPt<TFltNEDNet> PFltNEDNet;

typedef TNodeEDatNet<TStr, TFlt> TStrFltNEDNet;
typedef TPt<TStrFltNEDNet> PStrFltNEDNet;

// pairwise transmission models
typedef enum {
	EXP, // exponential
	POW, // powerlaw
	RAY, // rayleigh
} TModel;

// Hit info (timestamp, first parent) about a node in a cascade
class THitInfo {
public:
  TInt NId, Parent;
  TFlt Tm;
public:
  THitInfo(const int& NodeId=-1, const double& HitTime=0) : NId(NodeId), Parent(-1), Tm(HitTime) { }
  THitInfo(const int& NodeId, const int& Prnt, const double& HitTime) : NId(NodeId), Parent(Prnt), Tm(HitTime) { }
  THitInfo(TSIn& SIn) : NId(SIn), Parent(SIn), Tm(SIn) { }
  void Save(TSOut& SOut) const { NId.Save(SOut); Parent.Save(SOut); Tm.Save(SOut); }
  bool IsParent() const { return (NId==-1); }
  bool operator < (const THitInfo& Hit) const {
    return Tm < Hit.Tm; }
};

// Cascade
class TCascade {
public:
  TInt CId; // cascade id
  THash<TInt, THitInfo> NIdHitH; // infection times & first parents
  PFltNEDNet Network; // graph cascade (with infection times in nodes, transmission times in edges)
public:
  TCascade() : CId(0), NIdHitH() { Network = TFltNEDNet::New(); }
  TCascade(const int &cid) : CId(cid), NIdHitH() { Network = TFltNEDNet::New(); }
  TCascade(TSIn& SIn) : CId(SIn), NIdHitH(SIn), Network(SIn) { }
  void Save(TSOut& SOut) const  { CId.Save(SOut); NIdHitH.Save(SOut); Network.Save(SOut); }
  void Clr() { NIdHitH.Clr(); Network->Clr(); }
  int GetId() { return CId; }
  int Len() const { return NIdHitH.Len(); }
  int GetNode(const int& i) const { return NIdHitH.GetKey(i); }
  THash<TInt, THitInfo>::TIter BegI() const { return NIdHitH.BegI(); }
  THash<TInt, THitInfo>::TIter EndI() const { return NIdHitH.EndI(); }
  double GetTm(const int& NId) const { return NIdHitH.GetDat(NId).Tm; }
  double FirstTm() const { return NIdHitH[0].Tm.Val; }
  double LastTm() const { return NIdHitH[NIdHitH.Len()-1].Tm.Val; }
  void Add(const int& NId, const double& HitTm) { NIdHitH.AddDat(NId, THitInfo(NId, HitTm)); }
  void Add(const int& NId, const int& Prnt, const double& HitTm) { NIdHitH.AddDat(NId, THitInfo(NId, Prnt, HitTm)); }
  void Del(const int& NId) { NIdHitH.DelKey(NId); }
  bool IsNode(const int& NId) const { return NIdHitH.IsKey(NId); }
  void Sort() { NIdHitH.SortByDat(true); }
};

// Node info (name and number of cascades)
class TNodeInfo {
public:
  TStr Name;
  TInt Vol;
public:
  TNodeInfo() { }
  TNodeInfo(const TStr& NodeNm, const int& Volume) : Name(NodeNm), Vol(Volume) { }
  TNodeInfo(TSIn& SIn) : Name(SIn), Vol(SIn) { }
  void Save(TSOut& SOut) const { Name.Save(SOut); Vol.Save(SOut); }
};


// MLSS class
class TMLSS {
public:
  THash<TInt, TCascade> CascH; // cascades

  PStrFltNEDNet Network; // network (beta or alpha in edges)
  THash<TInt, TNodeInfo> NodeNmH; // node info

  TIntV SourcesV; // source set

  TFlt Window; // time horizon for continuous time independent cascade model
  TFlt Delta;   // delta for power-law

  TModel Model; // transmission model for CTIC

public:
  TMLSS( ) { Network = TStrFltNEDNet::New(); }
  TMLSS(TSIn& SIn) : CascH(SIn), Network(SIn), NodeNmH(SIn), SourcesV(SIn), Window(SIn), Delta(SIn) { }
  void Save(TSOut& SOut) const { CascH.Save(SOut); Network.Save(SOut); NodeNmH.Save(SOut); SourcesV.Save(SOut); Window.Save(SOut); Delta.Save(SOut); }

  void LoadNetworkTxt(TSIn& SIn, const double& Prob=-1);

  void GenerateNetwork(const int& TNetwork, const int& NNodes, const int& NEdges, const TStr& NetworkParams, const bool& verbose);
  void GenerateAlphas(const double& MinAlpha, const double& MaxAlpha);
  void GenerateBetas(const double& MinBeta, const double& MaxBeta);
  void ConvertAlphasToBetas(const double& T);

  int GetNodes() { return Network->GetNodes(); }
  void AddNodeNm(const int& NId, const TNodeInfo& Info) { NodeNmH.AddDat(NId, Info); }
  TStr GetNodeNm(const int& NId) const { return NodeNmH.GetDat(NId).Name; }
  TNodeInfo GetNodeInfo(const int& NId) const { return NodeNmH.GetDat(NId); }
  bool IsNodeNm(const int& NId) const { return NodeNmH.IsKey(NId); }
  bool IsSource(const int& NId) const { return SourcesV.IsIn(NId); }

  void GenCascadeIC(TCascade& C, const bool& verbose=false);
  void GenCascadeCTIC(TCascade& C, const bool& verbose=false);

  TCascade & GetCasc(int c) { return CascH[c]; }
  int GetCascs() { return CascH.Len(); }

  void SetModel(const int& m) { Model = (TModel)m; }
  void SetWindow(const double& window) { Window = window; }
  void SetDelta(const double& delta) { Delta = delta; }

  void SetSources(const TStrV& SourcesStrV);
  void SetSources(const int& NSources);
  void AddSource(const int& NId) { SourcesV.AddUnique(NId); }
  void DelSources() { SourcesV.Clr(); }
  void DelSource(const int& NId) { if (SourcesV.IsIn(NId)) SourcesV.DelIfIn(NId); }

  void SaveCascadesTxt(const TStr& OutFNm);
  void SaveCascadesGephi(const TStr& OutFNm, const int&NCascades, const bool& PrintSources=true, const double& TimeSep=1);
  void SaveCascadesMatlab(const TStr& OutFNm);

  void SaveNetworkTxt(const TStr& OutFNm);
  void SaveNetworkGephi(const TStr& OutFNm);
  void SaveNetworkMatlab(const TStr& OutFNm);
};

#endif
