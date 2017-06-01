#include "stdafx.h"
#include "mlss.h"

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nMLSS'17: Generate random networks. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try

  // network type: kronecker or forest fire
  const int TNetwork = Env.GetIfArgPrefixInt("-t:", -1, "Network to generate\n0:kronecker, 1:forest fire\n");

  // network params for kronecker/forest-fire
  // kronecker format, e.g., 0.5 0.5; 0.5 0.5
  // forest fire format, e.g., 1;0.2;0.17;1;0;0
  const TStr NetworkParams = Env.GetIfArgPrefixStr("-g:", TStr("0.5 0.5; 0.5 0.5"), "Parameters for the network\nKronecker format: 0.5 0.5; 0.5 0.5\nForest fire format: 1;0.2;0.17;1;0;0\n");

  // nodes, edges
  const int NNodes = Env.GetIfArgPrefixInt("-n:", 512, "Number of nodes\n");
  const int NEdges = Env.GetIfArgPrefixInt("-e:", 1024, "Number of edges (for -t:0)\n");

  const int MNetwork = Env.GetIfArgPrefixInt("-m:", 0, "Network type\n0:discrete time, 1:continuous time\n");

  // edge min/max weight
  const double MinParameter = Env.GetIfArgPrefixFlt("-la:", 0, "Minimum edge weight\n");
  const double MaxParameter = Env.GetIfArgPrefixFlt("-ua:", 1.0, "Maximum edge weight\n");

  // network name
  const TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "random-network", "Network name");

  // flags for saving
  const int SaveTxt = Env.GetIfArgPrefixInt("-st:", 1, "Save network txt format (0:no, 1:yes)");
  const int SaveGephi = Env.GetIfArgPrefixInt("-sg:", 0, "Save network Gephi format (0:no, 1:yes)");
  const int SaveMatlab = Env.GetIfArgPrefixInt("-sm:", 0, "Save network Matlab format (0:no, 1:yes)");

  // verbose
  const int Verbose = Env.GetIfArgPrefixInt("-v:", 0, "Print network (0:no, 1:yes)");

  // by default, print help only
  if (TNetwork==-1) { return (-1); }

  // MLSS class
  TMLSS MLSS;

  // generate random network
  MLSS.GenerateNetwork(TNetwork, NNodes, NEdges, NetworkParams, Verbose);

  // generate edge weights
  if (MNetwork==0) { MLSS.GenerateBetas(MinParameter, MaxParameter); }
  else { MLSS.GenerateAlphas(MinParameter, MaxParameter); }

  // print network number of nodes and edges
  printf("\nDirected network has been generated: %d nodes, %d edges!\n", MLSS.Network->GetNodes(), MLSS.Network->GetEdges());

  // save network to text format
  if (SaveTxt==1) {
	  MLSS.SaveNetworkTxt(TStr::Fmt("%s.txt", OutFNm.CStr()));
	  printf("Network in txt format: %s\n", TStr::Fmt("%s.txt", OutFNm.CStr()).CStr());
  }

  // save network to gephi format
  if (SaveGephi==1) {
	  MLSS.SaveNetworkGephi(TStr::Fmt("%s.gexf", OutFNm.CStr()));
	  printf("Network in Gephi format: %s\n", TStr::Fmt("%s.gexf", OutFNm.CStr()).CStr());
  }

  // save network to matlab format
  if (SaveMatlab==1) {
  	  MLSS.SaveNetworkMatlab(TStr::Fmt("%s.mat", OutFNm.CStr()));
  	  printf("Network in ascii Matlab format: %s\n", TStr::Fmt("%s.mat", OutFNm.CStr()).CStr());
  }

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
