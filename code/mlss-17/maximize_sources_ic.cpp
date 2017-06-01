#include "stdafx.h"
#include "mlss.h"

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nMLSS'17: Maximize influence by source set (indep cascade model). build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try

  // input network
  const TStr InputNetworkFileName = Env.GetIfArgPrefixStr("-i:", TStr("network.txt"), "Filename of the input network\n");

  // simulation type
  const int TOptimization = Env.GetIfArgPrefixInt("-t:", 0, "Optimization type\n\
  0:out-degree\n");

  // number of sources for random sources (TMode:0)
  const int NumSources = Env.GetIfArgPrefixFlt("-ns:", 1, "Number of random sources in the network (for TMode:0)\n");

  // number of cascades
  const int NCascades = Env.GetIfArgPrefixInt("-c:", 100, "Number of cascades\n");

  // network name
  const TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "cascades", "Cascades name");

  // if we use a network without edge weights, fix infection probability
  const double Prob = Env.GetIfArgPrefixFlt("-p:", -1.0, "Infection probability");

  TMLSS MLSS;  // MLSS class
  double AvCascSz = 0.0, StdError = 0.0; // average cascade size / standard error

  // load network
  TFIn FIn(InputNetworkFileName);
  if (Prob>0) { MLSS.LoadNetworkTxt(FIn, Prob); }
  else { MLSS.LoadNetworkTxt(FIn); }

  TInt::Rnd.Randomize();

  // find optimal sources
  if (TOptimization==0) {
	  TIntH NodeOutDegree;
	  for (TStrFltNEDNet::TNodeI NI = MLSS.Network->BegNI(); NI < MLSS.Network->EndNI(); NI++) {
		  NodeOutDegree.AddDat(NI.GetId()) = NI.GetOutDeg();
	  }
	  NodeOutDegree.SortByDat(false); // source nodes by descending out degree

	  for (int i=0; i<TMath::Mn(NumSources, NodeOutDegree.Len()); i++) {
		  MLSS.AddSource(NodeOutDegree.GetKey(i));
	  }
  }
  // TODO: add alternative methods to find optimal sources (TOptimization=1, TOptimization=2, etc..)
  // for independent cascade model
  else if (TOptimization==1) {
	  // other method to find optimal sources
  } else if (TOptimization==2) {
	  // another method to find optimal sources
  }
  // TODO: end

  // generate cascades for testing solution
  MLSS.CascH.Clr();
  for (int c=0; c<NCascades; c++) {
	  TCascade C;
	  MLSS.GenCascadeIC(C);

	  // add cascade to cascade set
	  MLSS.CascH.AddDat(MLSS.CascH.Len()) = C;

	  // cascade statistics udpate
	  AvCascSz += ((double)C.Len());
	  StdError += pow((double)C.Len(), 2.0);
  }
  AvCascSz /= (double)NCascades;
  StdError /= (double)NCascades;
  StdError -= AvCascSz;
  StdError = sqrt(StdError)/sqrt(NCascades);

  // optimal sources
  printf("Sources: %d (%s)", MLSS.SourcesV[0].Val, MLSS.GetNodeNm(MLSS.SourcesV[0]).CStr());
  for (int i=1;i<MLSS.SourcesV.Len(); i++) { printf(", %d (%s)", MLSS.SourcesV[i].Val, MLSS.GetNodeNm(MLSS.SourcesV[i]).CStr()); }
  printf("\n");

  // print out average influence
  printf("Average influence (cascade size): %f (std_error:%f)\n", AvCascSz, StdError);

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}

