#include "stdafx.h"
#include "mlss.h"

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nMLSS'17: Maximize influence by source set (continuous time indep cascade model). build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try

  // input network
  const TStr InputNetworkFileName = Env.GetIfArgPrefixStr("-i:", TStr("network"), "Name of the input network (default:network)\n");

  // simulation type
  const int TOptimization = Env.GetIfArgPrefixInt("-t:", 0, "Optimization type\n\
  0:edges to lowest out-degree\n");

  // number of sources for random sources (TMode:0)
  const int NumSources = Env.GetIfArgPrefixFlt("-ns:", 1, "Number of random sources in the network (for TMode:0)\n");

  // observation window
  const double Window = Env.GetIfArgPrefixFlt("-w:", 1.0, "Time window\n");

  // pairwise transmission model
  const TModel Model = (TModel)Env.GetIfArgPrefixInt("-m:", 0, "Transmission model\n0:exponential, 1:power law, 2:rayleigh, 3:weibull (default:0)\n");

  // delta for power law transmission model
  const double Delta = Env.GetIfArgPrefixFlt("-d:", 0.01, "Delta for power law transmission model (-m:1)\n");

  // number of cascades
  const int NCascades = Env.GetIfArgPrefixInt("-c:", 100, "Number of cascades\n");

  // network name
  const TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "cascades", "Cascades name");

  TMLSS MLSS;  // MLSS class
  double AvCascSz = 0.0, StdError = 0.0; // average cascade size / standard error
  TIntPrV OptimalEdges;

  // load network
  TFIn FIn(InputNetworkFileName);
  MLSS.LoadNetworkTxt(FIn);

  // set tx model
  MLSS.SetModel(Model);

  // set time window
  MLSS.SetWindow(Window);

  // delta for power law transmission rate
  MLSS.SetDelta(Delta);

  TInt::Rnd.Randomize();

  // find optimal sources
  if (TOptimization==0) {
	  TIntH NodeOutDegree;
	  for (TStrFltNEDNet::TNodeI NI = MLSS.Network->BegNI(); NI < MLSS.Network->EndNI(); NI++) {
		  NodeOutDegree.AddDat(NI.GetId()) = NI.GetOutDeg();
	  }
	  NodeOutDegree.SortByDat(false); // source nodes by ascending out degree

	  for (int i=0; i<TMath::Mn(NumSources, NodeOutDegree.Len()); i++) {
		  int DstNode =  MLSS.Network->GetRndNId();
		  while (DstNode==NodeOutDegree.GetKey(i)) { DstNode = MLSS.Network->GetRndNId(); } // avoid self loops

		  MLSS.Network->AddEdge(NodeOutDegree.GetKey(i), DstNode);
		  OptimalEdges.Add(TIntPr(NodeOutDegree.GetKey(i), DstNode));
	  }
  }
  // TODO: add alternative methods to find optimal edges to add (TOptimization=1, TOptimization=2, etc..)
  // for continuous time independent cascade model
  else if (TOptimization==1) {
	  // other method to find optimal sources
  } else if (TOptimization==2) {
	  // another method to find optimal sources
  }
  // TODO: end

  // generate cascades
  MLSS.CascH.Clr();
  for (int c=0; c<NCascades; c++) {
	  // generate random source set
	  MLSS.SetSources(NumSources);

	  TCascade C;
	  MLSS.GenCascadeCTIC(C);

	  // add cascade to cascade set
	  MLSS.CascH.AddDat(MLSS.CascH.Len()) = C;

	  // cascade statistics udpate
	  AvCascSz += ((double)C.Len());
  }
  AvCascSz /= (double)NCascades;
  StdError /= (double)NCascades;
  StdError -= AvCascSz;
  StdError = sqrt(StdError)/sqrt(NCascades);

  // print out average influence
  printf("Average influence (cascade size): %f (std_error:%f)\n", AvCascSz, StdError);

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}

