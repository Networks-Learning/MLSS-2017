#include "stdafx.h"
#include "mlss.h"

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nMLSS'17: Simulate independent cascade model. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try

  // input network
  const TStr InputNetworkFileName = Env.GetIfArgPrefixStr("-i:", TStr("network.txt"), "Filename of the input network\n");

  // simulation type
  const int TMode = Env.GetIfArgPrefixInt("-t:", -1, "Source types\n\
      	  0:a fixed source set, 1:a single random source set for all cascades, 2:a random source set for each cascade (default:0)\n");

  // number of sources for random sources (TMode:0)
  const int NumSources = Env.GetIfArgPrefixFlt("-ns:", 1, "Number of random sources in the network (for TMode:1 or TMode:2)\n");

  // sources (TMod:0)
  const TStr Sources = Env.GetIfArgPrefixStr("-s:", TStr("1"), "Sources (source_id1;source_id2;...) (for TMode:0)\n");

  // number of cascades
  const int NCascades = Env.GetIfArgPrefixInt("-c:", 100, "Number of cascades\n");

  // number of cascades to save to gephi format
  const int NCascadeGephi = Env.GetIfArgPrefixInt("-cg:", 100, "Number of cascades to save to Gephi format (for -sg:1)\n");

  // flags for saving
  const int SaveTxt = Env.GetIfArgPrefixInt("-st:", 1, "Save cascades txt format (0:no, 1:yes)");
  const int SaveGephi = Env.GetIfArgPrefixInt("-sg:", 0, "Save cascades Gephi format (0:no, 1:yes)");
  const int SaveMatlab = Env.GetIfArgPrefixInt("-sm:", 0, "Save cascades Matlab format (0:no, 1:yes)");

  // network name
  const TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "cascades", "Cascades name");

  // verbose
  const int Verbose = Env.GetIfArgPrefixInt("-v:", 0, "Print cascade generation (0:no, 1:yes)");

  // if we use a network without edge weights, fix infection probability
  const double Prob = Env.GetIfArgPrefixFlt("-p:", -1.0, "Infection probability");

  printf("\n");

  // by default, print help only
  if (TMode==-1) { return (-1); }

  TMLSS MLSS;  // MLSS class
  double AvCascSz = 0.0, StdError = 0.0; // average cascade size
  TIntH CascSzDistr, NumNodeInf, NumNodeInfDistr; // cascade statistics

  // load network
  TFIn FIn(InputNetworkFileName);
  if (Prob>0) { MLSS.LoadNetworkTxt(FIn, Prob); }
  else { MLSS.LoadNetworkTxt(FIn); }

  TInt::Rnd.Randomize();

  // fill source set
  TStrV SourcesStrV;
  switch (TMode) {
  	  case 0:
  		  Sources.SplitOnAllCh(';', SourcesStrV);
  		  MLSS.SetSources(SourcesStrV);
  	  	  break;
  	  case 1:
  	  case 2:
  		  MLSS.SetSources(NumSources);
  		  break;
  	  default:
  		  printf("Unknown source type!\n");
  		  return (-1);
  		  break;
  }

  printf("\nSources: ");
  for (int i=0; i<MLSS.SourcesV.Len(); i++) { printf("%d ", MLSS.SourcesV[i].Val); }
  printf("\n");

  for (int c=0; c<NCascades; c++) {
	  if (Verbose==1) { printf("\nCascade %d\n", c); }

	  // generate cascade
	  TCascade C;
	  MLSS.GenCascadeIC(C, Verbose>0);

	  // add cascade to cascade set
	  MLSS.CascH.AddDat(MLSS.CascH.Len()) = C;

	  // cascade statistics udpate
	  AvCascSz += ((double)C.Len());
	  StdError += pow((double)C.Len(), 2.0);
	  if (!CascSzDistr.IsKey(C.Len())) { CascSzDistr.AddDat(C.Len()) = 0; }
	  CascSzDistr.GetDat(C.Len())++;
	  for (THash<TInt, THitInfo>::TIter NI = C.BegI(); NI < C.EndI(); NI++) {
		  if (MLSS.IsSource(NI.GetKey())) { continue; } // do not count infection of sources
		  if (!NumNodeInf.IsKey(NI.GetKey())) { NumNodeInf.AddDat(NI.GetKey()) = 0; }
		  NumNodeInf.GetDat(NI.GetKey())++;
	  }

	  // if TMode==2, generate new random source
	  if (TMode==2) {
		  MLSS.SetSources(NumSources);
		  printf("Sources: ");
		  for (int i=0; i<MLSS.SourcesV.Len(); i++) { printf("%d ", MLSS.SourcesV[i].Val); }
		  printf("\n");
	  }
  }

  // statistics
  AvCascSz /= (double)NCascades;
  StdError /= (double)NCascades;
  StdError -= AvCascSz;
  StdError = sqrt(StdError)/sqrt(NCascades);
  for (int i=0; i<NumNodeInf.Len(); i++) {
	  if (!NumNodeInfDistr.IsKey(NumNodeInf[i])) { NumNodeInfDistr.AddDat(NumNodeInf[i]) = 0; }
	  NumNodeInfDistr.GetDat(NumNodeInf[i])++;
  }
  CascSzDistr.SortByKey(true);
  NumNodeInfDistr.SortByKey(true);

  // print out average influence
  printf("Average influence (cascade size): %f (std_error:%f)\n", AvCascSz, StdError);

  // plot cascade size distribution
  TGnuPlot::PlotValCntH(CascSzDistr, TStr::Fmt("%s-cascsz", OutFNm.CStr()), "", "Cascade size", "N", gpsLog10XY, true, gpwPoints, false, false);

  // plot node infection distribution
  TGnuPlot::PlotValCntH(NumNodeInfDistr, TStr::Fmt("%s-nodeinf", OutFNm.CStr()), "", "Number of infections", "N", gpsLog10XY, false, gpwPoints, false, false);

  // save cascade to txt format
  if (SaveTxt==1) {
	  MLSS.SaveCascadesTxt(TStr::Fmt("%s.txt", OutFNm.CStr()));
	  printf("Cascades in txt format: %s\n", TStr::Fmt("%s.txt", OutFNm.CStr()).CStr());
  }

  // save cascades to gephi format
  if (SaveGephi==1) {
	  MLSS.SaveCascadesGephi(TStr::Fmt("%s.gexf", OutFNm.CStr()), NCascadeGephi, (TMode!=2));
	  printf("Cascades (%d) in Gephi format: %s\n", NCascadeGephi, TStr::Fmt("%s.gexf", OutFNm.CStr()).CStr());
  }

  // save cascades to ascii matlab format
  if (SaveMatlab==1) {
	  MLSS.SaveCascadesMatlab(TStr::Fmt("%s.mat", OutFNm.CStr()));
	  printf("Cascades in ASCII Matlab format: %s\n", TStr::Fmt("%s.mat", OutFNm.CStr()).CStr());
  }

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
