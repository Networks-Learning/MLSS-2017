#include "stdafx.h"
#include "mlss.h"

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nMLSS'17: Compute properties. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try

  // input network
  const TStr InputNetworkFileName = Env.GetIfArgPrefixStr("-i:", TStr("network.txt"), "Input network filename\n");

  // output name
  const TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "network", "Output prefix filename");

  // MLSS class
  TMLSS MLSS;

  // load network
  TFIn FIn(InputNetworkFileName);
  MLSS.LoadNetworkTxt(FIn);

  // compute property
  TIntPrV CntV;

  // outdegree distribution
  TSnap::GetOutDegCnt(MLSS.Network, CntV);
  TGnuPlot::PlotValV(CntV, TStr::Fmt("%s-out-deg", OutFNm.CStr()), "", "Out-degree", "N", gpsLog10XY, true, gpwPoints);
  printf("\nOut-degree distribution: %s\n", TStr::Fmt("%s-out-deg.eps", OutFNm.CStr()).CStr());

  // indegree distribution
  TSnap::GetInDegCnt(MLSS.Network, CntV);
  TGnuPlot::PlotValV(CntV, TStr::Fmt("%s-in-deg", OutFNm.CStr()), "", "In-degree", "N", gpsLog10XY, true, gpwPoints);
  printf("In-degree distribution: %s\n", TStr::Fmt("%s-in-deg.eps", OutFNm.CStr()).CStr());

  // weakly connected component size distribution
  TSnap::GetWccSzCnt(MLSS.Network, CntV);
  TGnuPlot::PlotValV(CntV, TStr::Fmt("%s-wcc", OutFNm.CStr()), "", "Component size", "N", gpsAuto, false, gpwPoints);
  printf("Weakly component size distribution: %s\n", TStr::Fmt("%s-wcc.eps", OutFNm.CStr()).CStr());

  // number of triangles
  printf("Number of triangles: %d\n", (int)TSnap::GetTriads(MLSS.Network));

  // clustering coefficient
  printf("Clustering coefficient: %f\n", TSnap::GetClustCf(MLSS.Network));

  // TODO: compute extra properties

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
