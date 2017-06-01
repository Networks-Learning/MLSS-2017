#include "stdafx.h"
#include "mlss.h"

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nMLSS'17: Convert from txt to gephi. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try

  // input network
  const TStr InputNetworkFileName = Env.GetIfArgPrefixStr("-i:", TStr("network.txt"), "Filename of the input network\n");

  // network name
  const TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "network", "Name of the output network");

  TMLSS MLSS;  // MLSS class

  // load network
  TFIn FIn(InputNetworkFileName);
  MLSS.LoadNetworkTxt(FIn);

  MLSS.SaveNetworkGephi(TStr::Fmt("%s.gexf", OutFNm.CStr()));

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
