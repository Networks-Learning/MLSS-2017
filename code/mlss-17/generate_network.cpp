#include "stdafx.h"
#include "mlss.h"

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nMLSS'17: Generate network. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try

  // directed graph (smart pointer to TNGraph object)
  PNGraph PGraph = TNGraph::New(); TIntStrH GraphNodeLabels;

  // TODO: modify code to create double star network on the board
  for (int i=0; i<10; i++) { PGraph->AddNode(i); GraphNodeLabels.AddDat(i) = TStr::Fmt("%d", i); }
  for (int i=1; i<10; i++) { PGraph->AddEdge(i, 0); } // star
  // TODO: END

  printf("\nDirected graph (%d nodes, %d edges)\n", PGraph->GetNodes(), PGraph->GetEdges());

  // print node id's
  for (TNGraph::TNodeI NI = PGraph->BegNI(); NI < PGraph->EndNI(); NI++) {
	  printf("Node %d\n", NI.GetId());
  }

  // print edges
  for (TNGraph::TEdgeI EI = PGraph->BegEI(); EI < PGraph->EndEI(); EI++) {
	  printf("Edge %d -> %d\n", EI.GetSrcNId(), EI.GetDstNId());
  }

  // generate graphviz
  TSnap::DrawGViz(PGraph, gvlDot, TStr("graph.png"), TStr("Graph"), GraphNodeLabels);

  // directed network with an integer value in each node, and a float value in each edge
  // (smart pointer to TIntFltNEDNet object)
  PIntFltNEDNet PNetwork = TIntFltNEDNet::New(); TIntStrH NetworkNodeLabels;

  // TODO: modify code to only connect nodes i -> j when value stored at node i is smaller than value in node j
  for (int i=0; i<10; i++) {
	  PNetwork->AddNode(i, TInt::Rnd.GetUniDevInt(100)); // Store a random value at node i
	  NetworkNodeLabels.AddDat(i) = TStr::Fmt("%d (%d)", i, PNetwork->GetNDat(i).Val); }
  for (int i=0; i<10; i++) {
	  for (int j=i+1; j<10; j++) { PNetwork->AddEdge(i, j, TFlt::Rnd.GetUniDev()); } // DAG
  }
  // TODO: END

  printf("\nDirected network (%d nodes, %d edges)\n", PNetwork->GetNodes(), PNetwork->GetEdges());

  // print nodes
  for (TIntFltNEDNet::TNodeI NI = PNetwork->BegNI(); NI < PNetwork->EndNI(); NI++) {
  	  printf("Node %d (val:%d)\n", NI.GetId(), NI.GetDat().Val);
  }

  // print edges
  for (TIntFltNEDNet::TEdgeI EI = PNetwork->BegEI(); EI < PNetwork->EndEI(); EI++) {
  	  printf("Edge %d -> %d (val:%f)\n", EI.GetSrcNId(), EI.GetDstNId(), EI.GetDat().Val);
  }

  // generate graphviz
  TSnap::DrawGViz(PNetwork, gvlDot, TStr("network.png"), TStr("Network"), NetworkNodeLabels);

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
