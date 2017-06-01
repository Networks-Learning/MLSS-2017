#include "stdafx.h"
#include "mlss.h"
#include "kronecker.h"

void TMLSS::LoadNetworkTxt(TSIn& SIn, const double& Prob) {
	bool verbose = false;
	TStr Line;

	Network->Clr(); // clear network (if any)

	// add nodes
	SIn.GetNextLn(Line);
	while (!SIn.Eof() && Line != "") {
		TStrV NIdV; Line.SplitOnAllCh(',', NIdV);
		Network->AddNode(NIdV[0].GetInt(), NIdV[1]);
		if (!IsNodeNm(NIdV[0].GetInt())) {
			AddNodeNm(NIdV[0].GetInt(), TNodeInfo(NIdV[1], 0));
		}

		if (verbose) { printf("Node %d (%s)\n", NIdV[0].GetInt(), NIdV[1].CStr()); }
		SIn.GetNextLn(Line);
	}

	if (verbose) { printf("\n"); }

	// add edges
	while (!SIn.Eof()) {
		SIn.GetNextLn(Line);
		TStrV FieldsV; Line.SplitOnAllCh(',', FieldsV);

		if (FieldsV.Len()<2) { continue; }
		if (FieldsV.Len()==2) {
			Network->AddEdge(FieldsV[0].GetInt(), FieldsV[1].GetInt(), Prob);
			if (verbose) { printf("Edge %d -> %d: %f\n", FieldsV[0].GetInt(), FieldsV[1].GetInt(), Prob); }
		} else {
			Network->AddEdge(FieldsV[0].GetInt(), FieldsV[1].GetInt(), FieldsV[2].GetFlt());
			if (verbose) { printf("Edge %d -> %d: %f\n", FieldsV[0].GetInt(), FieldsV[1].GetInt(), FieldsV[2].GetFlt()); }
		}
	}

	printf("\nNetwork loaded -> nodes:%d edges:%d\n", Network->GetNodes(), Network->GetEdges());
}

void TMLSS::GenerateNetwork(const int& TNetwork, const int& NNodes, const int& NEdges, const TStr& NetworkParams, const bool& verbose) {
	PNGraph Graph;
	PUNGraph UGraph;

	// set random seed
	TInt::Rnd.Randomize();

	// create random network structure
	switch (TNetwork) {
	// 2-dimension kronecker network
	case 0:
		{
			printf("Kronecker graph\n");
			TStr MtxNm;
			TKronMtx SeedMtx = TKronMtx::GetMtx(NetworkParams.CStr()); // e.g. 0.5 0.5; 0.5 0.5

			printf("\n*** Seed matrix:\n");
			SeedMtx.Dump();

			// generate directed kronecker graph
			Graph = TKronMtx::GenFastKronecker(SeedMtx, (int)TMath::Log2(NNodes), NEdges, true, 0);

			break;
		}
	// forest fire network
	case 1:
		{
			printf("Forest Fire graph\n");
			TStrV NetworkParamsV; NetworkParams.SplitOnAllCh(';', NetworkParamsV);

			// generate directed forest fire graph (e.g., ForwBurnProb (0.2), BackBurnProb (0.17))
			Graph = TSnap::GenForestFire(NNodes, NetworkParamsV[0].GetFlt(), NetworkParamsV[1].GetFlt());

			break;
		}
	// TODO: other type of random graph
	case 2:
		{
			printf("Other random graph\n");
			break;
		}
	// TODO END
	}

	if (verbose) { printf("\n"); }

	// fill network structure with graph (do not allow self loops!)
	for (TNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
		Network->AddNode(NI.GetId(), TStr::Fmt("%d", NI.GetId()));
		AddNodeNm(NI.GetId(), TNodeInfo(TStr::Fmt("%d", NI.GetId()), 0));
		if (verbose) { printf("Node %d\n", NI.GetId()); }
	}

	if (verbose) { printf("\n"); }

	int selfloops = 0;
	for (TNGraph::TEdgeI EI = Graph->BegEI(); EI < Graph->EndEI(); EI++) {
		if (EI.GetSrcNId()==EI.GetDstNId()) { selfloops++; continue; }
		Network->AddEdge(EI.GetSrcNId(),EI.GetDstNId());
		if (verbose) { printf("Edge %d -> %d\n", EI.GetSrcNId(), EI.GetDstNId()); }
	}

	if (verbose) { printf("\n"); }

	if (verbose) {
		printf("Network structure has been generated successfully!\n");
		printf("%d self-loops has been removed\n", selfloops);
	}
}

void TMLSS::GenerateAlphas(const double& MinAlpha, const double& MaxAlpha) {
	TInt::Rnd.Randomize();

	IAssert(MinAlpha>=0.0 && MaxAlpha>=MinAlpha);

	// generate alphas uniformly at random
	for (TStrFltNEDNet::TEdgeI EI = Network->BegEI(); EI<Network->EndEI(); EI++) {
		EI.GetDat() = MinAlpha + (MaxAlpha-MinAlpha)*TFlt::Rnd.GetUniDev();
	}
}

void TMLSS::GenerateBetas(const double& MinBeta, const double& MaxBeta) {
	TInt::Rnd.Randomize();

	IAssert(MinBeta>=0 && MaxBeta>=MinBeta && MaxBeta<=1);

	// generate betas uniformly at random
	for (TStrFltNEDNet::TEdgeI EI = Network->BegEI(); EI<Network->EndEI(); EI++) {
		EI.GetDat() = MinBeta + (MaxBeta-MinBeta)*TFlt::Rnd.GetUniDev();
	}
}

void TMLSS::ConvertAlphasToBetas(const double& T) {
	TIntPrV EdgesToDelete;

	for (TStrFltNEDNet::TEdgeI EI = Network->BegEI(); EI<Network->EndEI(); EI++) {
		switch (Model) {
		case EXP:
			EI.GetDat() = 1-exp(-EI.GetDat()*T);
			break;
		case POW:
			EI.GetDat() = 1-pow(T, -EI.GetDat());
			break;
		case RAY:
			EI.GetDat() = 1-exp(0.5*EI.GetDat()*pow(T,2.0));
			break;
		default:
			break;
		}

		if (EI.GetDat()<0.01) { EdgesToDelete.Add(TIntPr(EI.GetSrcNId(), EI.GetDstNId())); }
	}

	for (int i=0; i<EdgesToDelete.Len(); i++) { Network->DelEdge(EdgesToDelete[i].Val1, EdgesToDelete[i].Val2); }
}

// it generates a cascade starting from Sources using independent cascade model
void TMLSS::GenCascadeIC(TCascade& C, const bool& verbose) {
	TIntH InfectedNIdH, InfectedBy;
	int GlobalTime = 0; // cascade starts at 0; discrete time
	double beta = 0.0;
	int StartNId = -1;

	// if network is empty, stop
	if (Network->GetNodes() == 0) { return; }

	// reset cascade, and hash tables
	C.Clr();
	InfectedNIdH.Clr();
	InfectedBy.Clr();

	// choose first node (if more than one source, choose at random)
	if (SourcesV.Len()==0) { StartNId = Network->GetRndNId(); }
	else { StartNId = SourcesV[TInt::Rnd.GetUniDevInt(0, SourcesV.Len()-1)].Val; }

	// Infection time of all sources is 0
	for (int i=0; i<SourcesV.Len(); i++) { InfectedNIdH.AddDat(SourcesV[i]) = GlobalTime; }
	if (SourcesV.Len()==0) { InfectedNIdH.AddDat(StartNId) = GlobalTime; }

	while (true) {
		// sort by time & get the oldest node that did not run infection
		InfectedNIdH.SortByDat(true);
		const int& NId = InfectedNIdH.BegI().GetKey();
		GlobalTime = InfectedNIdH.BegI().GetDat();

		// all the nodes has run infection
		if ( GlobalTime >= TInt::Mx) { break; }

		// add node to graph if it does not exist
		if (!C.Network->IsNode(NId)) { C.Network->AddNode(NId); }

		// add infected node with first infection & first parent
		C.Add(NId, (InfectedBy.IsKey(NId)? InfectedBy.GetDat(NId).Val : -1), (double)GlobalTime);

		if (verbose) { printf("GlobalTime:%d, infected node:%d, first parent:%d\n", GlobalTime, NId, (InfectedBy.IsKey(NId)? InfectedBy.GetDat(NId).Val : -1)); }

		// run infection from the current oldest node
		TStrFltNEDNet::TNodeI NI = Network->GetNI(NId);
		for (int e = 0; e < NI.GetOutDeg(); e++) {
			const int DstNId = NI.GetOutNId(e);

			// do not infect the parent
			if (InfectedBy.IsKey(NId) && InfectedBy.GetDat(NId).Val == DstNId)
				continue;

			beta = Network->GetEDat(NId, DstNId);

			if (verbose) { printf("GlobalTime:%d, nodes:%d->%d, beta:%f\n", GlobalTime, NId, DstNId, beta); }

			if (TFlt::Rnd.GetUniDev() >= beta) { continue; } // if infection does not succeed, skip

			int t1 = GlobalTime + 1; // infections are always in discrete rounds

			if (!InfectedNIdH.IsKey(DstNId)) {
				// update cascade graph
				if (!C.Network->IsNode(DstNId)) { C.Network->AddNode(DstNId); }
				if (!C.Network->IsEdge(NId, DstNId)) { C.Network->AddEdge(NId, DstNId, (double)t1); }

				// update first parent
				InfectedNIdH.AddDat(DstNId) = t1;
				InfectedBy.AddDat(DstNId) = NId;
			} else {
				int t2 = InfectedNIdH.GetDat(DstNId);

				// update cascade graph
				if ( t2 != TInt::Mx ) {
					if (!C.Network->IsNode(DstNId)) { C.Network->AddNode(DstNId); }
					if (!C.Network->IsEdge(NId, DstNId)) { C.Network->AddEdge(NId, DstNId, (double)t1); }
				}
			}
		}

		// we cannot delete key (otherwise, we cannot sort), so we assign a big time (TInt::Mx)
		InfectedNIdH.GetDat(NId) = TInt::Mx;
	}

	if (verbose) { printf("%d nodes infected!\n", C.Len()); }

	C.Sort();
}

// it generates a cascade starting from Sources using continuous time independent cascade model
void TMLSS::GenCascadeCTIC(TCascade& C, const bool& verbose) {
	TIntFltH InfectedNIdH; TIntH InfectedBy;
	double GlobalTime = 0.0; // cascade starts at 0.0; continuous time
	double alpha = 0.0;
	int StartNId = -1;

	// if network is empty, stop
	if (Network->GetNodes() == 0) { return; }

	// if window is equal or less than 0, stop
	if (Window.Val <= 0) { return; }

	// reset cascade, and hash tables
	C.Clr();
	InfectedNIdH.Clr();
	InfectedBy.Clr();

	// choose first node (if more than one source, choose at random)
	if (SourcesV.Len()==0) { StartNId = Network->GetRndNId(); }
	else { StartNId = SourcesV[TInt::Rnd.GetUniDevInt(0, SourcesV.Len()-1)].Val; }

	// Infection time of all sources is 0
	for (int i=0; i<SourcesV.Len(); i++) { InfectedNIdH.AddDat(SourcesV[i]) = GlobalTime; }
	if (SourcesV.Len()==0) { InfectedNIdH.AddDat(StartNId) = GlobalTime; }

	while (true) {
		// sort by time & get the oldest node that did not run infection
		InfectedNIdH.SortByDat(true);
		const int& NId = InfectedNIdH.BegI().GetKey();
		GlobalTime = InfectedNIdH.BegI().GetDat();

		// all the nodes has run infection
		if ( GlobalTime >= Window )
			break;

		// add node to graph if it does not exist
		if (!C.Network->IsNode(NId)) { C.Network->AddNode(NId); }

		// add infected node with first infection & first parent
		C.Add(NId, (InfectedBy.IsKey(NId)? InfectedBy.GetDat(NId).Val : -1), GlobalTime);

		if (verbose) { printf("GlobalTime:%f, infected node:%d, first parent:%d\n", GlobalTime, NId, (InfectedBy.IsKey(NId)? InfectedBy.GetDat(NId).Val : -1)); }

		// run infection from the current oldest node
		TStrFltNEDNet::TNodeI NI = Network->GetNI(NId);
		for (int e = 0; e < NI.GetOutDeg(); e++) {
			const int DstNId = NI.GetOutNId(e);

			// do not infect the parent
			if (InfectedBy.IsKey(NId) && InfectedBy.GetDat(NId).Val == DstNId)
				continue;

			alpha = Network->GetEDat(NId, DstNId);

			if (verbose) { printf("GlobalTime:%f, nodes:%d->%d, alpha:%f\n", GlobalTime, NId, DstNId, alpha); }

			// draw delay from pairwise distribution
			double sigmaT = TFlt::Mx;
			switch (Model) {
			case EXP:
				// exponential with alpha parameter
				sigmaT = TInt::Rnd.GetExpDev(alpha);
				break;
			case POW:
				// power-law with alpha parameter
				sigmaT = Delta*TInt::Rnd.GetPowerDev(1+alpha);
				while (sigmaT < Delta) { sigmaT = Delta*TInt::Rnd.GetPowerDev(1+alpha); }
				break;
			case RAY:
				// rayleigh with alpha parameter
				sigmaT = TInt::Rnd.GetRayleigh(1/sqrt(alpha));
				break;
			default:
				sigmaT = 1;
				break;
			}

			IAssert(sigmaT >= 0);

			double t1 = GlobalTime + sigmaT;

			if (InfectedNIdH.IsKey(DstNId)) {
				double t2 = InfectedNIdH.GetDat(DstNId);

				// update cascade graph
				if ( t2 != Window ) {
					if (!C.Network->IsNode(DstNId)) { C.Network->AddNode(DstNId); }
					if (!C.Network->IsEdge(NId, DstNId)) { C.Network->AddEdge(NId, DstNId, t1); }
				}

				// update first parent
				if ( t2 > t1 && t2 != Window ) {
					InfectedNIdH.GetDat(DstNId) = t1;
					InfectedBy.GetDat(DstNId) = NId;
				}
			} else {
				// update cascade graph
				if (!C.Network->IsNode(DstNId)) { C.Network->AddNode(DstNId); }
				if (!C.Network->IsEdge(NId, DstNId)) { C.Network->AddEdge(NId, DstNId, t1); }

				// update first parent
				InfectedNIdH.AddDat(DstNId) = t1;
				InfectedBy.AddDat(DstNId) = NId;
			}
		}

		// we cannot delete key (otherwise, we cannot sort), so we assign a big time (window cut-off)
		InfectedNIdH.GetDat(NId) = Window;
	}

	C.Sort();
}

void TMLSS::SetSources(const TStrV& SourcesStrV) {
	SourcesV.Clr();
	Assert(SourcesStrV.Len()<=Network->GetNodes());
	for (int i=0; i<SourcesStrV.Len(); i++) {
		SourcesV.Add(SourcesStrV[i].GetInt());
	}
}

void TMLSS::SetSources(const int& NSources) {
	SourcesV.Clr();
	Assert(NSources<=Network->GetNodes());
	while (SourcesV.Len()<NSources) {
		SourcesV.AddUnique(Network->GetRndNId());
	}
}

void TMLSS::SaveCascadesTxt(const TStr& OutFNm) {
	TFOut FOut(OutFNm);

	// write nodes to file
	for (THash<TInt, TNodeInfo>::TIter NI = NodeNmH.BegI(); NI < NodeNmH.EndI(); NI++) {
		FOut.PutStr(TStr::Fmt("%d,%s\r\n", NI.GetKey().Val, NI.GetDat().Name.CStr()));
	}

	FOut.PutStr("\r\n");

	// write cascades to file, nodes in a cascade separated by ; and first number is Cascade ID
	for (THash<TInt, TCascade>::TIter CI = CascH.BegI(); CI < CascH.EndI(); CI++) {
		TCascade &C = CI.GetDat();
		int j = 0;
		for (THash<TInt, THitInfo>::TIter NI = C.NIdHitH.BegI(); NI < C.NIdHitH.EndI(); NI++, j++) {
			if (!NodeNmH.IsKey(NI.GetDat().NId)) { continue; }
			if (j > 0) { FOut.PutStr(TStr::Fmt(";%d,%f", NI.GetDat().NId.Val, NI.GetDat().Tm.Val)); }
			else { FOut.PutStr(TStr::Fmt("%d;%d,%f", CI.GetKey().Val, NI.GetDat().NId.Val, NI.GetDat().Tm.Val)); }
		}

		if (j >= 1) { FOut.PutStr(TStr::Fmt("\r\n")); }
	}
}

void TMLSS::SaveCascadesGephi(const TStr& OutFNm, const int& NCascades, const bool& PrintSources, const double& TimeSep) {
	TFOut FOut(OutFNm);

	// print headers for GEXF
	FOut.PutStr("<gexf xmlns=\"http://www.gexf.net/1.2draft\"\nxmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\nxmlns:viz=\"http://www.gexf.net/1.1draft/viz\"\nxsi:schemaLocation=\"http://www.gexf.net/1.1draft\"\nversion=\"1.1\">\n\n");
	FOut.PutStr("<graph mode=\"dynamic\" defaultedgetype=\"directed\" timeformat=\"float\">\n\n");

	// print nodes section (sources/normal nodes)
	FOut.PutStr("\t<nodes>\n\n");
	for (THash<TInt, TNodeInfo>::TIter NI = NodeNmH.BegI(); NI < NodeNmH.EndI(); NI++) {
		TStr NodeStr(TStr::Fmt("\t\t<node id=\"%d\" label=\"%s\">\n\n", NI.GetKey().Val, NI.GetDat().Name.CStr()));

		// red nodes are sources, black nodes are normal nodes
		if (PrintSources && IsSource(NI.GetKey())) { NodeStr += "\t\t<viz:color r=\"255\" g=\"0\" b=\"0\" a=\"1.0\" />\n"; }
		else { NodeStr += "\t\t<viz:color r=\"0\" g=\"0\" b=\"0\" a=\"1.0\" />\n"; }

		NodeStr += TStr::Fmt("\t\t<attvalues>\n\t\t\t<attvalue for=\"volume\" value=\"%d\"/>\n\t\t</attvalues>\n", NI.GetDat().Vol.Val);

		FOut.PutStr(NodeStr);
		FOut.PutStr("\t\t</node>\n");
	}
	FOut.PutStr("\n\t</nodes>\n\n");

	// print edges section
	FOut.PutStr("\t<edges>\n\n");

	for (TStrFltNEDNet::TEdgeI EI = Network->BegEI(); EI < Network->EndEI(); EI++) {
		if (!NodeNmH.IsKey(EI.GetSrcNId()) || !NodeNmH.IsKey(EI.GetDstNId())) { continue; }

		TStr EdgeStr(TStr::Fmt("\t\t<edge source=\"%d\" target=\"%d\">\n\t\t<spells>\n", EI.GetSrcNId(), EI.GetDstNId()));

		// sweet through cascades to find when the current edge was active
		int k = 0;
		double c_time = 0;
		for (int i=0; i<TMath::Mn(NCascades, CascH.Len()); i++) {
			if (CascH[i].Network->IsEdge(EI.GetSrcNId(), EI.GetDstNId())) {
				EdgeStr += TStr::Fmt("\t\t\t<spell start=\"%f\" ", CascH[i].Network->GetEDat(EI.GetSrcNId(), EI.GetDstNId()).Val+c_time);
				EdgeStr += TStr::Fmt("end=\"%f\" />\n", CascH[i].LastTm()+c_time+TimeSep);
				k++;
			}

			if (CascH[i].Len()==1) { c_time += TimeSep*2.0; }
			else { c_time += (CascH[i].LastTm() - CascH[i].FirstTm() + TimeSep*2.0); }
		}

		EdgeStr += "\t\t</spells>\n";

		if (k > 0) {
			FOut.PutStr(EdgeStr);
			FOut.PutStr("\t\t</edge>\n");
		}
	}

	FOut.PutStr("\n\t</edges>\n\n");

	// finish gexf file
	FOut.PutStr("</graph>\n</gexf>\n");
}

void TMLSS::SaveCascadesMatlab(const TStr& OutFNm) {
	TFOut FOut(OutFNm);

	// we do not write nodes to file, only write cascadeid, node, timestamp
	int j = 0;
	for (THash<TInt, TCascade>::TIter CI = CascH.BegI(); CI < CascH.EndI(); CI++, j++) {
		TCascade &C = CI.GetDat();
		for (THash<TInt, THitInfo>::TIter NI = C.NIdHitH.BegI(); NI < C.NIdHitH.EndI(); NI++) {
			FOut.PutStr(TStr::Fmt("%d\t%d\t%f\n", j+1, NodeNmH.GetKeyId(NI.GetDat().NId.Val)+1, NI.GetDat().Tm.Val));
		}
	}
}

void TMLSS::SaveNetworkTxt(const TStr& OutFNm) {
	TFOut FOut(OutFNm);

	// write nodes to file
	for (THash<TInt, TNodeInfo>::TIter NI = NodeNmH.BegI(); NI < NodeNmH.EndI(); NI++) {
		FOut.PutStr(TStr::Fmt("%d,%s\r\n", NI.GetKey().Val, NI.GetDat().Name.CStr()));
	}

	FOut.PutStr("\r\n");

	// write edges to file (not allowing self loops in the network)
	for (TStrFltNEDNet::TEdgeI EI = Network->BegEI(); EI < Network->EndEI(); EI++) {
		if (!NodeNmH.IsKey(EI.GetSrcNId()) || !NodeNmH.IsKey(EI.GetDstNId())) { continue; }

		TStr Line;
		FOut.PutStr(TStr::Fmt("%d,%d,%f\r\n", EI.GetSrcNId(), EI.GetDstNId(), EI().Val));
	}
}

void TMLSS::SaveNetworkGephi(const TStr& OutFNm) {
	TFOut FOut(OutFNm);

	// print headers for GEXF
	FOut.PutStr("<gexf xmlns=\"http://www.gexf.net/1.2draft\"\nxmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\nxmlns:viz=\"http://www.gexf.net/1.1draft/viz\"\nxsi:schemaLocation=\"http://www.gexf.net/1.1draft\"\nversion=\"1.1\">\n\n");
	FOut.PutStr("<graph mode=\"static\" defaultedgetype=\"directed\">\n\n");
	FOut.PutStr("\t<nodes>\n\n");
	for (THash<TInt, TNodeInfo>::TIter NI = NodeNmH.BegI(); NI < NodeNmH.EndI(); NI++) {
		// if (Network->GetNI(NI.GetKey()).GetDeg() < 1) { continue; } // uncomment to skip isolated nodes

		TStr NodeStr(TStr::Fmt("\t\t<node id=\"%d\" label=\"%s\">\n\n", NI.GetKey().Val, NI.GetDat().Name.CStr()));

		FOut.PutStr(NodeStr);
		FOut.PutStr("\t\t</node>\n");
	}

	FOut.PutStr("\n\t</nodes>\n\n");

	FOut.PutStr("\t<edges>\n\n");

	for (TStrFltNEDNet::TEdgeI EI = Network->BegEI(); EI < Network->EndEI(); EI++) {
		if (!NodeNmH.IsKey(EI.GetSrcNId()) || !NodeNmH.IsKey(EI.GetDstNId())) { continue; }

		FOut.PutStr(TStr::Fmt("\t\t<edge source=\"%d\" target=\"%d\" weight=\"%f\">\n", EI.GetSrcNId(), EI.GetDstNId(), EI.GetDat().Val));
		FOut.PutStr("\t\t</edge>\n");
	}

	FOut.PutStr("\n\t</edges>\n\n");
	FOut.PutStr("</graph>\n</gexf>\n");
}

void TMLSS::SaveNetworkMatlab(const TStr& OutFNm) {
	TFOut FOut(OutFNm);

	// write edges to file (not allowing self loops in the network)
	for (TStrFltNEDNet::TEdgeI EI = Network->BegEI(); EI < Network->EndEI(); EI++) {
		if (!NodeNmH.IsKey(EI.GetSrcNId()) || !NodeNmH.IsKey(EI.GetDstNId())) { continue; }

		FOut.PutStr(TStr::Fmt("%d\t%d\t%f\r\n", EI.GetSrcNId()+1, EI.GetDstNId()+1, EI.GetDat().Val));
	}
}
