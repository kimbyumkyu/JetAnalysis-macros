#include "BSHelper.cxx"
#include "DiJetCommon.h"
#include "Filipad2.h"
#include "TUnfold.h"

const float inf = 1e20;
using namespace std;

void DijetMass(Bool_t ispp = 1)
{

	Int_t colors[] = {kRed, kOrange, kGreen, kGreen + 3, kBlue, kBlue + 2, kViolet, kCyan + 2, kAzure, kMagenta, kMagenta + 3, kYellow + 3, kYellow + 2, kSpring - 7};
	double rangeMin = 10, rangeMax = 500;
	double ratioMin = 0.1, ratioMax = 1.2;
	double plotMin = 1e-7, plotMax = 1;

	TString fmc;
	TString fdata;
	TString mcname;
	TString dataname;

	if (ispp)
	{
		fmc = "DijetLHC17pMCAOD";
		fdata = fmc;
		mcname = fmc;
		//fdata = "DijetLHC17pAOD";
		//dataname = fdata;
		dataname = mcname;
	}
	else
	{
		fmc = "DijetLHC15oAODEmb";
		fdata = "DijetLHC15oAODEmb";
		mcname = fmc;
		dataname = "DijetLHC15oAOD";
	}

	const int npthardbin = 20; // LHC16j5 has 20 pt hard bins
	THnSparse *HDiJetPtPair_true = nullptr;
	THnSparse *HDiJetPtPair_mcrec = nullptr;
	THnSparse *HRes = nullptr;
	THnSparse *HMiss = nullptr;
	THnSparse *HFake = nullptr;

	//auto datalist = LoadDiJetResultList ("DijetLHC15oAODEmb", "DijetLHC15oAOD");
	//auto mclist = LoadDiJetResultList ("DijetLHC15oAODEmb", "DijetLHC15oAODEmb");
	auto mclist = LoadDiJetResultList(fmc.Data(), mcname.Data());
	//auto hnpthardbin = (TH1D *)mclist->FindObject("hPtHardBin"); // # of events of each pt hard bin

	auto HT = BSTHnSparseHelper::Load("hDiJetInvMPtPairTruth", mclist);
	auto HM = BSTHnSparseHelper::Load("hDiJetInvMPtPair", mclist);
	auto HR = BSTHnSparseHelper::Load("hDiJetInvMPtPairRes", mclist);
	auto HI = BSTHnSparseHelper::Load("hDiJetInvMPtPairMiss", mclist);
	auto HF = BSTHnSparseHelper::Load("hDiJetInvMPtPairFake", mclist);
	auto HN = BSTHnSparseHelper::Load("hPtHardBin", mclist);

	HT.PrintAxis("all");

	Int_t idijet = 10;			 //first bin
	Double1D binCent = {0, 10}; // second bin
	Int_t tpt = -1;				 //third bin  // fourth bin dijet mass
	Int_t ptpair = -1;			 //fifth bin // sisxth pthardbin
	HT.SetBin("Cent", binCent);
	HM.SetBin("Cent", binCent);
	HR.SetBin("Cent", binCent);
	HI.SetBin("Cent", binCent);
	HF.SetBin("Cent", binCent);
	HN.SetBin("Cent", binCent);


	Int_t centbin = 1;
	if (binCent[0]==0 && binCent[1] ==100) centbin = -1;
	if (ispp) centbin = -1;

	auto hnpthardbin = HN.GetTH1("hnpt",1,{centbin,-1}); // # of events of each pt hard bin
	TH1D *ht = nullptr;
	TH1D *hm = nullptr;
	TH2D *hr = nullptr;
	TH1D *hi = nullptr;
	TH1D *hf = nullptr;
	for (auto i = 0; i < 20; i++)
	{
		auto h = HT.GetTH1(Form("ht%d", i), 2, {idijet, centbin, -1, -1, i + 1});
		h->Scale(1. / hnpthardbin->GetBinContent(i + 1), "");
		if (i == 0)
			ht = (TH1D *)h->Clone();
		else
			ht->Add(h);

		h = HM.GetTH1(Form("hm%d", i), 2, {idijet, centbin, -1, -1, i + 1});
		h->Scale(1. / hnpthardbin->GetBinContent(i + 1), "");
		if (i == 0)
			hm = (TH1D *)h->Clone();
		else
			hm->Add(h);

		auto h2 = HR.GetTH2(Form("hr%d", i), 2, 3, {idijet, centbin, -1, -1, -1, -1, i + 1});
		h2->Scale(1. / hnpthardbin->GetBinContent(i + 1), "");
		if (i == 0)
			hr = (TH2D *)h2->Clone();
		else
			hr->Add(h2);

		h = HI.GetTH1(Form("hi%d", i), 2, {idijet, centbin, -1, -1, i + 1});
		h->Scale(1. / hnpthardbin->GetBinContent(i + 1), "");
		if (i == 0)
			hi = (TH1D *)h->Clone();
		else
			hi->Add(h);

		h = HF.GetTH1(Form("hf%d", i), 2, {idijet, centbin, -1, -1, i + 1});
		h->Scale(1. / hnpthardbin->GetBinContent(i + 1), "");
		if (i == 0)
			hf = (TH1D *)h->Clone();
		else
			hf->Add(h);
	}

	//hr -> Draw("colz");

	RooUnfoldResponse response(hm, ht);
	for (auto i = 0; i <= hr->GetNbinsX() + 1; i++)
	{ //invM
		for (auto j = 0; j <= hr->GetNbinsY() + 1; j++)
		{ //ptpair
			Double_t bincenx = hr->GetXaxis()->GetBinCenter(i);
			Double_t binceny = hr->GetYaxis()->GetBinCenter(j);
			Double_t bincont = hr->GetBinContent(i, j);
			response.Fill(bincenx, binceny, bincont);
		}
	}
	for (auto i = 0; i <= hi->GetNbinsX() + 1; i++)
	{ //invM
		Double_t bincenx = hi->GetXaxis()->GetBinCenter(i);
		Double_t bincont = hi->GetBinContent(i);
		response.Miss(bincenx, bincont);
	}
	for (auto i = 0; i <= hf->GetNbinsX() + 1; i++)
	{ //invM
		Double_t bincenx = hf->GetXaxis()->GetBinCenter(i);
		Double_t bincont = hf->GetBinContent(i);
		response.Fake(bincenx, bincont);
	}

	auto datalist = LoadDiJetResultList(fdata.Data(), dataname.Data());
	auto HD = BSTHnSparseHelper::Load("hDiJetInvMPtPair", datalist);
	HD.SetBin("Cent", binCent);
	auto hd = HD.GetTH1(Form("hdata"), 2, {idijet, centbin, -1, -1, -1});
	TH1D *NOE = (TH1D *)datalist->FindObject("hEventNumbers");
	Double_t noe = NOE->GetBinContent(3);

	//https://alice-notes.web.cern.ch/system/files/notes/public/711/2018-06-18-ALICE_public_note.pdf
	//https://alice-notes.web.cern.ch/system/files/notes/analysis/818/2018-Apr-27-analysis_note-AN_ChargedJetPbPb2018.pdf
	if (ispp) hd -> Scale(51.2/noe);
	else if (binCent[0] == 0 && binCent[1] == 10)
		hd -> Scale (1./(noe*0.1)/(26.32/2+20.56/2));
	else if (binCent[0] == 20 && binCent[1] == 50)
		hd -> Scale (1./(noe*0.3)/(4.98/3+2.659/3+8.675/3));
	else if (binCent[0] == 50 && binCent[1] == 80)
		hd -> Scale (1./(noe*0.3)/((0.4174*0.2+1.564*0.05+1.07*0.05)/0.3));
	else hd -> Scale(1./noe);

	//hd -> Draw();
	RooUnfoldBayes unfold(&response, hd, 4);

	Filipad2 *fpadIncl = new Filipad2(1, 2, 0.5, 100, 50, 0.7, 1, 0);
	fpadIncl->Draw();
	TPad *p = fpadIncl->GetPad(1); //upper pad
	p->SetTickx();
	p->SetLogy(1);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TH1D *hc = (TH1D *)unfold.Hreco();
	hc->SetMarkerColor(4);
	hc->SetLineColor(4);
	hc->Scale(1., "width");
	hc->GetXaxis()->SetRangeUser(40, 1000);
	hc->SetMarkerStyle(24);
	hc->SetMaximum(5e-3);
	hc->SetMinimum(1e-11);
	hset(*hc, "#it{m}_{jj}", "d^{2}#it{#sigma}/d#it{m}_{jj}d#it{#eta} (mb #it{c^{2}}/GeV)", 0.9, 1.1, 0.07, 0.06, 0.01, 0.001, 0.04, 0.05, 510, 510);

	hc->Draw("PZsame");
	hm->Scale(0.1, "width");
	hm->Draw("same");
	ht->Scale(0.1, "width");
	ht->Draw("same");
	hd->Scale(1., "width");
	hd->Draw("same");

	TLegend *lg = new TLegend(0.186603, 0.0330435, 0.490431, 0.252174);
	lg->SetTextSize(0.0521739);
	lg->SetBorderSize(0);
	vector<TH1D *> H = {hc, hd, ht, hm};
	vector<TString> S = {"Corrected data", "Raw data", "MC truth #times 0.1", "MC data #times 0.1"};
	vector<Int_t> C = {1, 1, 2, 2};
	vector<Int_t> M = {20, 24, 20, 24};
	int n = 0;
	for (auto h : H)
	{
		h->SetMarkerStyle(M.at(n));
		h->SetMarkerColor(C.at(n));
		h->SetLineColor(C.at(n));
		lg->AddEntry(h, S.at(n).Data(), "lp");
		n++;
	}
	lg->Draw();

	lg = new TLegend(0.578947, 0.634783, 0.95933, 0.933913);
	lg->SetTextSize(0.0521739);
	lg->SetBorderSize(0);
	lg->AddEntry((TObject *)0x0, "Charged jet", "");
	lg->AddEntry((TObject *)0x0, "R=0.4, Anti-kt", "");
	lg->AddEntry((TObject *)0x0, "#it{p}_{T, leading}^{track} >5 GeV/#it{c}", "");
	lg->Draw();

	p = fpadIncl->GetPad(2); //upper pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	auto r = (TH1D *)hc->Clone();
	auto htt = (TH1D *)ht->Clone();
	htt->Scale(10);
	r->Divide(htt);
	hset(*r, "#it{m}_{jj} (GeV/#it{c^{2}})", "Ratio (DATA/MC)", 1.0, 0.8, 0.08, 0.09, 0.01, 0.001, 0.07, 0.06, 10, 510);

	r->SetMinimum(0);
	r->SetMaximum(2);
	r->Draw();
}

/*

		ht -> Scale (1./hjetpttruth->Integral());
		ht -> GetXaxis()->SetRangeUser(0,400);
		ht -> SetMaximum(1000);
		ht -> Draw("PZ same");

	

		auto HD = BSTHnSparseHelper::Load( "hDiJetInvMPtPair", datalist );
		HD.SetBin("Cent",binCent);


    TLegend *legr = new TLegend(0.7,0.7,0.9,0.9,"pp #sqrt{#it{s}} = 5.02 TeV");
    legr->SetTextSize(0.028); legr->SetBorderSize(0);

		vector<TH1D*> H;
		for (int i=1; i<= binCent.size()-1; i++){
			auto hjetptdata = Hjetptdata.GetTH1(Form("hdata"),3,{idijet,i,tpt,-1,ptpair,-1});
			hjetptdata -> Scale(1,"width");
			hjetptdata->SetLineColor(colors[i-1]);
			hjetptdata->SetMarkerColor(colors[i-1]);
			hjetptdata -> Draw("PZ same");
			auto r = (TH1D*) hjetptdata -> Clone();
			r-> Scale(1./r->Integral());
			r -> Divide (hjetpttruth);
			H.push_back(r);
			legr -> AddEntry(hjetptdata,Form("Cent : %d - %d",int(binCent[i-1]),int(binCent[i])),"lp");
		}
		legr -> Draw();

    p = fpadIncl->GetPad(2); //upper pad
    p->SetTickx(); p->SetGridy(1); p->SetGridx(1); p->SetLogy(0); p->cd();
    gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
		p -> DrawFrame(0,0,400,3);
		for (auto h : H){
			h -> Draw("pzsame");

		}
*/
