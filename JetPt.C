#include "BSHelper.cxx"
#include "DiJetCommon.h"
#include "Filipad2.h"
#include "TUnfold.h"

const float inf = 1e20;
using namespace std;

void JetPt(Bool_t ispp = 1 )
{

	Int_t colors[] = {kRed, kOrange, kGreen, kGreen + 3, kBlue, kBlue + 2, kViolet, kCyan + 2, kAzure, kMagenta, kMagenta + 3, kYellow + 3, kYellow + 2, kSpring - 7};
	double rangeMin = 10, rangeMax = 500;
	double ratioMin = 0.1, ratioMax = 1.2;
	double plotMin = 1e-7, plotMax = 1;


	TString fmc;
	TString fdata;
	TString mcname;
	TString dataname;

	if (ispp) {
		fmc = "DijetLHC17pMCAOD";
		fdata = "DijetLHC17pAOD";
		mcname = fmc;
		dataname = fdata;
	}
	else {
		fmc = "DijetLHC15oAODEmb";
		fdata = "DijetLHC15oAODEmb";
		mcname = fmc;
		dataname = "DijetLHC15oAOD";

	}

	const int npthardbin = 20; // LHC16j5 has 20 pt hard bins

	//MC pt-hard-bin merging
	THnSparse *HDiJetPtPair_true = nullptr;
	THnSparse *HDiJetPtPair_mcrec = nullptr;
	THnSparse *HRes = nullptr;
	THnSparse *HMiss = nullptr;
	THnSparse *HFake = nullptr;

	auto mclist = LoadDiJetResultList(fmc.Data(), mcname.Data());

	auto HT = BSTHnSparseHelper::Load("hJetPtTruth", mclist);
	auto HM = BSTHnSparseHelper::Load("hJetPt", mclist);
	auto HR = BSTHnSparseHelper::Load("hJetPtRes", mclist);
	auto HI = BSTHnSparseHelper::Load("hJetPtMiss", mclist);
	auto HF = BSTHnSparseHelper::Load("hJetPtFake", mclist);
	auto HN = BSTHnSparseHelper::Load("hPtHardBin", mclist);

	HT.PrintAxis("all");

	Double1D binCent = {0, 100}; // second bin
	Int_t tpt = -1;				 //third bin  // fourth bin dijet mass
	Int_t ptpair = -1;			 //fifth bin // sisxth pthardbin
	HT.SetBin("Cent", binCent);
	HM.SetBin("Cent", binCent);
	HR.SetBin("Cent", binCent);
	HI.SetBin("Cent", binCent);
	HF.SetBin("Cent", binCent);
	HN.SetBin("Cent", binCent);

	Int_t centbin = 1;
	if (ispp) centbin = 1;
	if (binCent[0]==0 && binCent[1] ==100) centbin = -1;

	auto hnpthardbin = HN.GetTH1("hnpt",1,{centbin,-1}); // # of events of each pt hard bin

    //---------------------------------------------------------------------------------------------
	


	TH1D *ht = nullptr;
	TH1D *hm = nullptr;
	TH2D *hr = nullptr;
	TH1D *hi = nullptr;
	TH1D *hf = nullptr;
	for (auto i = 0; i < 20; i++)
	{
		auto h = HT.GetTH1(Form("ht%d", i), 1, {centbin, -1, i + 1});
		h->Scale(1. / hnpthardbin->GetBinContent(i + 1), "");
		if (i == 0)
			ht = (TH1D *)h->Clone();
		else
			ht->Add(h);

		h = HM.GetTH1(Form("hm%d", i), 1, {centbin, -1, i + 1});
		h->Scale(1. / hnpthardbin->GetBinContent(i + 1), "");
		if (i == 0)
			hm = (TH1D *)h->Clone();
		else
			hm->Add(h);

		auto h2 = HR.GetTH2(Form("hr%d", i), 1, 2, {centbin, -1, -1, i + 1});
		h2->Scale(1. / hnpthardbin->GetBinContent(i + 1), "");
		if (i == 0)
			hr = (TH2D *)h2->Clone();
		else
			hr->Add(h2);

		h = HI.GetTH1(Form("hi%d", i), 1, {centbin, -1, i + 1});
		h->Scale(1. / hnpthardbin->GetBinContent(i + 1), "");
		if (i == 0)
			hi = (TH1D *)h->Clone();
		else
			hi->Add(h);

		h = HF.GetTH1(Form("hf%d", i), 1, {centbin, -1, i + 1});
		h->Scale(1. / hnpthardbin->GetBinContent(i + 1), "");
		if (i == 0)
			hf = (TH1D *)h->Clone();
		else
			hf->Add(h);
	}

	//hr -> Draw("colz");

	RooUnfoldResponse response(hm, ht);
	for (auto i = 0; i <= hr->GetNbinsX() +1 ; i++)
	{ //invM
		for (auto j = 0; j <= hr->GetNbinsY() + 1; j++)
		{ //ptpair
			Double_t bincenx = hr->GetXaxis()->GetBinCenter(i);
			Double_t binceny = hr->GetYaxis()->GetBinCenter(j);
			//if (bincenx< 20 || binceny < 20) continue;
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
		if (ispp) response.Fake(bincenx, bincont);
	}

	auto datalist = LoadDiJetResultList (fdata.Data(), dataname.Data());
	auto HD = BSTHnSparseHelper::Load( "hJetPt", datalist );
	HD.SetBin("Cent",binCent);
	auto hd = HD.GetTH1(Form("hdata"),1,{centbin,-1,-1});
	TH1D* NOE = (TH1D*) datalist -> FindObject("hEventNumbers");
	Double_t noe = NOE->GetBinContent(3);
	if (ispp) hd -> Scale(51.2/noe);
	else if (binCent[0] == 0 && binCent[1] == 10)
		hd -> Scale (1./(noe*0.1)/(26.32/2+20.56/2));
	else if (binCent[0] == 50 && binCent[1] == 80)
		hd -> Scale (1./(noe*0.3)/((0.4174*0.2+1.564*0.05+1.07*0.05)/0.3));
	else hd -> Scale(1./noe);

	auto hdd = (TH1D*) hd -> Clone();

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
	hc->GetXaxis()->SetRangeUser(20, 200);
	hc->SetMaximum(2e4);
	hc->SetMinimum(1e-4);
	hset(*hc, "#it{p}_{T}^{jet}", "d^{2}#it{#sigma}/d#it{p}_{T}^{jet}d#it{#eta} (mb #it{c}/GeV)", 0.9, 1.1, 0.07, 0.06, 0.01, 0.001, 0.04, 0.05, 510, 510);
	//hc -> Scale(1. ,"width");
	hc->Scale(1., "width");
	hc->Draw("PZsame");
	hd->Scale(1., "width");

	hm->Scale(0.1, "width");
	hm->SetLineColor(2);
	hm->SetMarkerColor(2);

	ht->Scale(0.1, "width");
	ht->SetLineColor(1);
	ht->SetMarkerColor(1);


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
		h->Draw("PZsame");
		lg->AddEntry(h, S.at(n).Data(), "lp");
		n++;
	}
	lg->Draw();

	lg = new TLegend(0.578947, 0.634783, 0.95933, 0.933913);
	lg->SetTextSize(0.0521739);
	lg->SetBorderSize(0);
	lg->AddEntry((TObject *)0x0, "Charged dijet", "");
	lg->AddEntry((TObject *)0x0, "R=0.4, Anti-kt", "");
	//lg->AddEntry((TObject*)0x0,"#it{p}_{T, jet} > 20 GeV/#it{c}","");
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
	hset(*r, "#it{p}_{T}^{jet} (GeV/#it{c})", "Ratio (DATA/MC)", 1.0, 0.8, 0.08, 0.09, 0.01, 0.001, 0.07, 0.06, 10, 510);

	r->SetMinimum(0);
	r->SetMaximum(2);
	r->Draw();

	fpadIncl = new Filipad2(1, 2, 0.5, 100, 50, 0.7, 1, 1);
	fpadIncl->Draw();
	p = fpadIncl->GetPad(1); //upper pad
	p->SetTickx();
	p->SetLogy(1);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	auto hF = (TH1D *)hc->Clone();

	hF->Draw("PZ same");

	TFile::Open("../Results/ALICE_pp_Pub_5TeVLTCut.root");
	auto g = (TGraphAsymmErrors *)gROOT->FindObject("Graph1D_y3");
	g->SetLineColor(2);
	g->SetMarkerColor(2);
	g->Draw("P1same");

	TFile::Open("../Results/ATLAS_pp_full.root");
	auto atlas = (TGraphAsymmErrors *)gROOT->FindObject("Graph1D_y1");
	for (int i = 0; i < atlas->GetN(); i++)
	{
		atlas->GetY()[i] *= 1e-6;
		atlas->GetEYhigh()[i] *= 1e-6;
		atlas->GetEYlow()[i] *= 1e-6;
	}

 
	//atlas->Draw("P1same");


	Double_t x[]={22,26,30,35,41,47,54,62,71,81,93};
	Double_t y[]={1.69442*1e-3,0.810588*1e-3,0.409895*1e-3,0.206366*1e-3,0.0932333*1e-3,0.0526556*1e-3,0.0274628*1e-3,0.0127448*1e-3,0.0057247*1e-3,0.0029927*1e-3,0.0015268*1e-3};
	Double_t ex[]={2,2,2,2.5,3,3,3.5,4,4.5,5,6};
	Double_t ey[]={0.0345867*1e-3,0.0182654*1e-3,0.0108276*1e-3,0.0055*1e-3,0.0033823*1e-3,0.0024233*1e-3,0.0015585*1e-3,0.0010182*1e-3,0.0005915*1e-3,0.0004336*1e-3,0.0002891*1e-3};
	TGraphErrors *gr = new TGraphErrors(11,x,y,ex,ey);
	//gr->Draw("p1same");

	lg = new TLegend(0.186603, 0.0330435, 0.490431, 0.252174);
	lg->SetTextSize(0.0521739);
	lg->SetBorderSize(0);
	lg->AddEntry(hF, "Corrected result at 5 TeV", "lp");
	lg->AddEntry(g, "ALICE published result", "lp");
	lg->Draw();

	p = fpadIncl->GetPad(2); //upper pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	auto rr = CalcRatio(hF,g);
	TH1D* r2 = (TH1D*) r -> Clone();
  hset(*r2, "#it{p}_{T}^{jet} (GeV/#it{c})", "Ratio (Published/Corrected)"
      ,1.0,0.8, 0.08,0.09, 0.01,0.001, 0.07,0.06, 10,510);
	r2->SetMaximum(2);
	r2->SetMinimum(0);
	r2->Reset();
	r2->Draw();
	rr -> Draw("p1same");

}