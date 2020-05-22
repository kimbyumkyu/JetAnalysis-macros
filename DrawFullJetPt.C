#include "BSHelper.cxx"
#include "JetCommon.h"
#include "Filipad2.h"
#include "TUnfold.h"

const float inf = 1e20;
using namespace std;
Int_t colors[] = {kRed, kOrange, kGreen, kGreen + 3, kBlue, kBlue + 2, kViolet, kCyan + 2, kAzure, kMagenta, kMagenta + 3, kYellow + 3, kYellow + 2, kSpring - 7};

TH1D *getHisto(TFile* fIn, TString sHistoName);
TGraphErrors *CalculateRatio( TGraphErrors* gr1, TGraphErrors* gr2, double xshift=0.0, bool invert=false);

TH1D* DrawFullJetPt(TString ndata = "JtLHC13cAOD", TString nmc = "JtLHC13cMCAOD", Bool_t scale = true)
{

	gStyle->SetErrorX(0.001);

	Double1D binCent = {0, 100}; // centbin merging
	const int npthardbin = 10; // LHC13b4_plus has 10 pt hard bins

	auto datalist = LoadJetResultList (ndata.Data(), ndata.Data());
	auto DPT = BSTHnSparseHelper::Load( "hFullJetPt", datalist ); // data jet pt
	TH1D* NOE = (TH1D*) datalist -> FindObject("hEventNumbers");
	Double_t noe = NOE->GetBinContent(3);

	// Basic idea : corrected jet pt = DPT* (TPT/MPT), 
	auto mclist = LoadJetResultList (nmc.Data(), nmc.Data());
	auto TPT = BSTHnSparseHelper::Load( "hTrueFullJetPt", mclist ); // MC true jet pt
	auto MPT = BSTHnSparseHelper::Load( "hFullJetPt", mclist ); // MC rec jet pt
	auto RES = BSTHnSparseHelper::Load( "hFullJetPtRes", mclist ); // Response matrix 
	auto FAKE = BSTHnSparseHelper::Load( "hFullJetPtFake", mclist ); // Fake jet pt in response matrix 
	auto MISS = BSTHnSparseHelper::Load( "hFullJetPtMiss", mclist ); // Missing jet pt in response matrix 
	auto HN = BSTHnSparseHelper::Load( "hPtHardBin", mclist ); // # of events for each pt hard bins

	// Set centrality range {0,100} 
	DPT.SetBin("Cent",binCent);
	TPT.SetBin("Cent",binCent);
	MPT.SetBin("Cent",binCent);
	RES.SetBin("Cent",binCent);
	FAKE.SetBin("Cent",binCent);
	MISS.SetBin("Cent",binCent);
	HN.SetBin("Cent",binCent);
	auto hn = HN.GetTH1("hn",1,{-1,-1}); // # of events of each pt hard bin by projection
	auto dpt = DPT.GetTH1("hdpt",1,{-1,-1, 1}); // # of events of each pt hard bin by projection

	//-----------------------------------------------------------
    //pT hard bin scaling for all MC histograms
	TH1D *tpt = nullptr;
	TH1D *mpt = nullptr;
	TH2D *res = nullptr;
	TH1D *fake = nullptr;
	TH1D *miss = nullptr;
	for (auto i = 0; i < npthardbin; i++)
	{
		auto h = TPT.GetTH1("htpt",1, {-1,-1,i+1});
		h->Scale(1. / hn->GetBinContent(i + 1), "");
		if (i == 0)
			tpt = (TH1D *) h->Clone();
		else
			tpt->Add(h);
		
		h = MPT.GetTH1("hmpt",1, {-1,-1,i+1});
		h->Scale(1. / hn->GetBinContent(i + 1), "");
		if (i == 0)
			mpt = (TH1D *)h->Clone();
		else
			mpt->Add(h);

		h = FAKE.GetTH1("hfake",1, {-1,-1,i+1});
		h->Scale(1. / hn->GetBinContent(i + 1), "");
		if (i == 0)
			fake = (TH1D *)h->Clone();
		else
			fake->Add(h);
			
		h = MISS.GetTH1("hmiss",1, {-1,-1,i+1});
		h->Scale(1. / hn->GetBinContent(i + 1), "");
		if (i == 0)
			miss = (TH1D *)h->Clone();
		else
			miss->Add(h);
		
		auto h2 = RES.GetTH2("hres",1,2, {-1,-1,-1,i+1});
		h2->Scale(1. / hn->GetBinContent(i + 1), "");
		if (i == 0)
			res = (TH2D *)h2->Clone();
		else
			res->Add(h2);

	}

	// Make a dummy response matrix
	RooUnfoldResponse response(mpt, tpt);
	// Filling the response matrix
	for (auto i = 0; i <= res->GetNbinsX() +1 ; i++)
	{ 
		for (auto j = 0; j <= res->GetNbinsY() + 1; j++)
		{ //ptpair
			Double_t bincenx = res->GetXaxis()->GetBinCenter(i);
			Double_t binceny = res->GetYaxis()->GetBinCenter(j);
			Double_t bincont = res->GetBinContent(i, j);
			response.Fill(bincenx, binceny, bincont);
		}
	}
	for (auto i = 0; i <= miss->GetNbinsX() + 1; i++)
	{ 
		Double_t bincenx = miss->GetXaxis()->GetBinCenter(i);
		Double_t bincont = miss->GetBinContent(i);
		response.Miss(bincenx, bincont);
	}
	for (auto i = 0; i <= fake->GetNbinsX() + 1; i++)
	{ 
		Double_t bincenx = fake->GetXaxis()->GetBinCenter(i);
		Double_t bincont = fake->GetBinContent(i);
		response.Fake(bincenx, bincont);
	} 

	RooUnfoldBayes unfold(&response, dpt, 4);
	auto hfinal = (TH1D*) unfold.Hreco();

	// Let's draw the results
	Filipad2 *fpadIncl = new Filipad2(1, 2, 0.5, 100, 50, 0.7, 1, 0);
	fpadIncl->Draw();
	TPad *p = fpadIncl->GetPad(1); //upper pad
	p->SetTickx();
	p->SetLogy(1);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	hfinal->GetXaxis()->SetRangeUser(0, 150);
	hset(*hfinal, "#it{p}_{T}^{jet}", "d^{2}#it{#sigma}/d#it{p}_{T}^{jet}d#it{#eta} (mb #it{c}/GeV)", 0.9, 1.1, 0.07, 0.06, 0.01, 0.001, 0.04, 0.05, 510, 510);
	hfinal -> SetMarkerStyle (20);
	hfinal -> SetMarkerColor (1);
	hfinal -> SetLineColor (1);
	if (scale) hfinal -> Scale(1./noe*0.5/0.3,"width");
	hfinal -> Draw();
	dpt->GetXaxis()->SetRangeUser(0, 150);
	if (scale) dpt -> Scale(1./noe*0.5/0.3,"width");
	dpt -> SetMarkerColor(1);
	dpt -> SetLineColor(1);
	dpt -> SetMarkerStyle(24);
	dpt -> Draw("same");

	auto lg = new TLegend(0.578947, 0.634783, 0.95933, 0.933913);
	lg->SetTextSize(0.0521739);
	lg->SetBorderSize(0);
	lg->AddEntry((TObject *)0x0, "Full jet", "");
	lg->AddEntry((TObject *)0x0, "#it{R} = 0.4, Anti-kt", "");
	lg->AddEntry((TObject *)0x0, "|#it{#eta}_{jet}|<0.3", "");
	lg->Draw();

	lg = new TLegend(0.26555,0.711304,0.569378,0.930435);
	lg->SetTextSize(0.0521739);
	lg->SetBorderSize(0);
	lg->AddEntry(hfinal, "Corrected data","lp");
	lg->AddEntry(dpt, "Raw data");
	lg->Draw();

	p = fpadIncl->GetPad(2); //upper pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	auto r = (TH1D *) dpt->Clone();
	r->Divide(hfinal);
	hset(*r, "#it{p}_{T}^{full jet} (GeV/#it{c})", "Ratio (raw/final)", 1.0, 0.8, 0.08, 0.09, 0.01, 0.001, 0.07, 0.06, 10, 510);
	r->Draw();
	r->Draw("same");
	return hfinal;

	/*
	pt->Scale(0.5/0.3, "width");
	tpt->Scale(0.5/0.3, "width");
	tptm->Scale(0.5/0.3, "width");
	ptc->Scale(1., "width");
	tptc->Scale(1., "width");

	auto lg = new TLegend(0.578947, 0.634783, 0.95933, 0.933913);
	lg->SetTextSize(0.0521739);
	lg->SetBorderSize(0);
	lg->AddEntry((TObject *)0x0, "Full jet", "");
	lg->AddEntry((TObject *)0x0, "R=0.4, Anti-kt", "");
	lg->AddEntry((TObject *)0x0, "|#it{#eta}_{jet}|<0.3", "");
	//------------------------------------------------------------


	pt->SetMarkerStyle(20);
	pt->SetMarkerColor(1);
	hset(*pt, "#it{p}_{T}^{full jet}", "d^{2}#it{#sigma}/d#it{p}_{T}^{jet}d#it{#eta} (mb #it{c}/GeV)", 0.9, 1.1, 0.07, 0.06, 0.01, 0.001, 0.04, 0.05, 510, 510);
	pt->Draw();
	tpt->SetMarkerStyle(24);
	tpt->SetMarkerColor(2);
	tpt->Draw("same");

	ptc->SetMarkerStyle(21);
	ptc->Draw("same");
	tptc->SetMarkerStyle(25);
	tptc->SetMarkerColor(2);

	tptc->Draw("same");

	lg->Draw();

	lg = new TLegend(0.186603, 0.0330435, 0.490431, 0.252174);
	lg->SetTextSize(0.0521739);
	lg->SetBorderSize(0);
	lg->AddEntry(pt, "rec full jet", "lp");
	lg->AddEntry(tpt, "true full jet", "lp");
	lg->AddEntry(ptc, "rec ch jet", "lp");
	lg->AddEntry(tptc, "true ch jet", "lp");
	lg->Draw();

	

	new TCanvas;
	jes->Draw("colz");

	Double1D binpt = {20, 22, 24, 26 ,28, 30 ,32, 36, 40, 42,44,46,48, 50,52,54,56,58, 60,65,70, 80,90, 100,110, 120, 130,140,150,160, 200,300};
	Double_t xbins[] = {20, 22, 24, 26 ,28, 30 ,32, 36, 40, 42,44,46,48, 50,52,54,56,58, 60,65,70, 80,90, 100,110, 120, 130,140,150,160, 200,300};
	auto *hjes = new TH1D("hjes", "hjes", binpt.size()-1, xbins);
	auto *hjesrms = new TH1D("hjesrms", "hjesrms", binpt.size()-1, xbins);
	vector<TH1D *> H;
	for (auto i = 0; i < binpt.size() - 1; i++)
	{
		auto h = jes->ProjectionY(Form("proy%d", i), jes->GetXaxis()->FindBin(binpt[i] + 0.0001), jes->GetXaxis()->FindBin(binpt[i + 1] - 0.0001), "e");
		hjes->SetBinContent(i + 1, h->GetMean());
		hjes->SetBinError(i + 1, h->GetMeanError());
		hjesrms->SetBinContent(i + 1, h->GetRMS());
		hjesrms->SetBinError(i + 1, h->GetRMSError());
		H.push_back(h);
	}
	new TCanvas;
	hjes->SetMarkerColor(1);
	hjes->SetLineColor(1);
	hjes->Draw("e1");

	new TCanvas;
	hjesrms->SetMarkerColor(2);
	hjesrms->SetLineColor(2);
	hjesrms->Draw("e1");

	fpadIncl = new Filipad2(1, 2, 0.5, 100, 50, 0.7, 1, 1);
	fpadIncl->Draw();
	p = fpadIncl->GetPad(1); //upper pad
	p->SetTickx();
	p->SetLogy(1);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	tpt->GetXaxis()->SetRangeUser(0, 300);
	tptm->GetXaxis()->SetRangeUser(0, 300);
	tpt->SetMarkerStyle(20);
	tpt->SetMarkerColor(1);
	hset(*tpt, "#it{p}_{T}^{full jet}", "d^{2}#it{#sigma}/d#it{p}_{T}^{jet}d#it{#eta} (mb #it{c}/GeV)", 0.9, 1.1, 0.07, 0.06, 0.01, 0.001, 0.04, 0.05, 510, 510);
	tpt->Draw();
	tptm->SetMarkerStyle(24);
	tptm->SetMarkerColor(2);
	tptm->Draw("same");

	p = fpadIncl->GetPad(2); //upper pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	auto mateff = (TH1D *)tptm->Clone();
	mateff->Divide(tpt);
	hset(*mateff, "#it{p}_{T}^{full jet} (GeV/#it{c})", "Matching eff.", 1.0, 0.8, 0.08, 0.09, 0.01, 0.001, 0.07, 0.06, 10, 510);
	mateff->Draw();

	fpadIncl = new Filipad2(1, 2, 0.5, 100, 50, 0.7, 1, 2);
	fpadIncl->Draw();
	p = fpadIncl->GetPad(1); //upper pad
	p->SetTickx();
	p->SetLogy(1);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	tm->GetXaxis()->SetRangeUser(0, 10);
	m->GetXaxis()->SetRangeUser(0, 10);
	m->SetMarkerStyle(20);
	m->SetMarkerColor(1);
	hset(*tm, "#it{p}_{T}^{full jet}", "d^{2}#it{#sigma}/d#it{p}_{T}^{jet}d#it{#eta} (mb #it{c}/GeV)", 0.9, 1.1, 0.07, 0.06, 0.01, 0.001, 0.04, 0.05, 510, 510);
	m->Draw();
	tm->SetMarkerStyle(24);
	tm->SetMarkerColor(2);
	tm->Draw("same");

	p = fpadIncl->GetPad(2); //upper pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	auto feff = (TH1D *)m->Clone();
	feff->Divide(tm);
	hset(*feff, "# jets", "Finding eff.", 1.0, 0.8, 0.08, 0.09, 0.01, 0.001, 0.07, 0.06, 10, 510);
	feff->Draw();



	TFile *f = new TFile("../histograms.root", "recreate");
	r->SetName("KinematicEfficiencyFullJet");
	r->Write();
	rc->SetName("KinematicEfficiencyChJet");
	rc->Write();
	mateff->SetName("MatchingEfficiency");
	mateff->Write();
	feff->SetName("FulljetFindingEfficiency");
	hjes->SetName("JES");
	hjes->Write();
	hjesrms->SetName("JER");
	hjesrms->Write();

	for (auto i = 0; i < binpt.size() - 1; i++){
		auto h = H.at(i);
		h->SetName(Form("pt%d_%d",Int_t(binpt.at(i)),Int_t(binpt.at(i+1))));
		h->Write();
	}
	f->Close();
	*/
}

TH1D *getHisto(TFile* fIn, TString sHistoName) {
    if(fIn==0) ErrorExit("getHisto: Input file pointer null!");
    TH1D *hReturn;
    hReturn=(TH1D*)fIn->Get(sHistoName.Data());
    if(hReturn==0) ErrorExit(Form("getHisto: histo=%s pointer null! Check histo name.",sHistoName.Data()));
    return hReturn;
}
//
// Calculate ratio of two TGraphErrors. If bin size is different, gr2 will be evaluated by linear interpolation.
TGraphErrors *CalculateRatio( TGraphErrors* gr1, TGraphErrors* gr2, double xshift=0.0, bool invert=false){
    int NPoint = gr1->GetN();
    TGraphErrors *gr_ratio = new TGraphErrors( NPoint);
    TGraph ger( gr2->GetN(), gr2->GetX(), gr2->GetEY() ); // << Err estimation of gr2
    double x, y1, ey1, y2, ey2, ex, ratio;
    cout << "x, y1, y2, ratio, ey1, ey2, error" << endl;//ey1/ey2 << endl;
    for(int i=0; i<NPoint; i++){
        x = gr1->GetX()[i];
        y1 = gr1->GetY()[i];
        ey1 = gr1->GetEY()[i];
        y2 = gr2->Eval(x);
        ey2 = ger.Eval(x);
        ex = gr1->GetErrorX(i);
        if (y1==0 || y2==0) continue; // To prevent division by zero
        if(invert) ratio = TMath::Abs( y1 / y2 );
        else ratio = TMath::Abs( y2 / y1 );
        gr_ratio->SetPoint( i, x+xshift, ratio);
        gr_ratio->SetPointError( i, ex, ratio * TMath::Sqrt( ey1*ey1/y1/y1 + ey2*ey2/y2/y2 )); 
        cout << x <<"\t"<< y1 <<"\t"<< y2 <<"\t"<< ratio << "\t" << ey1 << "\t" << ey2 << "\t" << ratio * TMath::Sqrt( ey1*ey1/y1/y1 + ey2*ey2/y2/y2 ) << endl;//ey1/ey2 << endl;
    }   
    return gr_ratio;
}
