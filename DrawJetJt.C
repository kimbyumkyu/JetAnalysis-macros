#include "BSHelper.cxx"
#include "JetCommon.h"
#include "Filipad2.h"
#include "TUnfold.h"

const float inf = 1e20;
using namespace std;
Int_t colors[] = {kRed, kOrange, kGreen, kGreen + 3, kBlue, kBlue + 2, kViolet, kCyan + 2, kAzure, kMagenta, kMagenta + 3, kYellow + 3, kYellow + 2, kSpring - 7};

TH1D *getHisto(TFile* fIn, TString sHistoName);
TGraphErrors *CalculateRatio( TGraphErrors* gr1, TGraphErrors* gr2, double xshift=0.0, bool invert=false);
TH1D *FullJetPt(TString ndata, TString nmc, Int_t type); // type : 1 inclusive, 2 perpendicular jet

void DrawJetJt(TString ndata = "JtLHC13cAOD", TString nmc = "JtLHC13cMCAOD")
{

	auto hfulljetpt =  FullJetPt( "JtLHC13cAOD", "JtLHC13cMCAOD", 1);

	gStyle->SetErrorX(0.001);

	Double1D binCent = {0, 100}; // centbin merging
	Double1D binz = {0, 1};		 // inclusive for z
	const int npthardbin = 10;	 // LHC13b4_plus has 10 pt hard bins

	auto datalist = LoadJetResultList (ndata.Data(), ndata.Data());
	auto DJT = BSTHnSparseHelper::Load( "hJetJt", datalist ); // data jet pt
	auto NOE = (TH1D*) datalist -> FindObject("hEventNumbers");
	Double_t noe = NOE->GetBinContent(3);

	// Basic idea : corrected jet jt = DJT* (TJT/MJT), 
	auto mclist = LoadJetResultList (nmc.Data(), nmc.Data());
	auto TJT = BSTHnSparseHelper::Load( "hJetJtTruth", mclist ); // MC true jet jt
	auto MJT = BSTHnSparseHelper::Load( "hJetJt", mclist ); // MC rec jet jt
	auto RES2 = BSTHnSparseHelper::Load( "hJetJtRes", mclist ); // Response matrix 
	auto FAKE = BSTHnSparseHelper::Load( "hJetJtFake", mclist ); // Fake jet jt in response matrix 
	auto MISS = BSTHnSparseHelper::Load( "hJetJtMiss", mclist ); // Missing jet jt in response matrix 
	auto HN = BSTHnSparseHelper::Load( "hPtHardBin", mclist ); // # of events for each pt hard bins

	// Set centrality range {0,100} 
	DJT.SetBin("Cent",binCent);
	TJT.SetBin("Cent",binCent);
	MJT.SetBin("Cent",binCent);
	RES2.SetBin("Cent",binCent);
	FAKE.SetBin("Cent",binCent);
	MISS.SetBin("Cent",binCent);


	DJT.SetBin("zbin",binz);
	TJT.SetBin("zbin",binz);
	MJT.SetBin("zbin",binz);
	RES2.SetBin("zbin",binz);
	FAKE.SetBin("zbin",binz);
	MISS.SetBin("zbin",binz);



	HN.SetBin("Cent",binCent);
	auto hn = HN.GetTH1("hn",1,{-1,-1}); // # of events of each pt hard bin by projection
	auto djt = DJT.GetTH2("hdpt",1,3,{-1,-1,-1,-1, 1}); // # data jt weighted distribution 

	//-----------------------------------------------------------
    //pT hard bin scaling for all MC histograms
	TH2D *tjt = nullptr;
	TH2D *mjt = nullptr;
	THnSparse *RESTEMP = nullptr;
	TH2D *fake = nullptr;
	TH2D *miss = nullptr;
	for (auto i = 0; i < npthardbin; i++)
	{
		auto h = TJT.GetTH2("htjt",1,3, {-1,-1,-1,-1,i+1});
		h->Scale(1. / hn->GetBinContent(i + 1), "");
		if (i == 0)
			tjt = (TH2D *) h->Clone();
		else
			tjt->Add(h);
		
		h = MJT.GetTH2("hmjt",1, 3, {-1,-1,-1,-1,i+1});
		h->Scale(1. / hn->GetBinContent(i + 1), "");
		if (i == 0)
			mjt = (TH2D *)h->Clone();
		else
			mjt->Add(h);

		h = FAKE.GetTH2("hfake",1, 3, {-1,-1,-1,-1,i+1});
		h->Scale(1. / hn->GetBinContent(i + 1), "");
		if (i == 0)
			fake = (TH2D *)h->Clone();
		else
			fake->Add(h);
			
		h = MISS.GetTH2("hmiss",1,3, {-1,-1,-1,-1,i+1});
		h->Scale(1. / hn->GetBinContent(i + 1), "");
		if (i == 0)
			miss = (TH2D *)h->Clone();
		else
			miss->Add(h);
		
		auto H = RES2.GetTHnSparse("HRES",{1,2,4,5}, {-1,-1,-1,-1,-1,-1,i+1});
		H->Scale(1. / hn->GetBinContent(i + 1));
		if (i == 0)
			RESTEMP = (THnSparse *)H->Clone();
		else
			RESTEMP->Add(H);

	}

	auto RES = BSTHnSparseHelper(RESTEMP);
	RES.PrintAxis("all");
    RooUnfoldResponse response (mjt, tjt);

	for (auto i = 0; i <= RESTEMP->GetAxis(0)->GetNbins() + 1; i++)
	{ //recjetpt
		for (auto j = 0; j <= RESTEMP->GetAxis(1)->GetNbins() + 1; j++)
		{ //truthjetpt
			for (auto k = 0; k <= RESTEMP->GetAxis(2)->GetNbins() + 1; k++)
			{ //rec jt
				for (auto l = 0; l <= RESTEMP->GetAxis(3)->GetNbins() + 1; l++)
				{ //truth jt
					Int_t bins[4] = {i, j, k, l};
					Double_t bincenters[4] = {
						RESTEMP->GetAxis(0)->GetBinCenter(i),
						RESTEMP->GetAxis(1)->GetBinCenter(j),
						RESTEMP->GetAxis(2)->GetBinCenter(k),
						RESTEMP->GetAxis(3)->GetBinCenter(l)};
					response.Fill(bincenters[0], bincenters[2], bincenters[1], bincenters[3], RESTEMP->GetBinContent(bins));
				}
			}
		}
	}

	for (auto i = 0; i <= miss->GetNbinsX() + 1; i++)
	{ //invM
		for (auto j = 0; j <= miss->GetNbinsY() + 1; j++)
		{ //ptpair
			Double_t bincenx = miss->GetXaxis()->GetBinCenter(i);
			Double_t binceny = miss->GetYaxis()->GetBinCenter(j);
			Double_t bincont = miss->GetBinContent(i, j);
			response.Miss(bincenx, binceny, bincont);
		}
	}
	for (auto i = 0; i <= fake->GetNbinsX() + 1; i++)
	{ //invM
		for (auto j = 0; j <= fake->GetNbinsY() + 1; j++)
		{ //ptpair
			Double_t bincenx = fake->GetXaxis()->GetBinCenter(i);
			Double_t binceny = fake->GetYaxis()->GetBinCenter(j);
			Double_t bincont = fake->GetBinContent(i, j);
			response.Fake(bincenx, binceny, bincont);
		}
	}

    RooUnfoldBayes  unfold (&response, djt,4);
    auto hfinal= (TH2D*) unfold.Hreco();


	// Let's draw the results
	Filipad2 *fpad = new Filipad2(1, 2, 0.5, 100, 50, 0.7, 1, 1);
	fpad->Draw();
	TPad *p = fpad->GetPad(1); //upper pad
	p->SetTickx();
	p->SetLogy(1);
	p->SetLogx(1);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	auto hdummy = hfinal->ProjectionY("dum", 0, -1, "");
	hdummy->Reset();
	hdummy->GetXaxis()->SetRangeUser(0.01, 10);
	hdummy->Draw();

	Double1D binJetPt = {40, 60, 80, 100, 150}; // centbin merging
	for (auto i = 0; i < binJetPt.size()-1; i++){
		auto lbin = hfinal->GetXaxis()->FindBin(binJetPt.at(i)+0.0001);
		auto rbin = hfinal->GetXaxis()->FindBin(binJetPt.at(i+1)-0.0001);
		auto h = hfinal->ProjectionY(Form("hjt%d", i), lbin, rbin, "e");
		auto lbin2 = hfulljetpt->GetXaxis()->FindBin(binJetPt.at(i)+0.0001);
		auto rbin2 = hfulljetpt->GetXaxis()->FindBin(binJetPt.at(i+1)-0.0001);
		auto numberofjets = hfulljetpt->Integral(lbin2, rbin2);
		h->Scale(1./numberofjets, "width");
		auto g = new TGraph;
		for (auto j = 1; j <= h->GetNbinsX(); j++ ){
			g->SetPoint(g->GetN(), h->GetXaxis()->GetBinCenter(j), h->GetBinContent(j));
		}
		g->SetMarkerStyle(20);
		g->SetMarkerColor(colors[i]);
		g->Draw("PZsame");
	}

	/*

	hfinal->GetXaxis()->SetRangeUser(0, 150);
	hset(*hfinal, "#it{j}_{T}^{jet}", "d^{2}#it{#sigma}/d#it{p}_{T}^{jet}d#it{#eta} (mb #it{c}/GeV)", 0.9, 1.1, 0.07, 0.06, 0.01, 0.001, 0.04, 0.05, 510, 510);
	hfinal -> SetMarkerStyle (20);
	hfinal -> SetMarkerColor (1);
	hfinal -> SetLineColor (1);
	hfinal -> Scale(1./noe*0.5/0.3,"width");
	hfinal -> Draw();
	djt->GetXaxis()->SetRangeUser(0, 150);
	djt -> Scale(1./noe*0.5/0.3,"width");
	djt -> SetMarkerColor(1);
	djt -> SetLineColor(1);
	djt -> SetMarkerStyle(24);
	djt -> Draw("same");

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
	lg->AddEntry(djt, "Raw data");
	lg->Draw();

	p = fpadIncl->GetPad(2); //upper pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	auto r = (TH1D *) djt->Clone();
	r->Divide(hfinal);
	hset(*r, "#it{p}_{T}^{full jet} (GeV/#it{c})", "Ratio (raw/final)", 1.0, 0.8, 0.08, 0.09, 0.01, 0.001, 0.07, 0.06, 10, 510);
	r->Draw();
	r->Draw("same");
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

TH1D* FullJetPt(TString ndata , TString nmc, Int_t type) // type : 1 inclusive, 2 perpendicular jet
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
	auto dpt = DPT.GetTH1("hdpt",1,{-1,-1, 1, type}); // # of events of each pt hard bin by projection

	//-----------------------------------------------------------
    //pT hard bin scaling for all MC histograms
	TH1D *tpt = nullptr;
	TH1D *mpt = nullptr;
	TH2D *res = nullptr;
	TH1D *fake = nullptr;
	TH1D *miss = nullptr;
	for (auto i = 0; i < npthardbin; i++)
	{
		auto h = TPT.GetTH1("htpt",1, {-1,-1,i+1,type});
		h->Scale(1. / hn->GetBinContent(i + 1), "");
		if (i == 0)
			tpt = (TH1D *) h->Clone();
		else
			tpt->Add(h);
		
		h = MPT.GetTH1("hmpt",1, {-1,-1,i+1, type});
		h->Scale(1. / hn->GetBinContent(i + 1), "");
		if (i == 0)
			mpt = (TH1D *)h->Clone();
		else
			mpt->Add(h);

		h = FAKE.GetTH1("hfake",1, {-1,-1,i+1, type});
		h->Scale(1. / hn->GetBinContent(i + 1), "");
		if (i == 0)
			fake = (TH1D *)h->Clone();
		else
			fake->Add(h);
			
		h = MISS.GetTH1("hmiss",1, {-1,-1,i+1, type});
		h->Scale(1. / hn->GetBinContent(i + 1), "");
		if (i == 0)
			miss = (TH1D *)h->Clone();
		else
			miss->Add(h);
		
		auto h2 = RES.GetTH2("hres",1,2, {-1,-1,-1,i+1, type});
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
	auto hfinalbeforescale = (TH1D *)hfinal->Clone();
	hfinal->Scale(1. / noe * 0.5 / 0.3, "width");
	hfinal -> Draw();
	dpt->GetXaxis()->SetRangeUser(0, 150);
	dpt -> Scale(1./noe*0.5/0.3,"width");
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
	return hfinalbeforescale;
}