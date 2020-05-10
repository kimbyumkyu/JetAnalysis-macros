#include "BSHelper.cxx"
#include "JetCommon.h"
#include "Filipad2.h"
#include "TUnfold.h"

const float inf = 1e20;
using namespace std;
Int_t colors[] = {kRed, kOrange, kGreen, kGreen + 3, kBlue, kBlue + 2, kViolet, kCyan + 2, kAzure, kMagenta, kMagenta + 3, kYellow + 3, kYellow + 2, kSpring - 7};

TH1D *getHisto(TFile* fIn, TString sHistoName);
TGraphErrors *CalculateRatio( TGraphErrors* gr1, TGraphErrors* gr2, double xshift=0.0, bool invert=false);

void RawJetPt(Bool_t ispp = 1 )
{

	gStyle->SetErrorX(0.001);
	TString fdata("TestLHC13dAOD");

	Double1D binCent = {0, 100}; // centbin merging
	Double1D binJetPt = {40, 60, 80, 100, 150}; // centbin merging
	Double1D binz = {0, 100};

	auto datalist = LoadJetResultList (fdata.Data(), fdata.Data());
	auto JT = BSTHnSparseHelper::Load( "hJetJtWeight", datalist );
	JT.SetBin("Cent",binCent);
	JT.SetBin("binjetpt",binJetPt);
	JT.SetBin("zbin",binz);
	JT.PrintAxis("all");

	//-----------------------------------------------------------
    //Draw full jet pt distribution
	auto PT = BSTHnSparseHelper::Load( "hJetPtFullJet", datalist );
	PT.SetBin("Cent", binCent);
	auto pt = PT.GetTH1("fulljetpt", 1, {1, -1, 1}); // pthardbin = 1 always for data
	if(1){
		auto canvas = new TCanvas("canvas","canvas",1280,800);
		gPad->SetLogy(1);
		auto ptwidth = (TH1D *)pt->Clone();
		ptwidth->Scale(1, "width");
		ptwidth->SetMarkerStyle(20);
		ptwidth->GetXaxis()->SetRangeUser(0, 200);
		ptwidth->Draw();
	}
	//------------------------------------------------------------

	Filipad2 *fpadIncl = new Filipad2(1, 2, 0.5, 100, 50, 0.7, 1, 0);
	fpadIncl->Draw();
	TPad *p = fpadIncl->GetPad(1); //upper pad
	p->SetTickx();
	p->SetLogy(1);
	p->SetLogx(1);
	;
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TH1D *hdummy = new TH1D("dummy", "dummy", 40000, 0, 4);
	hdummy->GetXaxis()->SetRangeUser(0.08, 4);
	hset(*hdummy, "#it{j}_{T}", "d^{2}#it{#sigma}/d#it{p}_{T}^{jet}d#it{#eta} (mb #it{c}/GeV)", 0.9, 1.1, 0.07, 0.06, 0.01, 0.001, 0.04, 0.05, 510, 510);
	hdummy->Draw();

	TLegend *lg = new TLegend(0.186603, 0.0330435, 0.490431, 0.252174);
	lg->SetTextSize(0.0521739);
	lg->SetBorderSize(0);

	//for (auto i = 1; i <= binJetPt.size()-1; i++)
	for (auto i = 1; i <2; i++)
	{
		auto jt = JT.GetTH1(Form("h%d", i), 3, {1, i, 1, -1, 1});
		jt->SetMarkerStyle(20);
		jt->SetMarkerColor(colors[i-1]);
		Int_t lbin = pt->GetXaxis()->FindBin(binJetPt.at(i - 1) + 0.00001);
		Int_t rbin = pt->GetXaxis()->FindBin(binJetPt.at(i) - 0.00001);
		Double_t njets = pt->Integral(lbin, rbin);
		cout << i <<" "<< binJetPt.at(i - 1) << " " << binJetPt.at(i) << endl;
		jt->Scale(1. / njets, "width"); // bin width scaling
		//jt->Draw("hist same");
		auto g = new TGraph;
		//g->SetPoint(g->GetN(), 0, 0);
		for (auto j = 1; j <= jt->GetNbinsX(); j++ ){
			g->SetPoint(g->GetN(), jt->GetXaxis()->GetBinCenter(j), jt->GetBinContent(j));
		}
		g->SetMarkerStyle(20);
		g->Draw("PZsame");
	}
	//lg->Draw();

	lg = new TLegend(0.578947, 0.634783, 0.95933, 0.933913);
	lg->SetTextSize(0.0521739);
	lg->SetBorderSize(0);
	lg->AddEntry((TObject *)0x0, "Charged jet", "");
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

	//r->Divide(hd);
	//hset(*r, "#it{p}_{T}^{jet} (GeV/#it{c})", "Ratio (Oskari/Beomkyu)", 1.0, 0.8, 0.08, 0.09, 0.01, 0.001, 0.07, 0.06, 10, 510);



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
