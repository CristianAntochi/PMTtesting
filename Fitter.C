// Root Macro to fit gain spectra
//
// Bart
#include <utility>
#include <math.h>

Double_t fitGainfunc(Double_t *x, Double_t *par){
  Double_t fitval = par[0]*ROOT::Math::normal_pdf(x[0]-par[1], par[2],0)+
    par[3]*ROOT::Math::normal_pdf(x[0]-2*par[1],sqrt(2)*par[2],0)+
    par[4]*exp(-(x[0]-par[9])/par[5])+
    par[6]*ROOT::Math::normal_pdf(x[0]-par[7],par[8],0)+
    par[10]*exp(-x[0]/par[11])+
    par[12]*ROOT::Math::normal_pdf(x[0]-3*par[1],sqrt(3)*par[2],0)
    ;
  return fitval;
  
}

const pair<float, float> fitFile(const TString &fileName,int b,int channel,int v);

//void fitAll(const TString &led) {
  void    fitAll(float led) {
    /* fits all files matching: processed_wf[0-7]_Board[0-1]_Gain_[hv]V_[nameEnd]
     * in path
     * Make csv file with columns: board, ch, hv, mean, sigma
     */
    TString path = "/home/cristian/Suxess_files/data/run5/liquid/gain_study/Alfredo_method/600-300-600_i3-10/";  //"/media/suxess/DATADISK2/daq_data/run4/vacuum/gain171213/processed_test/";
    //TString outName = "Gain171213_fits_cuts.csv";
    TString outName = "Test_fits_cuts.csv";
    ofstream outFile;
    outFile.open (outName.Data(), ios_base::app);
   
    for (int board = 0; board < 2; ++board) {
        for (int ch = 0; ch < 7; ++ch) {
	  // for (int hv=1300; hv<=1350; hv=hv+50) {
            TString name = TString::Format("processed_wf%d_Gain%d_Test_standard_600_300_600.root", ch, board);//,led);


	   int hv=1400;
            const pair<float, float> fit = fitFile(path + name,board,ch,hv);
            outFile << board << "," << ch << ","/* << hv*/ << "," << fit.first << "," << fit.second << endl;
	    //	}
        }
    }

    outFile.close();
}

const pair<float, float> fitFile(const TString &fileName, int b, int channel, int v) {
    /* Fit a single file (single channel)
     * return a pair<mean, sigma> of the fit
     */
    TFile *f = TFile::Open(fileName);
    TTree *tree = (TTree*)f->Get("T1");
    TString title = TString::Format("Board%d  Channel %d voltage %d", b,channel, v);
    TCanvas* c = new TCanvas(title,title);

    
    TF1 *func = new TF1("fitGainfunc",fitGainfunc,50,800,13);//setting the fitting function

    
    c->cd();
    c->SetLogy();
    
    TH1D* h = new TH1D("h","h", 1000, 0, 1999);

    tree->Draw("Area>>h","Entropy<0.8");
    
    // Find 1pe peak maximum
    h->GetXaxis()->SetRange(70, 1000);
    double width = 40;
    if(v>1370){
      //h->GetXaxis()->SetRange(150, 700);
      width=70;
    }             
    if(v>1420){
      //h->GetXaxis()->SetRange(150, 700);
      width=100;
    }
    if(v>1520){
      //h->GetXaxis()->SetRange(300, 1300);
      width=300;
    }
    int binmax = h->GetMaximumBin();
    double mean = h->GetXaxis()->GetBinCenter(binmax);
    h->GetXaxis()->SetRange(5,800);

    TFitResultPtr res = h->Fit("gaus", "S", "", mean - width, mean + width);
    //double width = 600;
    //TFitResultPtr res = h->Fit("gaus", "S", "", mean - (0.5*v-width), mean + (0.5*v-width));
    
    

    // double chi2 = res->Chi2();
    // double ndf = res->Ndf();
    const double *pars = res->GetParams();
    const double *errs = res->GetErrors();
    
    func->SetParLimits(0,100,500000);
    func->SetParameter(0,pars[0]);
    func->SetParameter(1,pars[1]);
    func->SetParLimits(1,pars[1]-20,pars[1]+20);
    func->SetParameter(2,pars[2]);
    func->SetParLimits(2,pars[2]-30,pars[2]+60);
    func->SetParLimits(3,10,50000);
    func->SetParameter(3,1000);
    func->SetParLimits(4,100,500000);
    func->SetParameter(4,50000);
    func->SetParLimits(5,10,500);
    func->SetParameter(5,55);
    func->SetParLimits(6,10000,500000);
    func->SetParameter(6,50000);
    func->SetParLimits(7,-10,40);
    func->SetParameter(7,5);
    func->SetParLimits(8,0,50);
    func->SetParameter(8,10);
    func->SetParameter(9,50);
    func->SetParLimits(10,0,5000);
    func->SetParameter(10,500);
    func->SetParLimits(11,2,100);
    func->SetParameter(11,30);

    h->GetXaxis()->SetRange(5,800);
    TFitResultPtr res1 = h->Fit("fitGainfunc","S,M,E","");
    const double *pars1 = res1->GetParams();
    const double *errs1 = res1->GetErrors();
    
    func->SetParameter(0,pars1[0]);
    func->SetParameter(1,pars1[1]);
    func->SetParameter(2,pars1[2]);
    func->SetParameter(3,pars1[3]);
    func->SetParameter(4,25000);
    func->SetParameter(5,45);
    func->SetParameter(6,pars1[6]);
    func->SetParameter(7,pars1[7]);
    func->SetParameter(8,pars1[8]);
    func->SetParameter(9,pars1[9]);
    func->SetParameter(10,pars1[10]);
    func->SetParameter(11,pars1[11]);

    h->GetXaxis()->SetRange(5,800);
    TFitResultPtr res2 = h->Fit("fitGainfunc","S,M,E","");
    const double *pars2 = res2->GetParams();
    const double *errs2 = res2->GetErrors();
    
    double chi2 = res2->Chi2();
    double ndf = res2->Ndf();
    c->Update();
    
    /*if (pars[1] < 0) {
        res = h->Fit("gaus", "S", "", mean , mean + 70);
        chi2 = res->Chi2();
        ndf = res->Ndf();
        pars = res->GetParams();
        errs = res->GetErrors();
    }*/

    cout << "Fitted File: " << fileName << endl;
    cout << "chi2/ndf " << chi2 << "/" << ndf << endl;
    cout << "mean  sigma" << endl;
    cout << pars2[1] << "  " << pars2[2] << endl;
    return make_pair(pars2[1], pars2[2]);
}
