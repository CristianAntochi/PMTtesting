Double_t fitf(Double_t *x, Double_t *par)
{
   Double_t arg = 0;
   if (par[2] != 0) arg = (x[0] - par[1])/par[2];

   Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg)+par[3]*TMath::Exp(-(x[0]-2*par[1])*(x[0]-2*par[1])/(4*par[2]*par[2]))+par[4]*TMath::Exp(-(x[0]-3*par[1])*(x[0]-3*par[1])/(6*par[2]*par[2]));
   return fitval;
}


