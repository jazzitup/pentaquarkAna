#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
//#ifndef COMMONUTILITY_Yongsun_H
//#define COMMONUTILITY_Yongsun_H
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGaxis.h>

void easyRange(TH1 *h=0, double ratio=1.6, double minimum=0.001)
{
  double maxBinContent = h->GetMaximum();
  h->GetYaxis()->SetRangeUser(minimum, ratio * maxBinContent);


}
void handsomeTH1( TH1 *a=0, int color=1)
{
  //   a->GetYaxis()->SetTitleOffset(1.25);
   a->GetXaxis()->CenterTitle();
   a->GetYaxis()->CenterTitle();
   a->SetLineColor(color);
   a->SetMarkerColor(color);
   a->Sumw2();
}

void handsomeTH2( TH2 *a=0)
{
   a->GetYaxis()->SetTitleOffset(1.25);
   a->GetXaxis()->CenterTitle();
   a->GetYaxis()->CenterTitle();
}
void handsomeTH3( TH3 *a=0)
{
   a->GetYaxis()->SetTitleOffset(1.25);
   a->GetXaxis()->CenterTitle();
   a->GetYaxis()->CenterTitle();
}


void easyLeg( TLegend *l=0, char *title="")
{
  l->SetHeader(title,"C"); // option "C" allows to center the header
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  //  l->SetTextSize(20);
}
