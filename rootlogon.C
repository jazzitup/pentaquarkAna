{
   TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");

   // from ROOT plain style
   myStyle->SetCanvasBorderMode(0);
   myStyle->SetPadBorderMode(0);
   myStyle->SetPadColor(0);
   myStyle->SetCanvasColor(0);
   myStyle->SetTitleColor(1);
   myStyle->SetStatColor(0);

   myStyle->SetLabelSize(0.03,"xyz"); // size of axis values


   // default canvas positioning
   myStyle->SetCanvasDefX(900);
   myStyle->SetCanvasDefY(20);
   myStyle->SetCanvasDefH(550);
   myStyle->SetCanvasDefW(540);

   myStyle->SetPadBottomMargin(0.1);
   myStyle->SetPadTopMargin(0.1);
   myStyle->SetPadLeftMargin(0.1);
   myStyle->SetPadRightMargin(0.1);
   myStyle->SetPadTickX(1);
   myStyle->SetPadTickY(1);
   myStyle->SetFrameBorderMode(0);

   // Din letter
   myStyle->SetPaperSize(21, 28);
   

   myStyle->SetOptStat(0);
   //   myStyle->SetOptStat(111111);// Show overflow and underflow as well
   //   myStyle->SetOptFit(1011);
   myStyle->SetPalette(1); 

   // apply the new style
   gROOT->SetStyle("MyStyle"); //uncomment to set this style
   gROOT->ForceStyle(); // use this style, not the one saved in root files


    gStyle->SetTitleSize(20, "X"); // 0.05 corresponds to approximately font size 20 for X-axis
    gStyle->SetTitleSize(20, "Y"); // 0.05 corresponds to approximately font size 20 for Y-axis
    
    // Optional: Set default title font (e.g., Times New Roman)
    gStyle->SetTitleFont(43, "X");
    gStyle->SetTitleFont(43, "Y");

    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.13);
    
    gStyle->SetLabelSize(0.035, "X");
    gStyle->SetLabelSize(0.035, "Y");

    gStyle->SetMarkerSize(1);
    gStyle->SetMarkerStyle(20);
    gStyle->SetTitleOffset(.9, "X");
    gStyle->SetTitleOffset(.9, "Y");

    // gStyle->SetTitleAlign(23);
   
   printf("\n Beginning new ROOT session with private TStyle \n");

}
