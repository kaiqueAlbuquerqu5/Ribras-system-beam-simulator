#include<vector>

void help(){
	cout << "*********************************************" << endl;
	cout << "This program has the following functions:    " << endl;
	cout << "---> biparametric()                          " << endl;
	cout << "---> biparametric_contaminants()             " << endl;
	cout << "---> AngularHistogram(const char *file_name) " << endl;
	cout << "---> EnergyHistogram(const char *file_name)  " << endl;
	cout << "---> ReactionPoints()                        " << endl;
	cout << "---> DetectionPoints()                       " << endl;
	cout << "---> DetectionAndReactionPoints()            " << endl;
	cout << "*********************************************" << endl;
}

void biparametric_contaminants(){
   cout << "********************************************************************" << endl;
   cout << "              This function makes biparametric charts.              " << endl;
   Double_t Energy;
   Double_t Del_Energy;
   Double_t Max_E_x = 0, Max_E_y = 0;
   
   FILE *fp = fopen("Contaminants/BipaAllParticles.dat","r");
   
   vector<vector<double>> DATA;
   char line[80];
   while(fgets(line,80,fp)) {
      sscanf(&line[0],"%lf %lf", &Energy, &Del_Energy);
      if(Energy > Max_E_x){
         Max_E_x = Energy;
      }
      if(Del_Energy > Max_E_y){
         Max_E_y = Del_Energy;
      }   
      DATA.push_back({Energy,Del_Energy});
   } 
   Max_E_x = 10*(1+floor(Max_E_x/10));
   Max_E_y = 10*(1+floor(Max_E_y/10));
   
   TCanvas *Biparametric = new TCanvas("Biparametric Contaminants");
   TH2F *bipa = new TH2F("bipa",";E (MeV);#Delta E (MeV)",1500,0,Max_E_x,1500,0,Max_E_y);
   bipa->SetTitleSize(0.045,"X");
   bipa->SetTitleSize(0.045,"Y");
   for(int i = 0; i < DATA.size(); ++i){
      bipa->Fill(DATA[i][0],DATA[i][1]);
   }
  
   bipa->Draw("colz");
   Biparametric->Modified();
   Biparametric->Update();
//   gStyle->SetPalette(kRainbow);
   bipa->SetStats(0);
   fclose(fp);
   cout << "********************************************************************" << endl;
}

void biparametric(){
   cout << "********************************************************************" << endl;
   cout << "              This function makes biparametric charts.              " << endl;
   Double_t Energy;
   Double_t Del_Energy;
   Double_t Max_E_x = 0, Max_E_y = 0;
   
   FILE *fp = fopen("BiParametric.dat","r");
   
   vector<vector<double>> DATA;
   char line[80];
   while(fgets(line,80,fp)) {
      sscanf(&line[0],"%lf %lf", &Energy, &Del_Energy);
      if(Energy > Max_E_x){
         Max_E_x = Energy;
      }
      if(Del_Energy > Max_E_y){
         Max_E_y = Del_Energy;
      }   
      DATA.push_back({Energy,Del_Energy});
   } 
   Max_E_x = 10*(1+floor(Max_E_x/10));
   Max_E_y = 10*(1+floor(Max_E_y/10));
   
   TCanvas *Biparametric = new TCanvas("Biparametric");
   TH2F *bipa = new TH2F("bipa",";E (MeV);#Delta E (MeV)",1500,0,Max_E_x,1500,0,Max_E_y);
   bipa->SetTitleSize(0.045,"X");
   bipa->SetTitleSize(0.045,"Y");
   for(int i = 0; i < DATA.size(); ++i){
      bipa->Fill(DATA[i][0],DATA[i][1]);
   }
  
   bipa->Draw("colz");
   Biparametric->Modified();
   Biparametric->Update();
   gStyle->SetPalette(kRainbow);
   bipa->SetStats(0);
   fclose(fp);
   cout << "********************************************************************" << endl;
}

void AngularHistogram(const char *file_name){
   cout << "********************************************************************" << endl;
   cout << "    This function makes a histogram of an angular distribution.     " << endl;
   FILE *fp  = fopen(file_name,"r");
   
   Double_t Max = 0, Min=10000, Binage = 100;
   vector<double> DATA;
   Double_t A;
   char line[80];
   while(fgets(line,80,fp)) {
      sscanf(&line[0],"%lf",&A);
      if(A > Max){
         Max = A;
      }
      if(A < Min){
         Min = A;
      }
      DATA.push_back(A);
   }
   Max = 10*(1+floor(Max/10));
   Min = 10*(floor(Min/10));
   
   char BinChoice;
   TCanvas *c = new TCanvas();
   while(true){
      TH1D *Ahistogram = new TH1D("Angular histogram","Angular distribution;Angle;Counts",Binage,Min,Max);
       
      for(Int_t i = 0; i < DATA.size(); ++i){
         Ahistogram->Fill(DATA[i]);
      }
        
      Ahistogram->Draw();
      c->Modified();
      c->Update();
      cout << endl << "Do you want to change the binage (y/n)? ";
      cin >> BinChoice;
      if(BinChoice == 'n'){
         break;
      }
      cout << "Type the binage: ";
      cin >> Binage;
   }
   cout << "********************************************************************" << endl;
}

void EnergyHistogram(const char *file_name){
   cout << "********************************************************************" << endl;
   cout << "     This function makes a histogram of an energy distribution.     " << endl;
   FILE *fp  = fopen(file_name,"r");
   
   Double_t Max = 0, Min=10000, Binage = 100;
   vector<double> DATA;
   Double_t E;
   char line[80];
   while(fgets(line,80,fp)) {
      sscanf(&line[0],"%lf",&E);
      if(E > Max){
         Max = E;
      }
      if(E < Min){
         Min = E;
      }
      DATA.push_back(E);
   }
   Max = 10*(1+floor(Max/10));
   Min = 10*(floor(fabs(Min)/10));
   
   char BinChoice, LimitsChoice;
   TCanvas *c = new TCanvas();
   while(true){
      TH1D *Ehistogram = new TH1D("Energy histogram","Energy Distribution;Energy;Counts",Binage,Min,Max);
       
      for(Int_t i = 0; i < DATA.size(); ++i){
         Ehistogram->Fill(DATA[i]);
      }
        
      Ehistogram->Draw();
      c->Modified();
      c->Update();
      cout << endl << "Do you want to change the binage (y/n)? ";
      cin >> BinChoice;
      cout << "Do you want to change the axes limits (y/n)? ";
      cin >> LimitsChoice;
      if(BinChoice == 'n' && LimitsChoice == 'n'){
         break;
      }
      if(BinChoice != 'n'){
         cout << endl << "Type the binage: ";
         cin >> Binage;
      }
      if(LimitsChoice != 'n'){
      	 cout << "Type the x minimum: ";
         cin >> Min;
         cout << "Type the x maximum: ";
         cin >> Max;
      }
   }
   cout << "********************************************************************" << endl;
}

void ReactionPoints(){
   cout << "********************************************************************" << endl;
   cout << "              This function plots the reaction points               " << endl;
   FILE *fp  = fopen("ReactionPoints.dat","r");
   Double_t x,y,z;
   vector<double> ChartLimits = {10000,0,10000,0,10000,0};
   
   vector<vector<double>> DATA;
   char line[80];
   while(fgets(line,80,fp)) {
      sscanf(&line[0],"%lf %lf %lf", &x, &y, &z);
      if(x < ChartLimits[0]){
         ChartLimits[0] = x;
      }
      if(x > ChartLimits[1]){
         ChartLimits[1] = x;
      }   
      if(y < ChartLimits[2]){
         ChartLimits[2] = y;
      }
      if(y > ChartLimits[3]){
         ChartLimits[3] = y;
      }
      if(z < ChartLimits[4]){
         ChartLimits[4] = z;
      }
      if(z > ChartLimits[5]){
         ChartLimits[5] = z;
      }  
      DATA.push_back({x,y,z});
   }
   
   char LimitsChoice;
   TCanvas *c = new TCanvas();
   while(true){
      TH3D *Rpoints = new TH3D("Reaction points","Reaction Points",50,ChartLimits[4],ChartLimits[5],50,ChartLimits[2],ChartLimits[3],50,ChartLimits[0],ChartLimits[1]);
       
      for(Int_t i = 0; i < DATA.size(); ++i){
         Rpoints->Fill(DATA[i][2],DATA[i][1],DATA[i][0]);
      }
      Rpoints->Draw();
      c->Modified();
      c->Update();

      cout << "Change the axes limits (XYZ/n)? ";
      cin >> LimitsChoice;
      if(LimitsChoice == 'n'){
         break;
      }
      if(LimitsChoice == 'X'){
         cout << "Type the x minimum: ";
         cin >> ChartLimits[0];
         cout << "Type the x maximum: ";
         cin >> ChartLimits[1];
      }
      if(LimitsChoice == 'Y'){
         cout << "Type the y minimum: ";
         cin >> ChartLimits[2];
         cout << "Type the y maximum: ";
         cin >> ChartLimits[3];
      }
      if(LimitsChoice == 'Z'){
         cout << "Type the z minimum: ";
         cin >> ChartLimits[4];
         cout << "Type the z maximum: ";
         cin >> ChartLimits[5];
      }
   }
   fclose(fp);
   cout << "********************************************************************" << endl;
}

void DetectionPoints(){
   cout << "********************************************************************" << endl;
   cout << "             This function plots the detection points               " << endl;
   FILE *fp  = fopen("DetectionPoints.dat","r");
   Double_t x,y,z;
   vector<double> ChartLimits = {10000,0,10000,0,10000,0};
   
   vector<vector<double>> DATA;
   char line[80];
   while(fgets(line,80,fp)) {
      sscanf(&line[0],"%lf %lf %lf", &x, &y, &z);
      if(x < ChartLimits[0]){
         ChartLimits[0] = x;
      }
      if(x > ChartLimits[1]){
         ChartLimits[1] = x;
      }   
      if(y < ChartLimits[2]){
         ChartLimits[2] = y;
      }
      if(y > ChartLimits[3]){
         ChartLimits[3] = y;
      }
      if(z < ChartLimits[4]){
         ChartLimits[4] = z;
      }
      if(z > ChartLimits[5]){
         ChartLimits[5] = z;
      }  
      DATA.push_back({x,y,z});
   }
   
   char LimitsChoice;
   TCanvas *c = new TCanvas();
   while(true){
      TH3D *Dpoints = new TH3D("Detection points","Detection Points",50,ChartLimits[4],ChartLimits[5],50,ChartLimits[2],ChartLimits[3],50,ChartLimits[0],ChartLimits[1]);
       
      for(Int_t i = 0; i < DATA.size(); ++i){
         Dpoints->Fill(DATA[i][2],DATA[i][1],DATA[i][0]);
      }
        
      Dpoints->Draw();
      c->Modified();
      c->Update();
      
      cout << "Change the axes limits (XYZ/n)? ";
      cin >> LimitsChoice;
      if(LimitsChoice == 'n'){
         break;
      }
      if(LimitsChoice == 'X'){
         cout << "Type the x minimum: ";
         cin >> ChartLimits[0];
         cout << "Type the x maximum: ";
         cin >> ChartLimits[1];
      }
      if(LimitsChoice == 'Y'){
         cout << "Type the y minimum: ";
         cin >> ChartLimits[2];
         cout << "Type the y maximum: ";
         cin >> ChartLimits[3];
      }
      if(LimitsChoice == 'Z'){
         cout << "Type the z minimum: ";
         cin >> ChartLimits[4];
         cout << "Type the z maximum: ";
         cin >> ChartLimits[5];
      }
   }
   fclose(fp);
   cout << "********************************************************************" << endl;
}

void DetectionAndReactionPoints(){
   cout << "********************************************************************" << endl;
   cout << "       This function plots the detection and reaction points        " << endl;
   FILE *Dfp  = fopen("DetectionPoints.dat","r");
   FILE *Rfp  = fopen("ReactionPoints.dat","r");
   Double_t x,y,z;
   vector<double> ChartLimits = {10000,0,10000,0,10000,0};
   vector<vector<double>> DATA;
   char line[80];
   
   while(fgets(line,80,Dfp)) {
      sscanf(&line[0],"%lf %lf %lf", &x, &y, &z);
      DATA.push_back({x,y,z});
   }
   while(fgets(line,80,Rfp)) {
      sscanf(&line[0],"%lf %lf %lf", &x, &y, &z);
      DATA.push_back({x,y,z});
   }
   for(Int_t i = 0; i < DATA.size(); ++i){
      if(DATA[i][0] < ChartLimits[0]){
         ChartLimits[0] = DATA[i][0];
      }
      if(DATA[i][0] > ChartLimits[1]){
         ChartLimits[1] = DATA[i][0];
      }   
      if(DATA[i][1] < ChartLimits[2]){
         ChartLimits[2] = DATA[i][1];
      }
      if(DATA[i][1] > ChartLimits[3]){
         ChartLimits[3] = DATA[i][1];
      }
      if(DATA[i][2] < ChartLimits[4]){
         ChartLimits[4] = DATA[i][2];
      }
      if(DATA[i][2] > ChartLimits[5]){
         ChartLimits[5] = DATA[i][2];
      }
   }
   
   
   char LimitsChoice;
   TCanvas *c = new TCanvas();
   while(true){
      TH3D *DRpoints = new TH3D("Detection and Reaction points","Detection and Reaction Points",50,ChartLimits[4],ChartLimits[5],50,ChartLimits[2],ChartLimits[3],50,ChartLimits[0],ChartLimits[1]);
       
      for(Int_t i = 0; i < DATA.size(); ++i){
         DRpoints->Fill(DATA[i][2],DATA[i][1],DATA[i][0]);
      }
      DRpoints->Draw();
      c->Modified();
      c->Update();
      
      
      cout << "Change the axes limits (XYZ/n)? ";
      cin >> LimitsChoice;
      if(LimitsChoice == 'n'){
         break;
      }
      if(LimitsChoice == 'X'){
         cout << "Type the x minimum: ";
         cin >> ChartLimits[0];
         cout << "Type the x maximum: ";
         cin >> ChartLimits[1];
      }
      if(LimitsChoice == 'Y'){
         cout << "Type the y minimum: ";
         cin >> ChartLimits[2];
         cout << "Type the y maximum: ";
         cin >> ChartLimits[3];
      }
      if(LimitsChoice == 'Z'){
         cout << "Type the z minimum: ";
         cin >> ChartLimits[4];
         cout << "Type the z maximum: ";
         cin >> ChartLimits[5];
      }
   }
   fclose(Dfp);
   fclose(Rfp);
   cout << "********************************************************************" << endl;
}
