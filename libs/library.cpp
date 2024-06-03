#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>//used to convert stream
#include <fstream>
#include <unistd.h>//used to set directory
#include <math.h> 
#include <chrono>
#include <random>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <filesystem> 

#include "lib.hpp"

#define PI 3.14159265
#define N 3
#define RESET   "\033[0m"
#define BOLD    "\033[1m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define RED     "\033[31m"
#define GREEN     "\033[32m"
using namespace std;

void get_info_from_strig(string text, vector<string> &informations){
    istringstream stream(text);
    string token;
    while (getline(stream, token, '\t')) {
        informations.push_back(token);
    }
}

void get_parameters(const string filename, vector<string> &sec_beam_list, vector<string> &contaminant_list, vector<string> &bipa_list, vector<string> &geometric_list){
    string comment_tag = "--!-->", final = "END";

    fstream input_user;
    input_user.open(filename);

    //Leitura do sec_beam
    string line;
    while(getline(input_user, line)){
        if(line.substr(0,6) == comment_tag || line.empty()){
            continue;
        }else{
            if(line == final){
                break;
            }
            int posDoisPontos = line.find(':');
            string informacao = line.substr(posDoisPontos + 1);
            informacao.erase(0, informacao.find_first_not_of(" \t\n\r\f\v"));
            sec_beam_list.push_back(informacao);
        }
    }


    //Leitura dos contaminantes
    while(getline(input_user, line)){
        if(line.substr(0,6) == comment_tag || line.empty()){
            continue;
        }
        else{
            if(line == final){
                break;
            }
            contaminant_list.push_back(line);
        }
    }

    //Leitura do biparamétrico
    while(getline(input_user, line)){
        if(line.substr(0,6) == comment_tag || line.empty()){
            continue;
        }else{
            if(line == final){
                break;
            }
            int posDoisPontos = line.find(':');
            string informacao = line.substr(posDoisPontos + 1);
            informacao.erase(0, informacao.find_first_not_of(" \t\n\r\f\v"));
            bipa_list.push_back(informacao);
        }
    }

    //Leitura da geometria do solenoide
    while(getline(input_user, line)){
        if(line.substr(0,6) == comment_tag || line.empty()){
            continue;
        }else{
            if(line == final){
                break;
            }
            int posDoisPontos = line.find(':');
            string informacao = line.substr(posDoisPontos + 1);
            informacao.erase(0, informacao.find_first_not_of(" \t\n\r\f\v"));
            geometric_list.push_back(informacao);
        }
    }
}

void Load_bar(double progress){
  int barWidth = 67;
  cout << "[";
  int pos = barWidth * progress;
  for (int k = 0; k < barWidth; ++k){
    if (k < pos) cout << "=";
    else if (k == pos) cout << ">";
    else cout << " ";
  }
  cout << "] " << fixed << setprecision(0) << progress*100 << " %\r";
  cout.flush();
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Last_simulation(char opt, string pri_beam, string sec_beam, string targ){
  string File_name = "Apps/Last_SecBeam_Simulation";
  if(opt == 'w'){
    ifstream f(File_name);
    if(f.is_open()){
      string line;
      for(int i = 0; getline(f, line); ++i){}
      cout << "********************************************************************" << RESET << endl;
      cout << YELLOW << "Warning: " << RESET << "the last Secundary beam simulation was: " << line << endl;
      cout << "********************************************************************" << RESET << endl;
    }else{
      cout << RED << "ERROR: " << RESET << "The file " << File_name << " did not existed in the current folder" << endl;
      cout << GREEN << "SOLVED: " << RESET << "The file " << File_name << " was criated and filled." << endl;
      opt = 'c';
    }
  }
  if(opt == 'c'){
    ofstream f_out(File_name, ofstream::trunc);
    if(f_out.is_open()){
      f_out << targ << "(" << pri_beam << "," << sec_beam << ")";
      f_out.close();
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////

double e_broh(double a, double q, double brho)
{
  double e;
  e=pow(brho*q,2)/(2*a);
  return e;
}

///////////////////////////////////////////////////////////////////////////////////////////////

void maxValue(vector<vector<double>> Vdata, double &max){
  max=0;
  for(long i=0; i<Vdata.size(); ++i){
    if(Vdata[i][1]>max){
      max = Vdata[i][1];
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////

double PDF_Normal(double mean, double stdev)
{
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  normal_distribution<double> distribution (mean,stdev);
  return distribution(generator);
}// construct a trivial random generator engine from a time-based seed NormalDistributionNumber:

///////////////////////////////////////////////////////////////////////////////////////////////

double PDF_RealNumber()
{
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  uniform_real_distribution<double> distribution(0.0,1.0);
  return distribution(generator);
}// construct a trivial random generator engine from a time-based seed: RealNumber

///////////////////////////////////////////////////////////////////////////////////////////////

void Particle_Mass(string particle_tag, int &particle_mass){
  if(isdigit(particle_tag[0])){
   	particle_mass = atoi(particle_tag.c_str());
  }else{
  	cout << RED << "ERROR: " << RESET << "The particle symbol was mistyped" << endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////

void interpolation(vector<vector<double>> Vdata, double xint, double &fxint ){

  double identidade[N][N];
  double a11, a12, a13, a21, a22, a23, a31, a32, a33;
  double a0, a1, a2;
  double sumY, sumXY, sumX2Y;
  double X[N], Y[N];
  
  for(long i=0; i < Vdata.size(); ++i){
    if(xint >= Vdata[i][0] && xint <= Vdata[i+2][0]){
      X[0]=Vdata[i][0];
      X[1]=Vdata[i+1][0];      
      X[2]=Vdata[i+2][0];
      Y[0]=Vdata[i][1];
      Y[1]=Vdata[i+1][1];      
      Y[2]=Vdata[i+2][1];
      i=Vdata.size();
      
    }
  }
  
  a11 = N;
  a12 = X[0] + X[1] + X[2];
  a13 = pow(X[0],2) + pow(X[1],2) + pow(X[2],2);
  a21 = a12;
  a22 = a13;
  a23 = pow(X[0],3) + pow(X[1],3) + pow(X[2],3);
  a31 = a22;
  a32 = a23;
  a33 = pow(X[0],4) + pow(X[1],4) + pow(X[2],4);

  sumY   = Y[0]+Y[1]+Y[2];
  sumXY  = X[0]*Y[0]+X[1]*Y[1]+X[2]*Y[2];
  sumX2Y = pow(X[0],2)*Y[0]+pow(X[1],2)*Y[1]+pow(X[2],2)*Y[2];
  
  double Matriz_A[N][N] = {{a11, a12, a13},{a21, a22, a23},{a31, a32, a33}};//{{n, x, x2},{x, x2, x3},{x2, x3, x4}}

  int linha, coluna, k = 0;
  double pivo = 0, m = 0;	
  for(linha = 0; linha < N; linha++){
    for(coluna = 0; coluna < N; coluna++){
      if(linha == coluna){
    	identidade[linha][coluna] = 1;
      }
      else{
        identidade[linha][coluna] = 0;     
      } 
    }    
  }

  for(coluna = 0; coluna < N; coluna++){
    pivo = Matriz_A[coluna][coluna];
    for(k = 0; k < N; k++){
      Matriz_A[coluna][k] = (Matriz_A[coluna][k])/(pivo); //L1 -> L1/pivo , L2 -> L2/pivo, L3 -> L3/pivo
      identidade[coluna][k] = (identidade[coluna][k])/(pivo); //Ex: 1 0 0 / pivo  , 0 1 0 / pivo,   0 0 1/ pivo
    }
    for(linha = 0; linha < N; linha++){
      if(linha != coluna){
      m = Matriz_A[linha][coluna];
        for(k = 0; k < N; k++){
          Matriz_A[linha][k] = (Matriz_A[linha][k]) - (m*Matriz_A[coluna][k]); //Ex: L2 -> L2 - ( 1"m" - L1) 
	  identidade[linha][k] = (identidade[linha][k]) - (m*identidade[coluna][k]);
        }
      }
    }  
  }
  a0=identidade[0][0]*sumY+identidade[1][0]*sumXY+identidade[2][0]*sumX2Y;
  a1=identidade[0][1]*sumY+identidade[1][1]*sumXY+identidade[2][1]*sumX2Y;
  a2=identidade[0][2]*sumY+identidade[1][2]*sumXY+identidade[2][2]*sumX2Y;  
  fxint = a0 + a1*xint + a2*pow(xint,2);
}

///////////////////////////////////////////////////////////////////////////////////////////////

string PointCheck(double x, double y, vector<vector<double>> data){
  //Verify if a point is under or above a determined curve
  string result;
  long posAng;
  double fx;
  interpolation(data,x,fx);
  if(y < fx){
    result = "down";
  }
  else{
    result = "above";
  }
  return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Save_Vector_Data(vector<double> V, string file_name){
  ofstream File;
  File.open(file_name,ios::app);
  for(int i = 0; i < V.size(); ++i){
      File << fixed << setprecision(10) << V[i] << "  ";
  }
  File << endl;
  File.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Primary_Target_Events(long nVector, double primaryTargThick, string CrossSection, double SigmaPrimaryBeam, double AngMin, double AngMax, vector<double> &Scatter_Angles, vector<double> &Phi_Angles, vector<vector<double>> &Reaction_Points, double &Nevents){
  //Importing the cross section data and getting the maximum value of the distribution
  vector<vector<double>> ThetaCross;
  double maxy;
  if(CrossSection == "yes" || CrossSection == "Yes" || CrossSection == "y"){
    Read_Matrix(ThetaCross,"Sec_Beam/InputFiles/crossSection.dat");
    maxValue(ThetaCross,maxy);
  }else{
    maxy = 10;
  }

  /*
  In this loop, the reaction points on the primary target 
  are obtained. For this, we consider a normal distribution 
  in the xy plane and a uniform distribution in the z axis.
  */
  cout << endl << "********************************************************************";
  cout << endl << "****    Simulanting the reaction points on the primary target   ****" << endl;
  for(long i = 0; i < nVector; ++i){
    double r, theta, z;
    vector<double> Reaction_Point;
    r = PDF_Normal(0,SigmaPrimaryBeam);
    theta = PDF_RealNumber() * 2 * PI;
    if(primaryTargThick < 0){
      z = PDF_RealNumber() * (abs(primaryTargThick))*0.01;
    }else{
      z = PDF_RealNumber() * 0.032;
    }
    Reaction_Points.push_back({r*cos(theta),r*sin(theta),z});

    Load_bar((1.0*i)/nVector);
  }

  //Printing the reaction points in the primary target
  ofstream ReactionPoints("Sec_Beam/OutputFiles/ReactionPoints.dat");
  for(long i = 0; i < Reaction_Points.size(); ++i){
    for(long j = 0; j < 3; ++j){
      ReactionPoints << Reaction_Points[i][j] << "  ";
    }
    ReactionPoints << endl;
  }

  /*
  Now, the scattering angles are obtanined from a vector
  isomorphic simulation. Besides that, the solid angle and
  the cross section are considered (this last one can be 
  disregarded). Furthermore, considering the reaction point
  and the direction of the propagation, It's verified if the
  secundary beam colids whit the first collimator or the 
  faraday cup.
  */                                                                  
  cout << endl << "\n********************************************************************";
  cout << endl << "**************        Simulating the scattering    *****************" << endl;
  vector<double> V(3,0);
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator(seed);
  uniform_real_distribution<double> VectorDistribution(-1.0,1.0);
  uniform_real_distribution<double> SigmaDistribuition(0,maxy);
  double ang, sigma;
  long j = 0;
  while(j < nVector){
    ++Nevents;
    V[0] = VectorDistribution(generator);//x
    V[1] = VectorDistribution(generator);//y
    V[2] = abs(VectorDistribution(generator));//z
    //Selecting the solid angle
    if(pow(V[0],2) + pow(V[1],2) + pow(V[2],2) <= 1){
      //Checking if scatter angle is between a deterministic range
      ang = atan(sqrt(pow(V[0],2) + pow(V[1],2))/V[2])*(180/PI);
      if(ang > AngMin && ang < AngMax){
        //Selecting point whit respect to the cross section (or dont)
        if(CrossSection == "yes" || CrossSection == "Yes" || CrossSection == "y"){
          sigma = SigmaDistribuition(generator);
          string AboveDown = PointCheck(ang,sigma,ThetaCross);
          if(AboveDown == "down"){
            ++j;
            Scatter_Angles.push_back(ang);
            if(V[1]>=0){
              Phi_Angles.push_back(acos(V[0]/sqrt(pow(V[0],2) + pow(V[1],2))));
            }
            else{
              Phi_Angles.push_back(2*PI-acos(V[0]/sqrt(pow(V[0],2) + pow(V[1],2))));
            }
          }
        }
        else{
          ++j;
          Scatter_Angles.push_back(ang);
          if(V[1]>=0){
            Phi_Angles.push_back(acos(V[0]/sqrt(pow(V[0],2) + pow(V[1],2))));
          }
          else{
            Phi_Angles.push_back(2*PI-acos(V[0]/sqrt(pow(V[0],2) + pow(V[1],2))));
          }
        }
      }
    }
    Load_bar((1.0*j)/nVector);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Reaction_Point(vector<double> &R_point, vector<double> U_vector, vector<double> I_point, double e_in_z){    
  R_point.clear();
  double t = (e_in_z)/U_vector[2];
  for(int i = 0; i < 3; ++i){
    double R_i = I_point[i] + U_vector[i] * t;
    R_point.push_back(R_i);
  }
  Save_Vector_Data(R_point,"Biparametric/Outputs/ReactionPoints.dat");
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Secundary_Target_Events(vector<double> I_point, double Theta_In, vector<double> U_vector, double e, double Detec_Radius, double Detec_Size_x, double Detec_Size_y, double Detec_Dist_Target, double Detec_Angle, double &e_in, double &e_out, double &Theta_Out, double &Scatter_Angle, double &e_out_z, vector<double> &V_vector){
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator(seed);

  //Simulating the first path on secundary target and correcting with the theta_in
  uniform_real_distribution<double> Thickness(0,e);
  double e_in_z = Thickness(generator);
  e_in = (e_in_z)/(fabs(cos(Theta_In)));


  //Calculating the reaction point (R)
  vector<double> R_point;
  Reaction_Point(R_point, U_vector, I_point, e_in_z);
  double Dx_line, Dy_line;
  //Simulating the incident point (D') on the detector
  if(Detec_Radius == 0){
    uniform_real_distribution<double> Pos_detec_x(-Detec_Size_x/200,Detec_Size_x/200);
    uniform_real_distribution<double> Pos_detec_y(-Detec_Size_y/200,Detec_Size_y/200);
    Dx_line = Pos_detec_x(generator);
    Dy_line = Pos_detec_y(generator);
  }else{
    while(true){
      uniform_real_distribution<double> Pos_detec_x(-Detec_Radius/100,Detec_Radius/100);
      uniform_real_distribution<double> Pos_detec_y(-Detec_Radius/100,Detec_Radius/100);
      double Dx_aux = Pos_detec_x(generator);
      double Dy_aux = Pos_detec_y(generator);
      if(pow(pow(Dx_aux,2) + pow(Dy_aux,2),0.5) < Detec_Radius/100){
        Dx_line = Dx_aux;
        Dy_line = Dy_aux;
        break;
      }
    }
  }

  //Converting the detection point from the detector to the secondary target referencial
  vector<double> D_detec = {Dx_line, Detec_Dist_Target * sin(Detec_Angle) - Dy_line * cos(Detec_Angle), Detec_Dist_Target * cos(Detec_Angle) + Dy_line * sin(Detec_Angle) + 2.64};
  Save_Vector_Data(D_detec,"Biparametric/Outputs/DetectionPoints.dat");

  //Exit director vector (v) and second path on the secundary target
  double V_module = sqrt(pow(D_detec[0] - R_point[0],2) + pow(D_detec[1] - R_point[1],2) + pow(D_detec[2] - R_point[2],2));
  for(int i = 0; i < 3; ++i){
      double component = (D_detec[i] - R_point[i])/V_module;
      V_vector.push_back(component);
  }
  Theta_Out = acos(V_vector[2]);
  Scatter_Angle = acos(U_vector[0]*V_vector[0] + U_vector[1]*V_vector[1] + U_vector[2]*V_vector[2]);
  if(Theta_Out < PI/2){
    e_out_z = e - e_in_z;
    e_out = e_out_z / fabs(cos(Theta_Out));
  }else{
    e_out = e_in_z / fabs(cos(Theta_Out));
  }
  
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Input_stopx(string target, string beam, double energy, double thickness)
{
  string in(".in");
  string titulo("stopx");
  string input;
  input=titulo+in;
  char *cstr = &input[0u];
  ofstream file(cstr);
  file << "proj"; file << " "; file << beam;      file <<"\r\n";
  file << "ea"; file << " ";  file << energy;     file <<"\r\n";
  file << "absb";                   	      file <<"\r\n";
  file << target; file << " "; file << fixed << setprecision(8) << thickness; file <<"\r\n";
                                                  file <<"\r\n";
  file << "elos"; 		       	      file <<"\r\n";
  file << "end"; 				      file <<"\r\n";
}//Input stopx

///////////////////////////////////////////////////////////////////////////////////////////////

void run_stopx()
{
  system("./Apps/stopx < stopx.in > stopx_out.dat");
}//run stopx

///////////////////////////////////////////////////////////////////////////////////////////////

double read_stopx()
{
  ifstream f("stopx_out.dat");
  if(!f){
    cerr << RED << "ERROR: " << RESET << "file not found: 'stopx_out.dat'" << endl;
    exit(0);
  }
  string line;
  long i;
  for (i = 0; getline(f, line); ++i){}
  f.close();
  ifstream g("stopx_out.dat");
  if(!g){
    cerr << RED << "ERROR: " << RESET << "file not found: 'stopx_out.dat'" << endl;
    exit(0);
  }
  for (long j=0; j<i-1; ++j)
  {
    getline(g, line);
  }
  g.close();
  double Energy_in,delEnergy,Energy_out;
  sscanf(line.c_str(),"%lf   %lf   %lf", &Energy_in, &delEnergy, &Energy_out);
  return Energy_out;
}// Energy_out stopx

///////////////////////////////////////////////////////////////////////////////////////////////

void Input_kineq(string target, string beam, string eject, double energy, double angle, double Q)
{
  string in(".in");
  string titulo("kineq");
  string input;
  input=titulo+in;
  char *cstr = &input[0u];
  ofstream file(cstr);
  file << target; file << "("; file << beam; file << ","; file << eject; file << ")"; file <<"\r\n";
  file <<"funa "; file << energy; file <<" "; file <<","; file<< angle; file <<","; file<< angle; file <<","; file <<"1"; file <<"\r\n";
  file <<"lev ";  file << Q;      file <<"\r\n";
  file <<"gos"; file <<"\r\n";
  file <<"end"; file <<"\r\n";
}//Input kineq

///////////////////////////////////////////////////////////////////////////////////////////////

void run_kineq()
{
  system("./Apps/kineq < kineq.in > kineq_out.dat");
}//run kineq

///////////////////////////////////////////////////////////////////////////////////////////////

double read_kineq()
{
  ifstream f("kineq_out.dat");
  if(!f){
    cerr << RED << "ERROR: " << RESET << "file not found: 'kineq_out.dat'" << endl;
    exit(0);
  }
  string line;
  long i;
  for (i = 0; getline(f, line); ++i){}
  f.close();
  ifstream g("kineq_out.dat");
  if(!g){
    cerr << RED << "ERROR: " << RESET << "file not found: 'kineq_out.dat'" << endl;
    exit(0);
  }
  for (long j=0; j<i-1; ++j)
  {
      getline(g, line);
  }
  g.close();
  double Energy_in,THLB,THCM,Energy_out;
  sscanf(line.c_str(),"%lf  %lf  %lf  %lf", &Energy_in, &THLB, &THCM, &Energy_out);
  return Energy_out;
}//Energy_out Kineq

///////////////////////////////////////////////////////////////////////////////////////////////

void Input_strag(string projectile, string target, double thickness, double energy)
{
  string in(".dat");
  string titulo("strag");
  string parameters;
  parameters = titulo+in;
  char *cstr = &parameters[0u];
  ofstream file(cstr);
  file << "proj"; file << " "; file << projectile; file << "\n";
  file << "absb"; file << "\n";
  file << target; file << " "; file << fixed << setprecision(8) << thickness; file << "\n";
  file << "\n";
  file << "ea"; file << " "; file << energy; file << "\n";
  file << "strag"; file << "\n"; 
  file << "end"; file << "\n";
}//straggling program input

///////////////////////////////////////////////////////////////////////////////////////////////

void run_strag()
{
  system("./Apps/stopx2 < strag.dat > strag_out.dat");
}//run straggling program

///////////////////////////////////////////////////////////////////////////////////////////////

double read_strag()
{
  ifstream f("strag_out.dat");
  if(!f){
    cerr << RED << "ERROR: " << RESET << "file not found: 'strag_out.dat'" << endl;
    exit(0);
  }
  string line;
  long i;
  for (i = 0; getline(f, line); ++i){}
  f.close();
  ifstream g("strag_out.dat");
  if(!g){
    cerr << RED << "ERROR: " << RESET << "file not found: 'strag_out.dat'" << endl;
    exit(0);
  }
  for (long j=0; j < i-3; ++j)
  {
    getline(g, line);
  }
  g.close();
  char txt;
  double value;
  sscanf(line.c_str(),"%s%lf",&txt, &value);
  if(value >= 0)
  {
    return value;
  }
  else
  {
    return 0;
  }
}//straggling programm output

///////////////////////////////////////////////////////////////////////////////////////////////

void cinema_input(string Proj_Tag, string Targ_Tag, string Eject_Tag, string Recoil_Tag, double Energy, double Angle, double Qreac){
  int ProjMass, TargMass, EjectilMass, RecoilMass;
  Particle_Mass(Proj_Tag, ProjMass);
  Particle_Mass(Targ_Tag, TargMass);
  Particle_Mass(Eject_Tag, EjectilMass);
  Particle_Mass(Recoil_Tag, RecoilMass);
  
  string in(".in");
  string titulo("cinema");
  string input;
  input=titulo+in;
  char *cstr = &input[0u];
  ofstream file(cstr);
  file << ProjMass;    file <<"\r\n";
  file << TargMass;    file <<"\r\n";  
  file << EjectilMass; file <<"\r\n";
  file << RecoilMass;  file <<"\r\n";
  file << Energy;      file <<"\r\n";
  file << Qreac;       file <<"\r\n";
  file << Angle;       file <<"\r\n";
  file << Angle;       file <<"\r\n";
  file << "1";         file <<"\r\n";
}

///////////////////////////////////////////////////////////////////////////////////////////////

void run_cinema(){
  system("./Apps/cinema < cinema.in > cinema.out");
}

///////////////////////////////////////////////////////////////////////////////////////////////

void read_cinema(double &En3, double &THLB4, double &En4)
{
  ifstream f("cinema.out");
  if(!f){
    cerr << RED << "ERROR: " << RESET << "file not found: 'cinema.out'" << endl;
    exit(0);
  }
  string line;
  long i;
  for (i = 0; getline(f, line); ++i){
    if(i==21)
    {
      getline(f, line);
      break;
    }
  }
  f.close();
  double THLB3, THCM3, jac;
  sscanf(line.c_str(),"%lf  %lf  %lf  %lf  %lf  %lf", &THLB3, &En3, &THCM3, &jac, &En4, &THLB4);
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Energy_Loss_stopx(double Thickness, double E_in, string Beam, string Target, string Straggling, double &E_out, double &sigE_dx1){
  Input_stopx(Target, Beam, E_in, Thickness);
  run_stopx();
  if(Straggling == "yes"){
      Input_strag(Beam, Target, Thickness, E_in);
      run_strag();
      sigE_dx1 = read_strag()/(2*sqrt(2*log(2)));
      E_out = PDF_Normal(read_stopx(), sigE_dx1);
    }//With straggling 
    else{
      E_out = read_stopx();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Energy_Loss_Reaction_Bipa(double E_in, double Theta, string Primary_Beam, string Target, string Sec_Beam, string Recoil, double Qreac, double &E_out_Sec_Beam, double &E_out_Recoil, double &Recoil_Angle){
  cinema_input(Primary_Beam, Target, Sec_Beam, Recoil, E_in, Theta, Qreac);
  run_cinema();
  read_cinema(E_out_Sec_Beam, Recoil_Angle, E_out_Recoil);
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Energy_Loss_Reaction_SB(double E_in, double Theta, string Primary_Beam, string Target, string Secundary_Beam, double Qreac, double &Energy_out){
  Input_kineq(Target, Primary_Beam, Secundary_Beam, E_in, Theta, Qreac);
  run_kineq();
  Energy_out = read_kineq();
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Current_Calculation(string secBeam, string fixcurrent, double currentFix, double massSecBeam, double Avarage_Energy, double qSecBeam, double &current){
  double Brho;
  if(secBeam == "8Li")
  {
    if(fixcurrent == "yes" || fixcurrent == "Yes" || fixcurrent == "YES")
    {
      current=currentFix;
      Brho=(current - 0.001)/49.60;
    }
    else
    {
      Brho=(0.101805)*(sqrt(2*massSecBeam*Avarage_Energy)/(qSecBeam));
      current=Brho*49.60 + 0.001;
    }
  }
  if(secBeam == "6He")
  {
    if(fixcurrent=="yes" || fixcurrent=="Yes" || fixcurrent=="YES")
    {
      current=currentFix;
      Brho=(current - 0.0346)/49.503;
    }
    else
    {
      Brho=(0.101805)*(sqrt(2*massSecBeam*Avarage_Energy)/(qSecBeam));
      current=Brho*49.503 + 0.0346;
    }
  }
  if(secBeam == "8B")
  {
    if(fixcurrent=="yes" || fixcurrent=="Yes" || fixcurrent=="YES")
    {
      current=currentFix;      
      Brho=(current - 0.329)/49.054;
    }
    else
    {
      Brho=(0.101805)*(sqrt(2*massSecBeam*Avarage_Energy)/(qSecBeam));
      current=Brho*49.054 + 0.329;
    }
  }
  if(secBeam != "8B" && secBeam != "8Li" && secBeam != "6He")
  {
    if(fixcurrent=="yes" || fixcurrent=="Yes" || fixcurrent=="YES")
    {
      current=currentFix;      
      Brho=(current + 0.024)/49.567;
    }
    else
    {
      Brho=(0.101805)*(sqrt(2*massSecBeam*Avarage_Energy)/(qSecBeam));
      current=Brho*49.567 - 0.024;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////

double avarage(vector<double> vec){
  double sum = 0;
  long nInteractions = vec.size();
  for(long i=0; i < nInteractions; ++i){
    sum = sum + vec[i];
  }  
  double Avarage = sum/nInteractions;
  return Avarage;
}//Avarage of a vector numbers

///////////////////////////////////////////////////////////////////////////////////////////////

double standard_deviation(vector<double> values){
  double med = avarage(values);
  double sum = 0;
  double diff = 0;
  double sigma_m = 0;
  long nInteractions = values.size();
  if(nInteractions > 1){
    for(long i=0; i < nInteractions; ++i){sum += pow(med - values[i],2);}
    sigma_m = sqrt((sum)/(nInteractions-1));
  }else{
    sigma_m = 0;
  }
  return sigma_m;
}

///////////////////////////////////////////////////////////////////////////////////////////////

double GetSmaller(vector<double> Values){
  double smaller = Values[0];
  for(long i=0; i < Values.size(); ++i){
    if(Values[i] < smaller){
        smaller = Values[i];
    }
  }
  return smaller;
}

///////////////////////////////////////////////////////////////////////////////////////////////

double GetLarger(vector<double> Values){
  double larger = Values[0];
  for(long i=0; i < Values.size(); ++i){
    if(Values[i] > larger){
        larger = Values[i];
    }
  }
  return larger;
}   

///////////////////////////////////////////////////////////////////////////////////////////////

double MoreRepeats(vector<double> Values, double smaller, double larger, double bin){
  double step = (larger - smaller)/bin;
  vector<double> repeats;
  while(smaller < larger){
    long n=0;
      for(long j=0; j < Values.size(); ++j){
        if(Values[j] >= smaller && Values[j] < smaller + step){
          ++n;
        }
      }
      repeats.push_back(n);
    smaller = smaller + step;
  }
  long maxy = GetLarger(repeats);
  return maxy;
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Copy_file(string file_origin, string file_target){
    string line_to_copy;
    char *in = &file_origin[0u];
    char *out = &file_target[0u];
    ifstream File_Origin(in);
    if(!File_Origin){
      cerr << RED << "ERROR: " << RESET << "file not found: " << in << endl;
      exit(EXIT_FAILURE);
    } 
    ofstream File_Target(out);
    for(int i = 0; getline(File_Origin, line_to_copy); ++ i){
      File_Target << line_to_copy << endl;
    }
    File_Origin.close();
    File_Target.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Program_Output(double Mean_Energy, double current, double Mean_energy_uncertainty, int NParticles_Target, int NParticles_Selected){
  cout << endl << "\n********************************************************************"           << endl;
  cout << endl << "Mean Energy on Target     : " << fixed << setprecision(2) << Mean_Energy << " MeV" << endl;
  cout << endl << "Standard dev mean energy  : " << fixed << setprecision(2) << Mean_energy_uncertainty << " MeV" << endl;
  cout << endl << "Solenoid Current          : " << fixed << setprecision(2) << current   << " A"    << endl;
  cout << endl << "# of particles Selected   : " << NParticles_Selected << endl;
  cout << endl << "# of particles on Target  : " << NParticles_Target << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////

double rad(double deg){
  return (PI*deg)/180;
}

///////////////////////////////////////////////////////////////////////////////////////////////

double deg(double rad){
  return (180*rad)/PI;
}

///////////////////////////////////////////////////////////////////////////////////////////////

double KinematicAngle_RecoilParticle(vector<double> V, vector<double> U, int Min, int Mout, double Ein, double Eout){
  vector<double> Pin, Precoil;
  double V_module, U_module, Precoil_module, angle, Pin_module;
  V_module = sqrt(pow(V[0],2) + pow(V[1],2) + pow(V[2],2));
  U_module = sqrt(pow(U[0],2) + pow(U[1],2) + pow(U[2],2));
  for(int i = 0; i < 3; ++i){
    double PrecoilI = sqrt(2*Min*Ein)*V[i]/V_module - sqrt(2*Mout*Eout)*U[i]/U_module;
    Precoil.push_back(PrecoilI);
  }
  Precoil_module = sqrt(pow(Precoil[0],2) + pow(Precoil[1],2) + pow(Precoil[2],2));
  angle = acos(Precoil[2]/Precoil_module);
  return angle;
}

void Read_Matrix(vector<vector<double>> &Vdata, string data){
    ifstream f(data);
    string line;
    vector<double> vec;
    double ux, uy, uz;
    if(f.is_open()){
      for(long i=0; getline(f, line); ++i){
          sscanf(line.c_str(),"%lf %lf %lf", &ux, &uy, &uz);
          vec = {ux, uy, uz};
          Vdata.push_back(vec);
          vec.clear();
      }
      f.close();
    }else{
      cout << RED << "ERROR: " << RESET << "The file " << data << " does not exist in the current folder" << endl;
      exit(0);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Read_Vector(vector<double> &vec, string data){
    ifstream f(data);
    string line;
    double component;
    if(f.is_open()){
      for(long i = 0; getline(f,line); ++i){
          sscanf(line.c_str(), "%lf", &component);
          vec.push_back(component);
      }
      f.close();
    }else{
      cout << RED << "ERROR: " << RESET << "The file " << data << " does not exist in the current folder" << endl;
      exit(0);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////

void ClearDirectorySB(){
  remove("fort.7");
  remove("fort.11");
  remove("output.dat");
  remove("rays2.dat");
  remove("rays.dat");
  remove("rays3.dat");
  remove("twsp.dat");
  remove("kineq.in");
  remove("kineq.log");
  remove("kineq_out.dat");
  remove("stopx.in");
  remove("stopx.log");
  remove("stopx_out.dat");
  remove("angEner.dat");
  remove("strag.dat");
  remove("strag_out.dat");
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Clear_directory_Bipa(){
  remove("stopx.in");
  remove("stopx.log");
  remove("stopx_out.dat");
  remove("kineq.in");
  remove("kineq.log");
  remove("kineq_out.dat");
  remove("strag.dat");
  remove("strag_out.dat");
  remove("cinema.in");
  remove("cinema.out");
  remove("CINEMA.DAT");
}

///////////////////////////////////////////////////////////////////////////////////////////////

void Join_Biparametrics(vector<string> files){
    // Output file
    ofstream outputFile("Biparametric/Outputs/Contaminants/BipaAllParticles.dat");

    // Read data from each file and write to the output file
    for (const auto& file : files) {
        ifstream inputFile(file);
        if (!inputFile) {
            cerr << RED << "ERROR: " << RESET << "opening input file: " << file << endl;
            continue;  // Skip to the next file
        }

        string line;
        while (getline(inputFile, line)) {
            istringstream iss(line);
            double col1, col2;
            if (iss >> col1 >> col2) {
                outputFile << col1 << " " << col2 << endl;
            } else {
                cerr << RED << "ERROR: " << RESET << "reading line from file: " << file << endl;
            }
        }

        inputFile.close();
    }

    outputFile.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////

void deleteDirectoryContents(const std::string& dir_path){
    for(const auto& entry : filesystem::directory_iterator(dir_path)) 
        filesystem::remove_all(entry.path());
}