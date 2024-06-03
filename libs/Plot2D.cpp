#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>//used to convert stream
#include <unistd.h>//used to set directory
#include <math.h> 
#include <chrono>
#include <random>

#include "lib.hpp"
#define RESET   "\033[0m"
#define BOLD    "\033[1m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define RED     "\033[31m"

using namespace std;
#define PI 3.14159265

///////////////////////////////////////////////////////////////////////////////////////////////

void Input_twsp(double current, double phi1, vector<double> Reaction_Point, double primaryTargThick, double qSecBeam, double massSecBeam, double EnerSecBeam, double Scatter_Angle)
{
  string in(".dat");
  string titulo("twsp");
  string input;
  input=titulo+in;
  char *cstr = &input[0u];
  ofstream file(cstr);
  file << "angEner.dat"; file <<"\r\n";
  file << phi1;    file <<"\r\n";
  file << "1.16";  file <<"\r\n";
  file << "2.80";  file <<"\r\n";
  file << Reaction_Point[0]; file << "\r\n";
  file << Reaction_Point[1]; file << "\r\n";
  file << Reaction_Point[2] - primaryTargThick/2; file << "\r\n";  
  file << "1";     file <<"\r\n";
  file << "2";     file <<"\r\n";
  file << current; file <<"\r\n"; 
  file << "n";     file <<"\r\n";
  file << "n";     file <<"\r\n";
  file << "n";     file <<"\r\n";
  file << "n";     file <<"\r\n";

  ofstream angEner;
  angEner.open("angEner.dat");
  angEner << fixed << setprecision(0) << 1 << endl;
  angEner << fixed << setprecision(4) << "  " << qSecBeam << "  " << massSecBeam<< "  " << EnerSecBeam << "  "<< Scatter_Angle << endl;
  angEner.close();
}//TWSP input

///////////////////////////////////////////////////////////////////////////////////////////////

void run_twsp()
{
  system("./Apps/twsp < twsp.dat");
}//run kineq

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void plotGraf(long nInteractions, double primaryTargThick, double Qreac, double current, vector<vector<double>> Reaction_Points, string secBeam, string primaryBeam, string primaryTarget, double alpha, double energyPrimaryBeam, double qSecBeam, double massSecBeam, double FirstCollimatorRadius, double FardayCupRadius, double LollipopRadius, double SecondCollimatorRadius, vector<double> vecEnerSecBeam, vector<double> Scatter_Angles_Temp, vector<double> Phi_Angles_Temp, vector<double> &VectorEnergySecTarget, vector<double> &Scatter_Angles, vector<double> &Phi_Angles, vector<double> &Sec_Beam_Energies, double lolli2, double Rlolli2, vector<vector<double>> &BlackPoints, vector<vector<double>> &RedPoints, vector<double> &Arrival_Angles_Target)
{
  double zplot0, xplot0, yplot0, rhoplot0, b, a, cblack, cred, cBlackTarget, cBlackFrame, cRedFrame, cRedTarget, ControladorFrame, tlolli, ttarget, A, gamma, AngTarget, LinTarget;
  double zplotCollimator1, xplotCollimator1, yplotCollimator1, rhoplotCollimator1, Collimator1Posz, Collimator1Posrho, angCollimator1, linCollimator1, Collimator1Control, tCollimator1, tTarget;
  double zplotFaradayCup, xplotFaradayCup, yplotFaradayCup, rhoplotFaradayCup, FaradayCupPosz, FaradayCupPosrho, angFaradayCup, linFaradayCup, FaradayCupControl, tFaradayCup;
  double zplotLollipop, xplotLollipop, yplotLollipop, rhoplotLollipop, LollipopPosz, LollipopPosrho, angLollipop, linLollipop, LollipopControl, tLollipop;
  double zplotCollimator2, xplotCollimator2, yplotCollimator2, rhoplotCollimator2, Collimator2Posz, Collimator2Posrho, angCollimator2, linCollimator2, Collimator2Control, TargetControl, tCollimator2;
  double auxThick;
  double xplot, yplot, zplot, xplotscatter, yplotscatter, zplotscatter;
  vector<double> xyzLeast, xyzAfter, point, VectorEnergy;
  string unitThick;
  
  if(primaryTargThick>0)
  {
    unitThick="mg/cm2";
    auxThick=primaryTargThick;
  }
  else
  { 
    unitThick="um";
    auxThick=primaryTargThick*(-10000);
  }

  ofstream raysPlot;
  ofstream rays3DPlot;
  ofstream SecEnergySecTarget;
  ofstream anglesOnTarget;
  ofstream BeamDirection;
  ofstream Points_on_Target;

  raysPlot.open("Sec_Beam/OutputCharts/raysPlot.agr");
  rays3DPlot.open("Sec_Beam/OutputFiles/rays3Dplot.dat");
  SecEnergySecTarget.open("Sec_Beam/OutputFiles/SecundaryBeamEnergyOnTarget.dat");
  anglesOnTarget.open("Sec_Beam/OutputFiles/anglesOnTarget.dat");
  BeamDirection.open("Sec_Beam/OutputFiles/BeamDirectionOnTarget.dat");
  Points_on_Target.open("Sec_Beam/OutputFiles/Points_on_Target.dat");

  raysPlot << "@version 50125"<<endl;
  raysPlot << "@page size 792, 612"<<endl;
  raysPlot << "@page scroll 5%"<<endl;
  raysPlot << "@page inout 5%"<<endl;
  raysPlot << "@link page off"<<endl;
  raysPlot << "@g0 on"<<endl;
  raysPlot << "@g0 hidden false"<<endl;
  raysPlot << "@g0 type XY"<<endl;
  raysPlot << "@with g0"<<endl;
  raysPlot << "@    world 0, 0, 3, 0.15"<<endl;
  raysPlot << "@    stack world 0, 0, 0, 0"<<endl;
  raysPlot << "@    znorm 1"<<endl;
  raysPlot << "@    view 0.150000, 0.150000, 1.150000, 0.850000"<<endl;
  raysPlot << "@    subtitle \"\"" <<endl;
  raysPlot << "@with string"<<endl;
  raysPlot << "@    string on"<<endl;
  raysPlot << "@    string loctype view"<<endl;
  raysPlot << "@    string 0.232843137255, 0.628676470588"<<endl;
  raysPlot << "@    string color 1"<<endl;
  raysPlot << "@    string rot 0"<<endl;
  raysPlot << "@    string font 60"<<endl;
  raysPlot << "@    string just 0"<<endl;
  raysPlot << "@    string char size 1.000000"<<endl;
  raysPlot << "@    string def \"Collimator 1\""<<endl;
  raysPlot << "@with string"<<endl;
  raysPlot << "@    string on"<<endl;
  raysPlot << "@    string loctype view"<<endl;
  raysPlot << "@    string 0.908088235294, 0.623774509804"<<endl;
  raysPlot << "@    string color 1"<<endl;
  raysPlot << "@    string rot 0"<<endl;
  raysPlot << "@    string font 60"<<endl;
  raysPlot << "@    string just 0"<<endl;
  raysPlot << "@    string char size 1.000000"<<endl;
  raysPlot << "@    string def \"Collimator 2\""<<endl;
  raysPlot << "@with string"<<endl;
  raysPlot << "@    string on"<<endl;
  raysPlot << "@    string loctype view"<<endl;
  raysPlot << "@    string 0.325980392157, 0.177696078431"<<endl;
  raysPlot << "@    string color 1"<<endl;
  raysPlot << "@    string rot 0"<<endl;
  raysPlot << "@    string font 60"<<endl;
  raysPlot << "@    string just 0"<<endl;
  raysPlot << "@    string char size 1.000000"<<endl;
  raysPlot << "@    string def \"Faraday cup\""<<endl;
  raysPlot << "@with string"<<endl;
  raysPlot << "@    string on"<<endl;
  raysPlot << "@    string loctype view"<<endl;
  raysPlot << "@    string 1.16544117647, 0.0526960784313"<<endl;
  raysPlot << "@    string color 2"<<endl;
  raysPlot << "@    string rot 0"<<endl;
  raysPlot << "@    string font 60"<<endl;
  raysPlot << "@    string just 0"<<endl;
  raysPlot << "@    string char size 0.750000"<<endl;
  raysPlot << "@    string def \"OCBSantos\""<<endl;
  raysPlot << "@with string"<<endl;
  raysPlot << "@    string on"<<endl;
  raysPlot << "@    string loctype view"<<endl;
  raysPlot << "@    string 1.14844117647, 0.0826960784313"<<endl;
  raysPlot << "@    string color 2"<<endl;
  raysPlot << "@    string rot 0"<<endl;
  raysPlot << "@    string font 60"<<endl;
  raysPlot << "@    string just 0"<<endl;
  raysPlot << "@    string char size 0.750000"<<endl;
  raysPlot << "@    string def \"KAlbuquerque\""<<endl;
  raysPlot << "@with string"<<endl;
  raysPlot << "@    string on"<<endl;
  raysPlot << "@    string loctype view"<<endl;
  raysPlot << "@    string 0.988970588235, 0.322303921568"<<endl;
  raysPlot << "@    string color 1"<<endl;
  raysPlot << "@    string rot 0"<<endl;
  raysPlot << "@    string font 60"<<endl;
  raysPlot << "@    string just 0"<<endl;
  raysPlot << "@    string char size 1.000000"<<endl;
  raysPlot << "@    string def \"Target\""<<endl;
  raysPlot << "@with string"<<endl;
  raysPlot << "@    string on"<<endl;
  raysPlot << "@    string loctype view"<<endl;
  raysPlot << "@    string 0.261, 0.91176"<<endl;
  raysPlot << "@    string color 1"<<endl;
  raysPlot << "@    string rot 0"<<endl;
  raysPlot << "@    string font 60"<<endl;
  raysPlot << "@    string just 0"<<endl;
  raysPlot << "@    string char size 1.000000"<<endl;
  raysPlot << "@    string def \"Pri.Beam: "<< primaryBeam <<" at Elab = "<< fixed << setprecision(2) << energyPrimaryBeam <<" MeV; Sec.Beam: "<< secBeam <<" at Elab = " << alpha << " MeV " <<"\""<<endl;
  raysPlot << "@with string"<<endl;
  raysPlot << "@    string on"<<endl;
  raysPlot << "@    string loctype view"<<endl;
  raysPlot << "@    string 0.33088, 0.875"<<endl;
  raysPlot << "@    string color 1"<<endl;
  raysPlot << "@    string rot 0"<<endl;
  raysPlot << "@    string font 60"<<endl;
  raysPlot << "@    string just 0"<<endl;
  raysPlot << "@    string char size 1.000000"<<endl;
  raysPlot << "@    string def \"Pri.Target: "<< primaryTarget <<" with thick = " << fixed << setprecision(2) << auxThick << " " << unitThick << "; Sol Current = " << current << " A\"" << endl;
  raysPlot << "@with line"<<endl;
  raysPlot << "@    line on"<<endl;
  raysPlot << "@    line loctype world"<<endl;
  raysPlot << "@    line g0"<<endl;
  raysPlot << "@    line 2.64, 0.035, 2.64, 0.025"<<endl;
  raysPlot << "@    line linewidth 2.0"<<endl;
  raysPlot << "@    line linestyle 1"<<endl;
  raysPlot << "@    line color 1"<<endl;
  raysPlot << "@    line arrow 2"<<endl;
  raysPlot << "@    line arrow type 0"<<endl;
  raysPlot << "@    line arrow length 2.000000"<<endl;
  raysPlot << "@    line arrow layout 1.000000, 1.000000"<<endl;
  raysPlot << "@line def"<<endl;
  raysPlot << "@with string"<<endl;
  raysPlot << "@    string on"<<endl;
  raysPlot << "@    string loctype view"<<endl;
  raysPlot << "@    string 0.772058823529, 0.177696078432"<<endl;
  raysPlot << "@    string color 1"<<endl;
  raysPlot << "@    string rot 0"<<endl;
  raysPlot << "@    string font 60"<<endl;
  raysPlot << "@    string just 0"<<endl;
  raysPlot << "@    string char size 1.000000"<<endl;
  raysPlot << "@    string def \"Lollipop\""<<endl;
  raysPlot << "@    yaxis  on"<<endl;
  raysPlot << "@    yaxis  tick on"<<endl;
  raysPlot << "@    yaxis  tick major 0.05"<<endl;
  raysPlot << "@    yaxis  tick minor ticks 1"<<endl;
  raysPlot << "@    yaxis  label \" \\xr\\f{} [m]\""<<endl;
  raysPlot << "@    xaxis  on"<<endl;
  raysPlot << "@    xaxis  tick on"<<endl;
  raysPlot << "@    xaxis  tick major 1.00"<<endl;
  raysPlot << "@    xaxis  tick minor ticks 4"<<endl;
  raysPlot << "@    xaxis  label \"z[m]\""<<endl;
  string line, energy;

  for(long i = 0; i < nInteractions; ++i)
  { 
    Input_twsp(current, Phi_Angles_Temp[i], Reaction_Points[i], primaryTargThick, qSecBeam, massSecBeam, vecEnerSecBeam[i], Scatter_Angles_Temp[i]);
    run_twsp();
    cblack = 0;
    cBlackTarget = 0;
    cBlackFrame = 0;
    cRedTarget = 0;
    cRedFrame = 0;
    cred = 0;
    a=0;
    b=0;
    Collimator1Control = 0;
    FaradayCupControl = 0;
    LollipopControl = 0;
    Collimator2Control = 0;
    TargetControl = 0;
    ControladorFrame = 0;
    Collimator1Posz = 0.204;
    FaradayCupPosz = 0.262;
    LollipopPosz = 2.2;
    Collimator2Posz = 2.238;
    ifstream f("rays.dat");
    if(!f){
      cerr << RED << "ERROR: " << RESET << "file not found: 'rays.dat'" << endl;
      exit(0);
    } 
    ifstream energies("Sec_Beam/OutputFiles/AllSecundaryBeamEnergies.dat");
    if(!energies){
      cerr << RED << "ERROR: " << RESET << "file not found: 'Sec_Beam/OutputFiles/AllSecundaryBeamEnergies.dat'" << endl;
      exit(0);
    } 
    VectorEnergy.clear();
    for(long j = 0; getline(energies,energy); ++j){
        double value;
        sscanf(energy.c_str(),"%lf",&value);
        VectorEnergy.push_back(value);
    }
    energies.close();
    if(Qreac==0)
    { 
      for (long j = 0; getline(f, line); ++j)
      {
        sscanf(line.c_str(),"%lf %lf %lf", &zplot, &xplot, &yplot);
        if(zplot <= Collimator1Posz){
          zplotCollimator1 = zplot;
          xplotCollimator1 = xplot;
          yplotCollimator1 = yplot;
          rhoplotCollimator1 = sqrt(pow(xplot,2)+pow(yplot,2));
        }
        if(zplot > Collimator1Posz && zplot < 0.23){
          angCollimator1=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplotCollimator1)/(zplot-zplotCollimator1);
          linCollimator1=(rhoplotCollimator1*zplot-zplotCollimator1*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplotCollimator1);
          tCollimator1 = (Collimator1Posz-zplotCollimator1)/(zplot - zplotCollimator1);
          if(angCollimator1*Collimator1Posz + linCollimator1 > FirstCollimatorRadius){
            raysPlot << fixed << setprecision(8) << Collimator1Posz <<"   "<< setprecision(12) << angCollimator1*Collimator1Posz+linCollimator1 << endl;
            rays3DPlot << fixed << setprecision(8) << Collimator1Posz << "   " << fixed << setprecision(8) << xplotCollimator1 + tCollimator1*(xplot-xplotCollimator1) << "   " << fixed << setprecision(8) << yplotCollimator1 + tCollimator1*(yplot-yplotCollimator1) << endl;
            Collimator1Control = 1;
          }
        }
        if(Collimator1Control == 0){
          if(zplot <= FaradayCupPosz){
            zplotFaradayCup = zplot;
            xplotFaradayCup = xplot;
            yplotFaradayCup = yplot;
            rhoplotFaradayCup = sqrt(pow(xplot,2)+pow(yplot,2));
          }
          if(zplot > FaradayCupPosz && zplot < 0.27){
            angFaradayCup=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplotFaradayCup)/(zplot-zplotFaradayCup);
            linFaradayCup=(rhoplotFaradayCup*zplot-zplotFaradayCup*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplotFaradayCup);
            tFaradayCup = (FaradayCupPosz-zplotFaradayCup)/(zplot - zplotFaradayCup);
            if(angFaradayCup*FaradayCupPosz + linFaradayCup < FardayCupRadius){
              raysPlot << fixed << setprecision(8) << FaradayCupPosz <<"   "<< setprecision(12) << angFaradayCup*FaradayCupPosz+linFaradayCup << endl;
              rays3DPlot << fixed << setprecision(8) << FaradayCupPosz << "   " << fixed << setprecision(8) << xplotFaradayCup + tFaradayCup*(xplot-xplotFaradayCup) << "   " << fixed << setprecision(8) << yplotFaradayCup + tFaradayCup*(yplot-yplotFaradayCup) << endl;
              FaradayCupControl = 1;
            }
          }
          if(FaradayCupControl == 0){
            if(zplot <= LollipopPosz){
              zplotLollipop = zplot;
              xplotLollipop = xplot;
              yplotLollipop = yplot;
              rhoplotLollipop = sqrt(pow(xplot,2)+pow(yplot,2));
            }
            if(zplot > LollipopPosz && zplot < 2.21){
              angLollipop=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplotLollipop)/(zplot-zplotLollipop);
              linLollipop=(rhoplotLollipop*zplot-zplotLollipop*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplotLollipop);
              tLollipop = (LollipopPosz-zplotLollipop)/(zplot - zplotLollipop);
              if(angLollipop*LollipopPosz + linLollipop < LollipopRadius){
                raysPlot << fixed << setprecision(8) << LollipopPosz <<"   "<< setprecision(12) << angLollipop*LollipopPosz+linLollipop << endl;
                rays3DPlot << fixed << setprecision(8) << LollipopPosz << "   " << fixed << setprecision(8) << xplotLollipop + tLollipop*(xplot-xplotLollipop) << "   " << fixed << setprecision(8) << yplotLollipop + tLollipop*(yplot-yplotLollipop) << endl;
                LollipopControl = 1;
              }
            }
            if(LollipopControl == 0){
              if(zplot <= Collimator2Posz){
                zplotCollimator2 = zplot;
                xplotCollimator2 = xplot;
                yplotCollimator2 = yplot;
                rhoplotCollimator2 = sqrt(pow(xplot,2)+pow(yplot,2));
              }
              if(zplot > Collimator2Posz && zplot < 2.25){
                angCollimator2=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplotCollimator2)/(zplot-zplotCollimator2);
                linCollimator2=(rhoplotCollimator2*zplot-zplotCollimator2*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplotCollimator2);
                tCollimator2 = (Collimator2Posz-zplotCollimator2)/(zplot - zplotCollimator2);
                if(angCollimator2*Collimator2Posz + linCollimator2 > SecondCollimatorRadius){
                  raysPlot << fixed << setprecision(8) << Collimator2Posz <<"   "<< setprecision(12) << angCollimator2*Collimator2Posz+linCollimator2 << endl;
                  rays3DPlot << fixed << setprecision(8) << Collimator2Posz << "   " << fixed << setprecision(8) << xplotCollimator2 + tCollimator2*(xplot-xplotCollimator2) << "   " << fixed << setprecision(8) << yplotCollimator2 + tCollimator2*(yplot-yplotCollimator2) << endl;
                  Collimator2Control = 1;
                }
              }
              if(Collimator2Control == 0){
                if(zplot <= lolli2)
                {
                  zplot0 = zplot;
                  xplot0 = xplot;
                  yplot0 = yplot;
                  rhoplot0 = sqrt(pow(xplot,2)+pow(yplot,2));
                  raysPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << sqrt(pow(xplot,2)+pow(yplot,2)) << endl;
                  rays3DPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << xplot << "   " << fixed << setprecision(8) << yplot << endl;
                }
                if(zplot > lolli2 && cblack==0)
                {
                  a=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplot0)/(zplot-zplot0);
                  b=(rhoplot0*zplot-zplot0*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplot0);
                  tlolli = (lolli2 - zplot0)/(zplot - zplot0);
                  raysPlot << fixed << setprecision(8) << lolli2 <<"   "<< setprecision(12) << a*lolli2+b << endl;
                  rays3DPlot << fixed << setprecision(8) << lolli2 << "   " << fixed << setprecision(8) << xplot0 + tlolli*(xplot-xplot0) << "   " << fixed << setprecision(8) << yplot0 + tlolli*(yplot-yplot0) << endl;
                  cblack=1;
                }
                if(cblack == 1 && a*lolli2+b > Rlolli2){
                  if(zplot <= 2.64){
                    zplot0 = zplot;
                    xplot0 = xplot;
                    yplot0 = yplot;
                    rhoplot0 = sqrt(pow(xplot,2)+pow(yplot,2));
                    raysPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << sqrt(pow(xplot,2)+pow(yplot,2)) << endl;
                    rays3DPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << xplot << "   " << fixed << setprecision(8) << yplot << endl;
                  }
                  if(zplot > 2.64 && zplot < 2.67){
                    AngTarget = (sqrt(pow(xplot,2)+pow(yplot,2))-rhoplot0)/(zplot-zplot0);
                    LinTarget = (rhoplot0*zplot-zplot0*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplot0);
                    tTarget = (2.64 - zplot0)/(zplot - zplot0);
                    if(AngTarget*2.64+LinTarget > 0.01 && AngTarget*2.64+LinTarget < 0.021 && cBlackTarget == 0){
                      raysPlot << fixed << setprecision(8) << 2.64 <<"   "<< setprecision(12) << AngTarget*2.64+LinTarget << endl;
                      rays3DPlot << fixed << setprecision(8) << 2.64 << "   " << fixed << setprecision(8) << xplot0 + tTarget*(xplot-xplot0) << "   " << fixed << setprecision(8) << yplot0 + tTarget*(yplot-yplot0) << endl; 
                      ControladorFrame = 1;
                      cBlackFrame = 1;
                    }
                    if(AngTarget*2.64+LinTarget < 0.021 && cBlackTarget == 0){
                      BlackPoints.push_back({xplot0 + tTarget*(xplot-xplot0),yplot0 + tTarget*(yplot-yplot0)});
                    }
                    if(AngTarget * 2.64 + LinTarget < 0.01 && cBlackTarget == 0){
                      raysPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << sqrt(pow(xplot,2)+pow(yplot,2)) << endl;
                      rays3DPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << xplot << "   " << fixed << setprecision(8) << yplot << endl;
                      Points_on_Target << fixed << setprecision(8) << xplot0 + tTarget*(xplot-xplot0) << "  " << yplot0 + tTarget*(yplot-yplot0) << "  " << 2.64 << endl;
                      SecEnergySecTarget << VectorEnergy[i] << endl;
                      VectorEnergySecTarget.push_back(VectorEnergy[i]);
                      A = sqrt(pow(xplot-xplot0,2) + pow(yplot-yplot0,2) + pow(zplot-zplot0,2));
                      BeamDirection << fixed << setprecision(10) << (xplot-xplot0)/A <<"  "<<(yplot-yplot0)/A <<"  "<< (zplot-zplot0)/A << endl;
                      gamma = (180/PI)*acos((zplot-zplot0)/A);
                      Arrival_Angles_Target.push_back(gamma);
                      anglesOnTarget << gamma <<endl;
                      cBlackTarget = 1;
                    }
                  }
                  if(ControladorFrame = 0 || cBlackTarget == 1 || cBlackFrame == 0){
                    raysPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << sqrt(pow(xplot,2)+pow(yplot,2)) << endl;
                    rays3DPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << xplot << "   " << fixed << setprecision(8) << yplot << endl;
                  }
                }
              }
            }
          }
        }
      }
      if(Collimator1Control == 0 && FaradayCupControl == 0){
        Scatter_Angles.push_back(Scatter_Angles_Temp[i]);
        Sec_Beam_Energies.push_back(vecEnerSecBeam[i]);
      }
      rays3DPlot << "&" << endl;
    }//only one type of trajectory

    if(Qreac>0)
    {
      if(i%2==0)
      {
        for (long j = 0; getline(f, line); ++j)
        {
          sscanf(line.c_str(),"%lf %lf %lf", &zplot, &xplot, &yplot);
          if(zplot <= Collimator1Posz){
            zplotCollimator1 = zplot;
            xplotCollimator1 = xplot;
            yplotCollimator1 = yplot;
            rhoplotCollimator1 = sqrt(pow(xplot,2)+pow(yplot,2));
          }
          if(zplot > Collimator1Posz && zplot < 0.23){
            angCollimator1=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplotCollimator1)/(zplot-zplotCollimator1);
            linCollimator1=(rhoplotCollimator1*zplot-zplotCollimator1*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplotCollimator1);
            tCollimator1 = (Collimator1Posz-zplotCollimator1)/(zplot - zplotCollimator1);
            if(angCollimator1*Collimator1Posz + linCollimator1 > 0.02125){
              raysPlot << fixed << setprecision(8) << Collimator1Posz <<"   "<< setprecision(12) << angCollimator1*Collimator1Posz+linCollimator1 << endl;
              rays3DPlot << fixed << setprecision(8) << Collimator1Posz << "   " << fixed << setprecision(8) << xplotCollimator1 + tCollimator1*(xplot-xplotCollimator1) << "   " << fixed << setprecision(8) << yplotCollimator1 + tCollimator1*(yplot-yplotCollimator1) << endl;
              Collimator1Control = 1;
            }
          }
          if(Collimator1Control == 0){
            if(zplot <= FaradayCupPosz){
              zplotFaradayCup = zplot;
              xplotFaradayCup = xplot;
              yplotFaradayCup = yplot;
              rhoplotFaradayCup = sqrt(pow(xplot,2)+pow(yplot,2));
            }
            if(zplot > FaradayCupPosz && zplot < 0.27){
              angFaradayCup=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplotFaradayCup)/(zplot-zplotFaradayCup);
              linFaradayCup=(rhoplotFaradayCup*zplot-zplotFaradayCup*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplotFaradayCup);
              tFaradayCup = (FaradayCupPosz-zplotFaradayCup)/(zplot - zplotFaradayCup);
              if(angFaradayCup*FaradayCupPosz + linFaradayCup < 0.013625){
                raysPlot << fixed << setprecision(8) << FaradayCupPosz <<"   "<< setprecision(12) << angFaradayCup*FaradayCupPosz+linFaradayCup << endl;
                rays3DPlot << fixed << setprecision(8) << FaradayCupPosz << "   " << fixed << setprecision(8) << xplotFaradayCup + tFaradayCup*(xplot-xplotFaradayCup) << "   " << fixed << setprecision(8) << yplotFaradayCup + tFaradayCup*(yplot-yplotFaradayCup) << endl;
                FaradayCupControl = 1;
              }
            }
            if(FaradayCupControl == 0){
              if(zplot <= LollipopPosz){
                zplotLollipop = zplot;
                xplotLollipop = xplot;
                yplotLollipop = yplot;
                rhoplotLollipop = sqrt(pow(xplot,2)+pow(yplot,2));
              }
              if(zplot > LollipopPosz && zplot < 2.21){
                angLollipop=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplotLollipop)/(zplot-zplotLollipop);
                linLollipop=(rhoplotLollipop*zplot-zplotLollipop*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplotLollipop);
                tLollipop = (LollipopPosz-zplotLollipop)/(zplot - zplotLollipop);
                if(angLollipop*LollipopPosz + linLollipop < LollipopRadius){
                  raysPlot << fixed << setprecision(8) << LollipopPosz <<"   "<< setprecision(12) << angLollipop*LollipopPosz+linLollipop << endl;
                  rays3DPlot << fixed << setprecision(8) << LollipopPosz << "   " << fixed << setprecision(8) << xplotLollipop + tLollipop*(xplot-xplotLollipop) << "   " << fixed << setprecision(8) << yplotLollipop + tLollipop*(yplot-yplotLollipop) << endl;
                  LollipopControl = 1;
                }
              }
              if(LollipopControl == 0){
                if(zplot <= Collimator2Posz){
                  zplotCollimator2 = zplot;
                  xplotCollimator2 = xplot;
                  yplotCollimator2 = yplot;
                  rhoplotCollimator2 = sqrt(pow(xplot,2)+pow(yplot,2));
                }
                if(zplot > Collimator2Posz && zplot < 2.25){
                  angCollimator2=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplotCollimator2)/(zplot-zplotCollimator2);
                  linCollimator2=(rhoplotCollimator2*zplot-zplotCollimator2*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplotCollimator2);
                  tCollimator2 = (Collimator2Posz-zplotCollimator2)/(zplot - zplotCollimator2);
                  if(angCollimator2*Collimator2Posz + linCollimator2 > SecondCollimatorRadius){
                    raysPlot << fixed << setprecision(8) << Collimator2Posz <<"   "<< setprecision(12) << angCollimator2*Collimator2Posz+linCollimator2 << endl;
                    rays3DPlot << fixed << setprecision(8) << Collimator2Posz << "   " << fixed << setprecision(8) << xplotCollimator2 + tCollimator2*(xplot-xplotCollimator2) << "   " << fixed << setprecision(8) << yplotCollimator2 + tCollimator2*(yplot-yplotCollimator2) << endl;
                    Collimator2Control = 1;
                  }
                }
                if(Collimator2Control == 0){
                  if(zplot <= lolli2)
                  {
                    zplot0 = zplot;
                    xplot0 = xplot;
                    yplot0 = yplot;
                    rhoplot0 = sqrt(pow(xplot,2)+pow(yplot,2));
                    raysPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << sqrt(pow(xplot,2)+pow(yplot,2)) << endl;
                    rays3DPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << xplot << "   " << fixed << setprecision(8) << yplot << endl;
                  }
                  if(zplot > lolli2 && cblack==0)
                  {
                    a=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplot0)/(zplot-zplot0);
                    b=(rhoplot0*zplot-zplot0*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplot0);
                    tlolli = (lolli2 - zplot0)/(zplot - zplot0);
                    raysPlot << fixed << setprecision(8) << lolli2 <<"   "<< setprecision(12) << a*lolli2+b << endl;
                    rays3DPlot << fixed << setprecision(8) << lolli2 << "   " << fixed << setprecision(8) << xplot0 + tlolli*(xplot-xplot0) << "   " << fixed << setprecision(8) << yplot0 + tlolli*(yplot-yplot0) << endl;
                    cblack=1;
                  }
                  if(cblack == 1 && a*lolli2+b > Rlolli2)
                  {
                    if(zplot <= 2.64){
                      zplot0 = zplot;
                      xplot0 = xplot;
                      yplot0 = yplot;
                      rhoplot0 = sqrt(pow(xplot,2)+pow(yplot,2));
                      raysPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << sqrt(pow(xplot,2)+pow(yplot,2)) << endl;
                      rays3DPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << xplot << "   " << fixed << setprecision(8) << yplot << endl;
                    }
                    if(zplot > 2.64 && zplot < 2.67){
                      AngTarget = (sqrt(pow(xplot,2)+pow(yplot,2))-rhoplot0)/(zplot-zplot0);
                      LinTarget = (rhoplot0*zplot-zplot0*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplot0);
                      tTarget = (2.64 - zplot0)/(zplot - zplot0);
                      if(AngTarget*2.64+LinTarget > 0.01 && AngTarget*2.64+LinTarget < 0.021 && cBlackTarget == 0){
                        raysPlot << fixed << setprecision(8) << 2.64 <<"   "<< setprecision(12) << AngTarget*2.64+LinTarget << endl;
                        rays3DPlot << fixed << setprecision(8) << 2.64 << "   " << fixed << setprecision(8) << xplot0 + tTarget*(xplot-xplot0) << "   " << fixed << setprecision(8) << yplot0 + tTarget*(yplot-yplot0) << endl; 
                        ControladorFrame = 1;
                        cBlackFrame = 1;
                      }
                      if(AngTarget*2.64+LinTarget < 0.021 && cBlackTarget == 0){
                        BlackPoints.push_back({xplot0 + tTarget*(xplot-xplot0),yplot0 + tTarget*(yplot-yplot0)});
                      }
                      if(AngTarget * 2.64 + LinTarget < 0.01 && cBlackTarget == 0){
                        raysPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << sqrt(pow(xplot,2)+pow(yplot,2)) << endl;
                        rays3DPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << xplot << "   " << fixed << setprecision(8) << yplot << endl;
                        Points_on_Target << fixed << setprecision(8) << xplot0 + tTarget*(xplot-xplot0) << "  " << yplot0 + tTarget*(yplot-yplot0) << "  " << 2.64 << endl;
                        SecEnergySecTarget << VectorEnergy[i] << endl;
                        VectorEnergySecTarget.push_back(VectorEnergy[i]);
                        A = sqrt(pow(xplot-xplot0,2) + pow(yplot-yplot0,2) + pow(zplot-zplot0,2));
                        BeamDirection << fixed << setprecision(10) << (xplot-xplot0)/A <<"  "<<(yplot-yplot0)/A <<"  "<< (zplot-zplot0)/A << endl;
                        gamma = (180/PI)*acos((zplot-zplot0)/A);
                        Arrival_Angles_Target.push_back(gamma);
                        anglesOnTarget << gamma <<endl;
                        cBlackTarget = 1;
                      }
                    }
                    if(ControladorFrame = 0 || cBlackTarget == 1 || cBlackFrame == 0){
                      raysPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << sqrt(pow(xplot,2)+pow(yplot,2)) << endl;
                      rays3DPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << xplot << "   " << fixed << setprecision(8) << yplot << endl;
                    }
                  }
                }
              }
            }
          }
        }
        if(Collimator1Control == 0 && FaradayCupControl == 0){
          Scatter_Angles.push_back(Scatter_Angles_Temp[i]);
          Sec_Beam_Energies.push_back(vecEnerSecBeam[i]);
        }
        rays3DPlot << "&" << endl;
      }  
      else{
        for (long j = 0; getline(f, line); ++j){
          sscanf(line.c_str(),"%lf %lf %lf", &zplot, &xplot, &yplot);
          if(zplot <= Collimator1Posz){
            zplotCollimator1 = zplot;
            xplotCollimator1 = xplot;
            yplotCollimator1 = yplot;
            rhoplotCollimator1 = sqrt(pow(xplot,2)+pow(yplot,2));
          }
          if(zplot > Collimator1Posz && zplot < 0.23){
            angCollimator1=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplotCollimator1)/(zplot-zplotCollimator1);
            linCollimator1=(rhoplotCollimator1*zplot-zplotCollimator1*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplotCollimator1);
            tCollimator1 = (Collimator1Posz-zplotCollimator1)/(zplot - zplotCollimator1);
            if(angCollimator1*Collimator1Posz + linCollimator1 > 0.02125){
              raysPlot << fixed << setprecision(8) << Collimator1Posz <<"   "<< setprecision(12) << angCollimator1*Collimator1Posz+linCollimator1 << endl;
              rays3DPlot << fixed << setprecision(8) << Collimator1Posz << "   " << fixed << setprecision(8) << xplotCollimator1 + tCollimator1*(xplot-xplotCollimator1) << "   " << fixed << setprecision(8) << yplotCollimator1 + tCollimator1*(yplot-yplotCollimator1) << endl;
              Collimator1Control = 1;
            }
          }
          if(Collimator1Control == 0){
            if(zplot <= FaradayCupPosz){
              zplotFaradayCup = zplot;
              xplotFaradayCup = xplot;
              yplotFaradayCup = yplot;
              rhoplotFaradayCup = sqrt(pow(xplot,2)+pow(yplot,2));
            }
            if(zplot > FaradayCupPosz && zplot < 0.27){
              angFaradayCup=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplotFaradayCup)/(zplot-zplotFaradayCup);
              linFaradayCup=(rhoplotFaradayCup*zplot-zplotFaradayCup*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplotFaradayCup);
              tFaradayCup = (FaradayCupPosz-zplotFaradayCup)/(zplot - zplotFaradayCup);
              if(angFaradayCup*FaradayCupPosz + linFaradayCup < 0.013625){
                raysPlot << fixed << setprecision(8) << FaradayCupPosz <<"   "<< setprecision(12) << angFaradayCup*FaradayCupPosz+linFaradayCup << endl;
                rays3DPlot << fixed << setprecision(8) << FaradayCupPosz << "   " << fixed << setprecision(8) << xplotFaradayCup + tFaradayCup*(xplot-xplotFaradayCup) << "   " << fixed << setprecision(8) << yplotFaradayCup + tFaradayCup*(yplot-yplotFaradayCup) << endl;
                FaradayCupControl = 1;
              }
            }
            if(FaradayCupControl == 0){
              if(zplot <= LollipopPosz){
                zplotLollipop = zplot;
                xplotLollipop = xplot;
                yplotLollipop = yplot;
                rhoplotLollipop = sqrt(pow(xplot,2)+pow(yplot,2));
              }
              if(zplot > LollipopPosz && zplot < 2.21){
                angLollipop=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplotLollipop)/(zplot-zplotLollipop);
                linLollipop=(rhoplotLollipop*zplot-zplotLollipop*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplotLollipop);
                tLollipop = (LollipopPosz-zplotLollipop)/(zplot - zplotLollipop);
                if(angLollipop*LollipopPosz + linLollipop < LollipopRadius){
                  raysPlot << fixed << setprecision(8) << LollipopPosz <<"   "<< setprecision(12) << angLollipop*LollipopPosz+linLollipop << endl;
                  rays3DPlot << fixed << setprecision(8) << LollipopPosz << "   " << fixed << setprecision(8) << xplotLollipop + tLollipop*(xplot-xplotLollipop) << "   " << fixed << setprecision(8) << yplotLollipop + tLollipop*(yplot-yplotLollipop) << endl;
                  LollipopControl = 1;
                }
              }
              if(LollipopControl == 0){
                if(zplot <= Collimator2Posz){
                  zplotCollimator2 = zplot;
                  xplotCollimator2 = xplot;
                  yplotCollimator2 = yplot;
                  rhoplotCollimator2 = sqrt(pow(xplot,2)+pow(yplot,2));
                }
                if(zplot > Collimator2Posz && zplot < 2.25){
                  angCollimator2=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplotCollimator2)/(zplot-zplotCollimator2);
                  linCollimator2=(rhoplotCollimator2*zplot-zplotCollimator2*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplotCollimator2);
                  tCollimator2 = (Collimator2Posz-zplotCollimator2)/(zplot - zplotCollimator2);
                  if(angCollimator2*Collimator2Posz + linCollimator2 > SecondCollimatorRadius){
                    raysPlot << fixed << setprecision(8) << Collimator2Posz <<"   "<< setprecision(12) << angCollimator2*Collimator2Posz+linCollimator2 << endl;
                    rays3DPlot << fixed << setprecision(8) << Collimator2Posz << "   " << fixed << setprecision(8) << xplotCollimator2 + tCollimator2*(xplot-xplotCollimator2) << "   " << fixed << setprecision(8) << yplotCollimator2 + tCollimator2*(yplot-yplotCollimator2) << endl;
                    Collimator2Control = 1;
                  }
                }
                if(Collimator2Control == 0){
                  if(zplot <= lolli2)
                  {
                    zplot0 = zplot;
                    xplot0 = xplot;
                    yplot0 = yplot;
                    rhoplot0 = sqrt(pow(xplot,2)+pow(yplot,2));
                    raysPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << sqrt(pow(xplot,2)+pow(yplot,2)) << endl;
                    rays3DPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << xplot << "   " << fixed << setprecision(8) << yplot << endl;
                  }
                  if(zplot > lolli2 && cred==0)
                  {
                    a=(sqrt(pow(xplot,2)+pow(yplot,2))-rhoplot0)/(zplot-zplot0);
                    b=(rhoplot0*zplot-zplot0*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplot0);
                    tlolli = (lolli2 - zplot0)/(zplot - zplot0);
                    raysPlot << fixed << setprecision(8) << lolli2 <<"   "<< setprecision(12) << a*lolli2+b << endl;
                    rays3DPlot << fixed << setprecision(8) << lolli2 << "   " << fixed << setprecision(8) << xplot0 + tlolli*(xplot-xplot0) << "   " << fixed << setprecision(8) << yplot0 + tlolli*(yplot-yplot0) << endl;
                    cred=1;
                  }
                  if(cred == 1 && a*lolli2+b > Rlolli2)
                  {
                    if(zplot <= 2.64){
                      zplot0 = zplot;
                      xplot0 = xplot;
                      yplot0 = yplot;
                      rhoplot0 = sqrt(pow(xplot,2)+pow(yplot,2));
                      raysPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << sqrt(pow(xplot,2)+pow(yplot,2)) << endl;
                      rays3DPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << xplot << "   " << fixed << setprecision(8) << yplot << endl;
                    }
                    if(zplot > 2.64 && zplot < 2.67){
                      AngTarget = (sqrt(pow(xplot,2)+pow(yplot,2))-rhoplot0)/(zplot-zplot0);
                      LinTarget = (rhoplot0*zplot-zplot0*sqrt(pow(xplot,2)+pow(yplot,2)))/(zplot-zplot0);
                      tTarget = (2.64 - zplot0)/(zplot - zplot0);
                      if(AngTarget*2.64+LinTarget > 0.01 && AngTarget*2.64+LinTarget < 0.021 && cRedTarget == 0){
                        raysPlot << fixed << setprecision(8) << 2.64 <<"   "<< setprecision(12) << AngTarget*2.64+LinTarget << endl;
                        rays3DPlot << fixed << setprecision(8) << 2.64 << "   " << fixed << setprecision(8) << xplot0 + tTarget*(xplot-xplot0) << "   " << fixed << setprecision(8) << yplot0 + tTarget*(yplot-yplot0) << endl; 
                        ControladorFrame = 1;
                        cRedFrame = 1;
                      }
                      if(AngTarget*2.64+LinTarget < 0.021 && cRedTarget == 0){
                        RedPoints.push_back({xplot0 + tTarget*(xplot-xplot0),yplot0 + tTarget*(yplot-yplot0)});
                      }
                      if(AngTarget * 2.64 + LinTarget < 0.01 && cRedTarget == 0){
                        raysPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << sqrt(pow(xplot,2)+pow(yplot,2)) << endl;
                        rays3DPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << xplot << "   " << fixed << setprecision(8) << yplot << endl;
                        Points_on_Target << fixed << setprecision(8) << xplot0 + tTarget*(xplot-xplot0) << "  " << yplot0 + tTarget*(yplot-yplot0) << "  " << 2.64 << endl;
                        SecEnergySecTarget << VectorEnergy[i] << endl;
                        VectorEnergySecTarget.push_back(VectorEnergy[i]);
                        A = sqrt(pow(xplot-xplot0,2) + pow(yplot-yplot0,2) + pow(zplot-zplot0,2));
                        BeamDirection << fixed << setprecision(10) << (xplot-xplot0)/A <<"  "<<(yplot-yplot0)/A <<"  "<< (zplot-zplot0)/A << endl;
                        gamma = (180/PI)*acos((zplot-zplot0)/A);
                        Arrival_Angles_Target.push_back(gamma);
                        anglesOnTarget << gamma <<endl;
                        cRedTarget = 1;
                      }
                    }
                    if(ControladorFrame = 0 || cRedTarget == 1 || cRedFrame == 0){
                      raysPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << sqrt(pow(xplot,2)+pow(yplot,2)) << endl;
                      rays3DPlot << fixed << setprecision(8) << zplot << "   " << fixed << setprecision(8) << xplot << "   " << fixed << setprecision(8) << yplot << endl;
                    }
                  }
                }
              }
            }
          }
        }
        if(Collimator1Control == 0 && FaradayCupControl == 0){
          Scatter_Angles.push_back(Scatter_Angles_Temp[i]);
          Sec_Beam_Energies.push_back(vecEnerSecBeam[i]);
        }
        rays3DPlot << "&" << endl;
      }
    }// two types of trajectories
  
    
    if(Qreac>0)
    {
      if(i%2==0)
      {
        raysPlot << "@    s"<<i<<" line color "<<1 << endl;
      }
      else
      {
        raysPlot << "@    s"<<i<<" line color "<<2 << endl;      
      }
    }
    else
    {
      raysPlot << "@    s"<< i <<" line color "<< 1 << endl;      
    }  
    raysPlot << "&" << endl;
    f.close();

    Load_bar((1.0*i)/nInteractions);
  }
  raysPlot << "2.2  0.00"<< endl;
  raysPlot << "2.2  \"" << LollipopRadius << "\"" << endl;
  raysPlot << "@    s"<<nInteractions <<" line linewidth 1.5"<<endl;
  raysPlot << "@    s"<<nInteractions <<" line color 3"<< endl;
  raysPlot << "&"<< endl;
  raysPlot << "0.262  0.000"<< endl;
  raysPlot << "0.262  0.013625"<< endl;
  raysPlot << "0.442  0.013625"<< endl;
  raysPlot << "0.442  0.0000"<< endl;
  raysPlot << "0.262  0.000"<< endl;
  raysPlot << "@    s"<<nInteractions+1<<" line linewidth 1.5"<<endl;
  raysPlot << "@    s"<<nInteractions+1<<" line color 3"<< endl;
  raysPlot << "@    s"<<nInteractions+1<<" fill type 1"<< endl;
  raysPlot << "@    s"<<nInteractions+1<<" fill rule 0"<< endl;
  raysPlot << "@    s"<<nInteractions+1<<" fill color 3"<< endl;
  raysPlot << "@    s"<<nInteractions+1<<" fill pattern 1"<< endl;
  raysPlot << "&"<< endl;
  raysPlot << "0.204	0.02125"<< endl;
  raysPlot << "0.204	0.15000"<< endl;
  raysPlot << "@    s"<<nInteractions+2<<" line linewidth 1.5"<<endl;
  raysPlot << "@    s"<<nInteractions+2<<" line color 3"<< endl;
  raysPlot << "&"<< endl;
  raysPlot << "2.238	0.03000"<< endl;
  raysPlot << "2.238	0.15000"<< endl;
  raysPlot << "@    s"<<nInteractions+3<<" line linewidth 1.5"<<endl;
  raysPlot << "@    s"<<nInteractions+3<<" line color 3"<< endl;
  raysPlot << "&"<< endl;
  raysPlot << "2.64	0.00"<< endl;
  raysPlot << "2.64	0.01"<< endl;
  raysPlot << "@    s"<<nInteractions+4<<" line linewidth 1.5"<<endl;
  raysPlot << "@    s"<<nInteractions+4<<" line color 3"<< endl;
  raysPlot << "&"<< endl;
  raysPlot << "2.64	0.01"<< endl;
  raysPlot << "2.64	0.021"<< endl;
  raysPlot << "@    s"<<nInteractions+5<<" line linewidth 2"<<endl;
  raysPlot << "@    s"<<nInteractions+5<<" line color 4"<< endl;
  raysPlot << "&"<< endl;
  raysPlot << lolli2 << "  0.000 "<< endl;
  raysPlot << lolli2 << "  "      << Rlolli2 << endl;//em Metros 
  raysPlot << "@    s"<<nInteractions+6<<" line linewidth 1.5"<<endl;
  raysPlot << "@    s"<<nInteractions+6<<" line color 3"<< endl;
  raysPlot.close();
  anglesOnTarget.close();
  BeamDirection.close();
}  

///////////////////////////////////////////////////////////////////////////////////////////////
