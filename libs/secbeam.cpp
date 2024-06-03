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
#include <cstdlib>
#include "lib.hpp"
#include <stdio.h>
#include <string.h>
#define PI 3.14159265
#define RESET   "\033[0m"
#define BOLD    "\033[1m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define RED     "\033[31m"

using namespace std;

void secbeam(long nVector, double primaryTargThick, double massSecBeam, double qSecBeam, double Primary_Beam_Avarage_Energy, double Qreac, double ScatterAngMin, double ScatterAngMax, double FirstCollimatorRadius, double FaradayCupRadius, double LollipopRadius, double SecondCollimatorRadius, double SigmaPrimaryBeam, double SigmaPrimaryEnergy, double lolli2, double Rloli2, double currentFix, string primaryBeam, string primaryTarget, string secBeam, string CrossSection, string fixcurrent, string straggling, vector<string> Contaminant_Tag, vector<double> Contaminant_Mass, vector<double> Contaminant_Charge){
    cout << endl << "********************************************************************" << endl;
    cout <<         "***************    Secundary beam in RIBRAS system   ***************" << endl;
    //############################################################################################
    //Simulating all events in the primary target as well as the scattering angle
    double Nevents;
    vector<double> Scatter_Angles_Temp, Phi_Angles_Temp;
    vector<vector<double>> Reaction_Points;

    Primary_Target_Events(nVector, primaryTargThick, CrossSection, SigmaPrimaryBeam, ScatterAngMin, ScatterAngMax, Scatter_Angles_Temp, Phi_Angles_Temp, Reaction_Points, Nevents);

    //############################################################################################
    cout << endl << "\n********************************************************************" << endl;
    cout <<         "**********    Calculating the secundary beam energies    ***********" << endl;
    //############################################################################################
    

    //############################################################################################   
    //Initializing energy variables and reading straggling option
    double Energy_out_dx1, Energy_out_kineq, Energy_out_dx2, Energy_out_havar1, Energy_out_havar2, dx1PrimaryTarg, dx2PrimaryTarg;
    double sigE_dx1, sigE_dx2, sigEhavar; 
    vector<double> strag1, strag2, stragT, allEnerSecBeam;
    ofstream Primary_Beam_Energies("Sec_Beam/OutputFiles/PrimaryBeamEnergies.dat");
    ofstream SecundaryBeamEnergies("Sec_Beam/OutputFiles/AllSecundaryBeamEnergies.dat");
    ofstream ScatterAnglesSecBeam("Sec_Beam/OutputFiles/AllScatterAnglesSecBeam.dat");

    //Calculates all energies loses in the target and reaction
    for(int i = 0; i < nVector; ++i){
        double Primary_Beam_Energy = PDF_Normal(Primary_Beam_Avarage_Energy, SigmaPrimaryEnergy);
        Primary_Beam_Energies << Primary_Beam_Energy << endl;
        double Q = (Qreac > 0 && i % 2 != 0) ? Qreac : 0;

        if(primaryTargThick < 0){//For solid target
            dx1PrimaryTarg = Reaction_Points[i][2]*(-100);
            dx2PrimaryTarg = -(fabs(primaryTargThick)-fabs(dx1PrimaryTarg))/cos((PI/180)*Scatter_Angles_Temp[i]);
            //Compute loss of energy in the first way before reaction
            Energy_Loss_stopx(dx1PrimaryTarg, Primary_Beam_Energy, primaryBeam, primaryTarget, straggling, Energy_out_dx1, sigE_dx1);
            strag1.push_back((2*sqrt(2*log(2)))*sigE_dx1);
            //Compute loss of energy in the reaction
            Energy_Loss_Reaction_SB(Energy_out_dx1, Scatter_Angles_Temp[i], primaryBeam, primaryTarget, secBeam, Q, Energy_out_kineq);
            //Compute loss of energy in the second way after reaction
            Energy_Loss_stopx(dx2PrimaryTarg, Energy_out_kineq, secBeam, primaryTarget, straggling, Energy_out_dx2, sigE_dx2);
            allEnerSecBeam.push_back(Energy_out_dx2);
            SecundaryBeamEnergies << Energy_out_dx2 << endl;
        }
        if(primaryTargThick > 0){//For gas target
            double In_Factor = (Reaction_Points[i][2]*(100))/3.2;
            double Out_Factor = (3.2 - Reaction_Points[i][2]*(100))/3.2;
            dx1PrimaryTarg = In_Factor * primaryTargThick;
            dx2PrimaryTarg = (Out_Factor * primaryTargThick)/cos((PI/180)*Scatter_Angles_Temp[i]);
            //Compute loss of energy in the havar leaves before entering the gas
            Energy_Loss_stopx(-0.0002, Primary_Beam_Energy, primaryBeam, "58Ni", straggling, Energy_out_havar1, sigEhavar);
            //Compute loss of energy in the first way before reaction
            Energy_Loss_stopx(dx1PrimaryTarg, Energy_out_havar1, primaryBeam, primaryTarget, straggling, Energy_out_dx1, sigE_dx1);
            strag1.push_back((2*sqrt(2*log(2)))*sigE_dx1);
            //Compute loss of energy in the reaction
            Energy_Loss_Reaction_SB(Energy_out_dx1, Scatter_Angles_Temp[i], primaryBeam, primaryTarget, secBeam, Q, Energy_out_kineq);
            //Compute loss of energy in the second way after reaction
            Energy_Loss_stopx(dx2PrimaryTarg, Energy_out_kineq, secBeam, primaryTarget, straggling, Energy_out_dx2, sigE_dx2);
            strag2.push_back((2*sqrt(2*log(2)))*sigE_dx2);
            //Compute loss of energy in the havar leaves before entering the gas
            Energy_Loss_stopx(-0.0002, Energy_out_dx2, secBeam, "58Ni", straggling, Energy_out_havar2, sigEhavar);
            allEnerSecBeam.push_back(Energy_out_havar2);
            SecundaryBeamEnergies << Energy_out_havar2 << endl;
        }

        ScatterAnglesSecBeam << Scatter_Angles_Temp[i] << endl;
        Load_bar((1.0*i)/nVector);
    }
    double Avarage_Energy = avarage(allEnerSecBeam);
    Primary_Beam_Energies.close();
    SecundaryBeamEnergies.close();
    ScatterAnglesSecBeam.close();

    //############################################################################################
    //Calculates the magnetic ridity in relation to the average energy (Avarage energy); Calculates the current and the adjusted current
    double current;
    Current_Calculation(secBeam, fixcurrent, currentFix, massSecBeam, Avarage_Energy, qSecBeam, current);

    //############################################################################################
    cout << endl << "\n********************************************************************";
    cout << endl << "****************    Calculating the trajectories    ****************" << endl;
    //############################################################################################   
    //Plot the (z , rho) trajectories whit grace
    vector<double> Arrival_Angles_Target, Scatter_Angles, Phi_Angles, Secundary_Beam_Energies, Secundary_Beam_Energies_Target;
    vector<vector<double>> BlackPoints, RedPoints;
    plotGraf(nVector, primaryTargThick, Qreac, current, Reaction_Points, secBeam, primaryBeam, primaryTarget, Avarage_Energy, Primary_Beam_Avarage_Energy, qSecBeam, massSecBeam, FirstCollimatorRadius, FaradayCupRadius, LollipopRadius, SecondCollimatorRadius, allEnerSecBeam, Scatter_Angles_Temp, Phi_Angles_Temp, Secundary_Beam_Energies_Target, Scatter_Angles, Phi_Angles, Secundary_Beam_Energies, lolli2, Rloli2, BlackPoints, RedPoints, Arrival_Angles_Target);
    ofstream ScatterAnglesSelected("Sec_Beam/OutputFiles/ScatterAnglesSelected.dat");
    ofstream SecundaryBeamEnergiesSelected("Sec_Beam/OutputFiles/SecundaryBeamEnergiesSelected.dat");
    for(long i = 0; i < Scatter_Angles.size(); ++ i){
        ScatterAnglesSelected << Scatter_Angles[i] << endl;
        SecundaryBeamEnergiesSelected << Secundary_Beam_Energies[i] << endl;
    }
    ScatterAnglesSelected.close();
    SecundaryBeamEnergiesSelected.close();
    double Mean_Energy = avarage(Secundary_Beam_Energies_Target);
    double Mean_Energy_dispersion = standard_deviation(Secundary_Beam_Energies_Target);
    //############################################################################################   
    //Calculating the energies of the contaminants
    deleteDirectoryContents("Biparametric/Inputs");
    deleteDirectoryContents("Sec_Beam/OutputFiles/Contaminants");

    for (int i = 0; i < Contaminant_Tag.size(); ++i) {
        string tag = Contaminant_Tag[i];
        double charge = Contaminant_Charge[i];
        stringstream formatted_charge;
        formatted_charge << fixed << setprecision(0) << charge;
        string file_name = "Sec_Beam/OutputFiles/Contaminants/" + tag + "_Q" + formatted_charge.str() + "_Energy.dat";

        ofstream ContaminantsEnergies(file_name);
        for(int j = 0; j < Secundary_Beam_Energies_Target.size(); ++j){
            double Energy_Contaminant = (massSecBeam/Contaminant_Mass[i]) * pow((Contaminant_Charge[i]/qSecBeam),2) * Secundary_Beam_Energies_Target[j];
            //scout << fixed << setprecision(8) << Contaminant_Tag[i] << "  " << Contaminant_Mass[i] << "  "  << Contaminant_Charge[i] << "  " << Energy_Contaminant << "  " << Secundary_Beam_Energies_Target[j] << endl;
            ContaminantsEnergies << Energy_Contaminant << endl;
        }

        string file_name_bipa = "Biparametric/Inputs/" + tag + "_Q" + formatted_charge.str() + "_Energy.dat";
        Copy_file(file_name,file_name_bipa);
    }

    //############################################################################################
    //Plot the scatter in the target
    Scatter(BlackPoints, RedPoints, secBeam, Mean_Energy);
    

    //############################################################################################
    //Plot the histogram of all secundary beam energies
    double bin, smaller, larger, maxy;
    if(Secundary_Beam_Energies.size() == 0){
        bin = 1;
        smaller = 1;
        larger = 1;
        maxy = 1;
    }else{
        bin = floor(11.5 + 2.8*pow(Secundary_Beam_Energies.size(),0.333333));
        smaller = floor(GetSmaller(Secundary_Beam_Energies));
        larger = ceil(GetLarger(Secundary_Beam_Energies));
        maxy = MoreRepeats(Secundary_Beam_Energies, smaller, larger, bin);
    }
    MakeHistogram(smaller, larger, maxy, bin, primaryTargThick, primaryBeam, primaryTarget, secBeam);
    //############################################################################################
    //Plot histogram of the energies on the target
    string energyTarget;
    vector<double> allEnerSecBeamOnTarg;
    ifstream SecundaryBeamEnergyOnTarget("Sec_Beam/OutputFiles/SecundaryBeamEnergyOnTarget.dat");
    if(!SecundaryBeamEnergyOnTarget){
        cerr << RED << "ERROR: " << RESET << "file not found: 'Sec_Beam/OutputFiles/SecundaryBeamEnergyOnTarget.dat'" << endl;
        exit(0);
    }
    for(int i = 0; getline(SecundaryBeamEnergyOnTarget, energyTarget); ++i){
    double EnergyValue;
    sscanf(energyTarget.c_str(),"%lf", &EnergyValue);
    allEnerSecBeamOnTarg.push_back(EnergyValue);
    }
    if(allEnerSecBeamOnTarg.size() != 0){
        double bintarget = floor(2.8*(pow(Secundary_Beam_Energies.size(),0.333333))+11.5);
        double smallertarget = floor(GetSmaller(Secundary_Beam_Energies));
        double largertarget = ceil(GetLarger(Secundary_Beam_Energies));
        double maxytarget = MoreRepeats(Secundary_Beam_Energies, smaller, larger, bin);
        MakeHistogramTarget(smallertarget, largertarget, maxytarget, bintarget, primaryTargThick, primaryBeam, primaryTarget, secBeam);
    }
    //############################################################################################
    //Histogram of the angles at which particles arrive at the target
    double AngMin, AngMax, binAng, angmaxy;
    if(Arrival_Angles_Target.size() == 0){
    AngMin = 0;
    AngMax = 0;
    binAng = 0;
    angmaxy = 0;
    }else{
    AngMin = floor(GetSmaller(Arrival_Angles_Target));
    AngMax = ceil(GetLarger(Arrival_Angles_Target));
    binAng = floor(2.8*(pow(Secundary_Beam_Energies.size(),0.333333))+11.5);
    angmaxy = MoreRepeats(Arrival_Angles_Target,AngMin,AngMax,binAng);
    }
    plotAngularDitribution(Arrival_Angles_Target, AngMin, AngMax, binAng, angmaxy);
    //############################################################################################
    //Coping the biparametric files from OutputFiles folder to the bipaFiles folder
    Copy_file("Sec_Beam/OutputFiles/anglesOnTarget.dat","Biparametric/Inputs/AnglesOnTarget.dat");
    Copy_file("Sec_Beam/OutputFiles/BeamDirectionOnTarget.dat","Biparametric/Inputs/BeamDirectionOnTarget.dat");
    Copy_file("Sec_Beam/OutputFiles/Points_on_Target.dat","Biparametric/Inputs/PointsOnTarget.dat");
    Copy_file("Sec_Beam/OutputFiles/SecundaryBeamEnergyOnTarget.dat","Biparametric/Inputs/SecEnergySecTarget.dat");

    //############################################################################################
    Program_Output(Mean_Energy, current, Mean_Energy_dispersion, Secundary_Beam_Energies_Target.size(), Scatter_Angles.size());

    //############################################################################################
    ClearDirectorySB();
}
        
