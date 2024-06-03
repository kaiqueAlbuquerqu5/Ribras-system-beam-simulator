#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>  //used to convert stream
#include <unistd.h> //used to set directory
#include <math.h>
#include <chrono>
#include <random>
#include <cstdlib>
#define RESET   "\033[0m"
#define BOLD    "\033[1m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define RED     "\033[31m"

using namespace std;

void MakeHistogram(double smaller, double larger, double maxy, double bin, double primaryTargThick, string primaryBeam, string primaryTarget, string secBeam)
{
    double auxThick;
    string unitThick;

    if (primaryTargThick > 0)
    {
        unitThick = "mg/cm2";
        auxThick = primaryTargThick;
    }
    else
    {
        unitThick = "um";
        auxThick = primaryTargThick * (-10000);
    }

    ofstream histogram;
    histogram.open("Sec_Beam/OutputCharts/histogram.agr");
    histogram << "@version 50125" << endl;
    histogram << "@page size 792, 612" << endl;
    histogram << "@page scroll 5%" << endl;
    histogram << "@page inout 5%" << endl;
    histogram << "@link page off" << endl;
    histogram << "@default linewidth 1.0" << endl;
    histogram << "@default linestyle 1" << endl;
    histogram << "@g0 on" << endl;
    histogram << "@g0 hidden false" << endl;
    histogram << "@g0 type XY" << endl;
    histogram << "@with g0" << endl;
    histogram << "@map font 0 to \"Times-Roman\", \"Times-Roman\"" << endl;
    histogram << "@map color 0 to (255, 255, 255), \"white\"" << endl;
    histogram << "@map color 2 to (255, 0, 0), \"red\"" << endl;
    histogram << "@map color 4 to (0, 0, 255), \"blue\"" << endl;
    histogram << "@page background fill on" << endl;
    histogram << "@r0 off" << endl;
    histogram << "@with string" << endl;
    histogram << "@    string on" << endl;
    histogram << "@    string loctype view" << endl;
    histogram << "@    string 1.1, 0.03" << endl;
    histogram << "@    string color 2" << endl;
    histogram << "@    string rot 0" << endl;
    histogram << "@    string font 0" << endl;
    histogram << "@    string just 0" << endl;
    histogram << "@    string char size 0.800000" << endl;
    histogram << "@    string def \"OCBSantos\"" << endl;
    histogram << "@with string" << endl;
    histogram << "@    string on" << endl;
    histogram << "@    string loctype view" << endl;
    histogram << "@    string 1.083, 0.066" << endl;
    histogram << "@    string color 2" << endl;
    histogram << "@    string rot 0" << endl;
    histogram << "@    string font 0" << endl;
    histogram << "@    string just 0" << endl;
    histogram << "@    string char size 0.800000" << endl;
    histogram << "@    string def \"KAlbuquerque\"" << endl;
    histogram << "@    world " << smaller << ", "
              << "0, " << larger << ", " << maxy + 2 << endl;
    histogram << "@    stack world 0, 0, 0, 0" << endl;
    histogram << "@    view 0.150000, 0.150000, 1.150000, 0.850000" << endl;
    histogram << "@    znorm 1" << endl;
    histogram << "@    view 0.150000, 0.150000, 1.150000, 0.850000" << endl;
    histogram << "@    title \"Secundary beam energies  \"" << endl;
    histogram << "@    title font 0" << endl;
    histogram << "@    title size 1.700000" << endl;
    histogram << "@    title color 1" << endl;
    histogram << "@    subtitle \"Primary beam: " << primaryBeam << "; Primary target: " << primaryTarget << " whit thickness = " << fixed << setprecision(2) << auxThick << " (" << unitThick << "); Secundary beam: " << secBeam << "\"" << endl;
    histogram << "@    subtitle font 0" << endl;
    histogram << "@    subtitle size 1.000000" << endl;
    histogram << "@    subtitle color 1" << endl;
    histogram << "@    xaxes scale Normal" << endl;
    histogram << "@    xaxis  on" << endl;
    histogram << "@    xaxis  bar color 1" << endl;
    histogram << "@    xaxis  bar linestyle 1" << endl;
    histogram << "@    xaxis  bar linewidth 1.0" << endl;
    histogram << "@    xaxis  label \"Energy (MeV)\"" << endl;
    histogram << "@    xaxis  label color 1" << endl;
    histogram << "@    xaxis  tick minor size 0.500000" << endl;
    histogram << "@    xaxis  ticklabel on" << endl;
    histogram << "@    yaxes scale Normal" << endl;
    histogram << "@    yaxis  on" << endl;
    histogram << "@    yaxis  bar color 1" << endl;
    histogram << "@    yaxis  bar linestyle 1" << endl;
    histogram << "@    yaxis  bar linewidth 1.0" << endl;
    histogram << "@    yaxis  label \"Count\"" << endl;
    histogram << "@    yaxis  label color 1" << endl;
    histogram << "@    yaxis  tick major 1" << endl;
    histogram << "@    yaxis  tick minor size 1.000000" << endl;
    if(maxy < 30){
        histogram << "@    yaxis  tick major 1" << endl;
    }
    else if(maxy < 70 && maxy >= 30){
        histogram << "@    yaxis  tick major 3" << endl;
    }
    else if(maxy >=70){
        histogram << "@    yaxis  tick major 6" << endl;
    }
    histogram << "@    yaxis  ticklabel on" << endl;
    histogram << "@    s1 hidden false" << endl;
    histogram << "@    s1 type xy" << endl;
    histogram << "@    s1 symbol 0" << endl;
    histogram << "@    s1 symbol size 0.950000" << endl;
    histogram << "@    s1 symbol pattern 1" << endl;
    histogram << "@    s1 symbol color 4" << endl;
    histogram << "@    s1 symbol fill color 4" << endl;
    histogram << "@    s1 symbol linewidth 1.0" << endl;
    histogram << "@    s1 symbol char 65" << endl;
    histogram << "@    s1 line type 3" << endl;
    histogram << "@    s1 line linestyle 1" << endl;
    histogram << "@    s1 line linewidth 1.5" << endl;
    histogram << "@    s1 line color 4" << endl;
    histogram << "@    s1 fill type 1" << endl;
    histogram << "@    s1 fill rule 0" << endl;
    histogram << "@    s1 fill color 4" << endl;
    histogram << "@    s1 fill pattern 14" << endl;
    histogram << "@    s1 comment \"Histogram from G0.S0\"" << endl;
    histogram << "@    s1 legend  \"\"" << endl;
    histogram << "@target G0.S0" << endl;
    histogram << "@type xy" << endl;

    string energy;
    vector<double> VectorEnergy;
    ifstream energies("Sec_Beam/OutputFiles/SecundaryBeamEnergiesSelected.dat");
    if(!energies){
      cerr << RED << "ERROR: " << RESET << "file not found: 'Sec_Beam/OutputFiles/SecundaryBeamEnergiesSelected.dat'" << endl;
      exit(0);
    } 
    long j;
    for (j = 0; getline(energies, energy); ++j)
    {
        double value;
        sscanf(energy.c_str(), "%lf", &value);
        VectorEnergy.push_back(value);
    }
    energies.close();

    histogram << "&" << endl;
    histogram << "@target G0.S1" << endl;
    histogram << "@type xy" << endl;

    double step = (larger - smaller) / bin;
    while (smaller < larger)
    {
        long n = 0;
        for (long j = 0; j < VectorEnergy.size(); ++j)
        {
            if (VectorEnergy[j] >= smaller && VectorEnergy[j] < smaller + step)
            {
                ++n;
            }
        }
        histogram << fixed << setprecision(5) << smaller - (step / 2) << " " << n << endl;
        smaller += step;
    }
    histogram.close();

    //remove("SecundaryEnergies.dat");
    remove("secBeamOut.dat");
}

void MakeHistogramTarget(double smallertarget, double largertarget, double maxytarget, double bintarget, double primaryTargThick, string primaryBeam, string primaryTarget, string secBeam)
{
    double auxThick;
    string unitThick;

    if (primaryTargThick > 0)
    {
        unitThick = "mg/cm2";
        auxThick = primaryTargThick;
    }
    else
    {
        unitThick = "um";
        auxThick = primaryTargThick * (-10000);
    }

    ofstream histogramtarget;
    histogramtarget.open("Sec_Beam/OutputCharts/histogramtarget.agr");
    histogramtarget << "@version 50125" << endl;
    histogramtarget << "@page size 792, 612" << endl;
    histogramtarget << "@page scroll 5%" << endl;
    histogramtarget << "@page inout 5%" << endl;
    histogramtarget << "@link page off" << endl;
    histogramtarget << "@default linewidth 1.0" << endl;
    histogramtarget << "@default linestyle 1" << endl;
    histogramtarget << "@g0 on" << endl;
    histogramtarget << "@g0 hidden false" << endl;
    histogramtarget << "@g0 type XY" << endl;
    histogramtarget << "@with g0" << endl;
    histogramtarget << "@map font 0 to \"Times-Roman\", \"Times-Roman\"" << endl;
    histogramtarget << "@map color 0 to (255, 255, 255), \"white\"" << endl;
    histogramtarget << "@map color 2 to (255, 0, 0), \"red\"" << endl;
    histogramtarget << "@map color 4 to (0, 0, 255), \"blue\"" << endl;
    histogramtarget << "@page background fill on" << endl;
    histogramtarget << "@r0 off" << endl;
    histogramtarget << "@with string" << endl;
    histogramtarget << "@    string on" << endl;
    histogramtarget << "@    string loctype view" << endl;
    histogramtarget << "@    string 1.1, 0.03" << endl;
    histogramtarget << "@    string color 2" << endl;
    histogramtarget << "@    string rot 0" << endl;
    histogramtarget << "@    string font 0" << endl;
    histogramtarget << "@    string just 0" << endl;
    histogramtarget << "@    string char size 0.800000" << endl;
    histogramtarget << "@    string def \"OCBSantos\"" << endl;
    histogramtarget << "@with string" << endl;
    histogramtarget << "@    string on" << endl;
    histogramtarget << "@    string loctype view" << endl;
    histogramtarget << "@    string 1.083, 0.066" << endl;
    histogramtarget << "@    string color 2" << endl;
    histogramtarget << "@    string rot 0" << endl;
    histogramtarget << "@    string font 0" << endl;
    histogramtarget << "@    string just 0" << endl;
    histogramtarget << "@    string char size 0.800000" << endl;
    histogramtarget << "@    string def \"KAlbuquerque\"" << endl;
    histogramtarget << "@    world " << smallertarget << ", "
              << "0, " << largertarget << ", " << maxytarget + 2 << endl;
    histogramtarget << "@    stack world 0, 0, 0, 0" << endl;
    histogramtarget << "@    view 0.150000, 0.150000, 1.150000, 0.850000" << endl;
    histogramtarget << "@    znorm 1" << endl;
    histogramtarget << "@    view 0.150000, 0.150000, 1.150000, 0.850000" << endl;
    histogramtarget << "@    title \"Secundary beam energies  \"" << endl;
    histogramtarget << "@    title font 0" << endl;
    histogramtarget << "@    title size 1.700000" << endl;
    histogramtarget << "@    title color 1" << endl;
    histogramtarget << "@    subtitle \"Primary beam: " << primaryBeam << "; Primary target: " << primaryTarget << " whit thickness = " << fixed << setprecision(2) << auxThick << " (" << unitThick << "); Secundary beam: " << secBeam << "\"" << endl;
    histogramtarget << "@    subtitle font 0" << endl;
    histogramtarget << "@    subtitle size 1.000000" << endl;
    histogramtarget << "@    subtitle color 1" << endl;
    histogramtarget << "@    xaxes scale Normal" << endl;
    histogramtarget << "@    xaxis  on" << endl;
    histogramtarget << "@    xaxis  bar color 1" << endl;
    histogramtarget << "@    xaxis  bar linestyle 1" << endl;
    histogramtarget << "@    xaxis  bar linewidth 1.0" << endl;
    histogramtarget << "@    xaxis  label \"Energy (MeV)\"" << endl;
    histogramtarget << "@    xaxis  label color 1" << endl;
    histogramtarget << "@    xaxis  tick minor size 0.500000" << endl;
    histogramtarget << "@    xaxis  ticklabel on" << endl;
    histogramtarget << "@    yaxes scale Normal" << endl;
    histogramtarget << "@    yaxis  on" << endl;
    histogramtarget << "@    yaxis  bar color 1" << endl;
    histogramtarget << "@    yaxis  bar linestyle 1" << endl;
    histogramtarget << "@    yaxis  bar linewidth 1.0" << endl;
    histogramtarget << "@    yaxis  label \"Count\"" << endl;
    histogramtarget << "@    yaxis  label color 1" << endl;
    histogramtarget << "@    yaxis  tick minor size 1.000000" << endl;
    if(maxytarget < 30){
        histogramtarget << "@    yaxis  tick major 1" << endl;
    }
    else if(maxytarget < 70 && maxytarget >= 30){
        histogramtarget << "@    yaxis  tick major 3" << endl;
    }
    else if(maxytarget >=70){
        histogramtarget << "@    yaxis  tick major 6" << endl;
    }
    histogramtarget << "@    yaxis  ticklabel on" << endl;
    histogramtarget << "@    s1 hidden false" << endl;
    histogramtarget << "@    s1 type xy" << endl;
    histogramtarget << "@    s1 symbol 0" << endl;
    histogramtarget << "@    s1 symbol size 0.950000" << endl;
    histogramtarget << "@    s1 symbol pattern 1" << endl;
    histogramtarget << "@    s1 symbol color 4" << endl;
    histogramtarget << "@    s1 symbol fill color 4" << endl;
    histogramtarget << "@    s1 symbol linewidth 1.0" << endl;
    histogramtarget << "@    s1 symbol char 65" << endl;
    histogramtarget << "@    s1 line type 3" << endl;
    histogramtarget << "@    s1 line linestyle 1" << endl;
    histogramtarget << "@    s1 line linewidth 1.5" << endl;
    histogramtarget << "@    s1 line color 4" << endl;
    histogramtarget << "@    s1 fill type 1" << endl;
    histogramtarget << "@    s1 fill rule 0" << endl;
    histogramtarget << "@    s1 fill color 4" << endl;
    histogramtarget << "@    s1 fill pattern 14" << endl;
    histogramtarget << "@    s1 comment \"Histogram from G0.S0\"" << endl;
    histogramtarget << "@    s1 legend  \"\"" << endl;
    histogramtarget << "@target G0.S0" << endl;
    histogramtarget << "@type xy" << endl;
    string energy;
    vector<double> VectorEnergy;
    ifstream energies("Sec_Beam/OutputFiles/SecundaryBeamEnergyOnTarget.dat");
    if(!energies){
      cerr << RED << "ERROR: " << RESET << "file not found: 'Sec_Beam/OutputFiles/SecundaryBeamEnergyOnTarget.dat'" << endl;
      exit(0);
    } 
    for (long j = 0; getline(energies, energy); ++j)
    {
        double value;
        sscanf(energy.c_str(), "%lf", &value);
        VectorEnergy.push_back(value);
    }
    energies.close();

    histogramtarget << "&" << endl;
    histogramtarget << "@target G0.S1" << endl;
    histogramtarget << "@type xy" << endl;

    double step = (largertarget - smallertarget) / bintarget;
    while (smallertarget < largertarget)
    {
        long n = 0;
        for (long j = 0; j < VectorEnergy.size(); ++j)
        {
            if (VectorEnergy[j] >= smallertarget && VectorEnergy[j] < smallertarget + step)
            {
                ++n;
            }
        }
        histogramtarget << fixed << setprecision(5) << smallertarget - (step / 2) << " " << n << endl;
        smallertarget += step;
    }
    histogramtarget.close();
//    remove("SecEnergySecTarget.dat");
}


void plotAngularDitribution(vector<double> angles, double AngleMin, double AngleMax, double binAngle, double angmaxy){
    ofstream AngDistTarget;
    AngDistTarget.open("Sec_Beam/OutputCharts/AnglesTarget.agr");
    AngDistTarget << "@version 50125" << endl;
    AngDistTarget << "@page size 792, 612" << endl;
    AngDistTarget << "@page scroll 5%" << endl;
    AngDistTarget << "@page inout 5%" << endl;
    AngDistTarget << "@link page off" << endl;
    AngDistTarget << "@default linewidth 1.0" << endl;
    AngDistTarget << "@default linestyle 1" << endl;
    AngDistTarget << "@map font 0 to \"Times-Roman\", \"Times-Roman\"" << endl;
    AngDistTarget << "@map color 0 to (255, 255, 255), \"white\"" << endl;
    AngDistTarget << "@map color 1 to (0, 0, 0), \"black\"" << endl;
    AngDistTarget << "@map color 2 to (255, 0, 0), \"red\"" << endl;
    AngDistTarget << "@map color 4 to (0, 0, 255), \"blue\"" << endl;
    AngDistTarget << "@map color 7 to (220, 220, 220), \"grey\"" << endl;
    AngDistTarget << "@map color 11 to (255, 165, 0), \"orange\"" << endl;
    AngDistTarget << "@g0 on" << endl;
    AngDistTarget << "@g0 hidden false" << endl;
    AngDistTarget << "@g0 type XY" << endl;
    AngDistTarget << "@with g0" << endl;
    AngDistTarget << "@    world 2, 0, 5.3, ";
    if(angmaxy < 30){
        AngDistTarget << angmaxy + 2 << endl;
    }
    else if(angmaxy < 70 && angmaxy >= 30){
        AngDistTarget << angmaxy+5 << endl;
    }
    else if(angmaxy >=70){
        AngDistTarget << angmaxy+10 << endl;
    }
    AngDistTarget << "@    stack world 0, 0, 0, 0" << endl;
    AngDistTarget << "@    znorm 1" << endl;
    AngDistTarget << "@    view 0.15000, 0.150000, 1.15000, 0.850000" << endl;
    AngDistTarget << "@page background fill on" << endl;
    AngDistTarget << "@with string" << endl;
    AngDistTarget << "@    string on" << endl;
    AngDistTarget << "@    string loctype view" << endl;
    AngDistTarget << "@    string 0.29, 0.936274509804" << endl;
    AngDistTarget << "@    string color 1" << endl;
    AngDistTarget << "@    string rot 0" << endl;
    AngDistTarget << "@    string font 0" << endl;
    AngDistTarget << "@    string just 0" << endl;
    AngDistTarget << "@    string char size 1.50000" << endl;
    AngDistTarget << "@    string def \"Counting the trajectory angles on the target\"" << endl;
    AngDistTarget << "@with string" << endl;
    AngDistTarget << "@    string on" << endl;
    AngDistTarget << "@    string loctype view" << endl;
    AngDistTarget << "@    string 1.08455882353, 0.080882352941" << endl;
    AngDistTarget << "@    string color 2" << endl;
    AngDistTarget << "@    string rot 0" << endl;
    AngDistTarget << "@    string font 0" << endl;
    AngDistTarget << "@    string just 0" << endl;
    AngDistTarget << "@    string char size 0.700000" << endl;
    AngDistTarget << "@    string def \"KAlbuquerque\"" << endl;
    AngDistTarget << "@with string" << endl;
    AngDistTarget << "@    string on" << endl;
    AngDistTarget << "@    string loctype view" << endl;
    AngDistTarget << "@    string 1.084, 0.054852941176" << endl;
    AngDistTarget << "@    string color 2" << endl;
    AngDistTarget << "@    string rot 0" << endl;
    AngDistTarget << "@    string font 0" << endl;
    AngDistTarget << "@    string just 0" << endl;
    AngDistTarget << "@    string char size 0.750000" << endl;
    AngDistTarget << "@    string def \"OCBSantos\"" << endl;
    AngDistTarget << "@with string" << endl;
    AngDistTarget << "@    string on" << endl;
    AngDistTarget << "@    string loctype world" << endl;
    AngDistTarget << "@    string g0" << endl;
    AngDistTarget << "@    string 1.085, 0.057598039216" << endl;
    AngDistTarget << "@    string color 4" << endl;
    AngDistTarget << "@    string rot 0" << endl;
    AngDistTarget << "@    string font 0" << endl;
    AngDistTarget << "@    string just 0" << endl;
    AngDistTarget << "@    string char size 0.700000" << endl;
    AngDistTarget << "@    string def \"1 cm\"" << endl;
    AngDistTarget << "@default sformat \"%.8g\"" << endl;
    AngDistTarget << "@background color 0" << endl;
    AngDistTarget << "@with g0" << endl;
    AngDistTarget << "@g0 hidden false" << endl;
    AngDistTarget << "@g0 type XY" << endl;
    AngDistTarget << "@    xaxes scale Normal" << endl;
    AngDistTarget << "@    yaxes scale Normal" << endl;
    AngDistTarget << "@    xaxis  on" << endl;
    AngDistTarget << "@    xaxis  label \"Trajectory angle on target\"" << endl;
    AngDistTarget << "@    xaxis  label place auto" << endl;
    AngDistTarget << "@    xaxis  tick on" << endl;
    AngDistTarget << "@    xaxis  tick major 0.5" << endl;
    AngDistTarget << "@    xaxis  tick minor ticks 1" << endl;
    AngDistTarget << "@    xaxis  tick default 6" << endl;
    AngDistTarget << "@    yaxis  on" << endl;
    AngDistTarget << "@    yaxis  type zero false" << endl;
    AngDistTarget << "@    yaxis  label \"Count\"" << endl;
    AngDistTarget << "@    yaxis  label layout para" << endl;
    AngDistTarget << "@    yaxis  label place auto" << endl;
    AngDistTarget << "@    yaxis  label char size 1.000000" << endl;
    AngDistTarget << "@    yaxis  label font 0" << endl;
    AngDistTarget << "@    yaxis  label color 1" << endl;
    AngDistTarget << "@    yaxis  label place normal" << endl;
    AngDistTarget << "@    yaxis  tick on" << endl;
    if(angmaxy < 30){
        AngDistTarget << "@    yaxis  tick major 1" << endl;
    }
    else if(angmaxy < 70 && angmaxy >= 30){
        AngDistTarget << "@    yaxis  tick major 3" << endl;
    }
    else if(angmaxy >=70){
        AngDistTarget << "@    yaxis  tick major 6" << endl;
    }
    AngDistTarget << "@    yaxis  tick minor ticks 1" << endl;
    AngDistTarget << "@    yaxis  tick default 6" << endl;
    AngDistTarget << "@    yaxis  tick place rounded true" << endl;
    AngDistTarget << "@    yaxis  tick in" << endl;
    AngDistTarget << "@    yaxis  tick major size 1.000000" << endl;
    AngDistTarget << "@    yaxis  tick major color 1" << endl;
    AngDistTarget << "@    yaxis  tick major linewidth 1.0" << endl;
    AngDistTarget << "@    yaxis  tick major linestyle 1" << endl;
    AngDistTarget << "@    yaxis  tick major grid off" << endl;
    AngDistTarget << "@    yaxis  ticklabel on" << endl;
    AngDistTarget << "@    yaxis  ticklabel format general" << endl;
    AngDistTarget << "@    yaxis  ticklabel prec 5" << endl;
    AngDistTarget << "@    yaxis  ticklabel place normal" << endl;
    AngDistTarget << "@    legend on" << endl;
    AngDistTarget << "@    legend loctype view" << endl;
    AngDistTarget << "@    legend 0.85, 0.8" << endl;
    AngDistTarget << "@    legend box color 1" << endl;
    AngDistTarget << "@    legend box pattern 1" << endl;
    AngDistTarget << "@    legend box linewidth 1.5" << endl;
    AngDistTarget << "@    legend box linestyle 1" << endl;
    AngDistTarget << "@    legend box fill color 0" << endl;
    AngDistTarget << "@    legend box fill pattern 1" << endl;
    AngDistTarget << "@    legend font 0" << endl;
    AngDistTarget << "@    legend char size 1.000000" << endl;
    AngDistTarget << "@    legend color 1" << endl;
    AngDistTarget << "@    legend length 4" << endl;
    AngDistTarget << "@    s0 hidden false" << endl;
    AngDistTarget << "@    s0 type xy" << endl;
    AngDistTarget << "@    s0 line type 2" << endl;
    AngDistTarget << "@    s0 line linestyle 1" << endl;
    AngDistTarget << "@    s0 line linewidth 2.0" << endl;
    AngDistTarget << "@    s0 line color 1" << endl;
    AngDistTarget << "@    s0 line pattern 1" << endl;
    AngDistTarget << "@    s0 fill type 2" << endl;
    AngDistTarget << "@    s0 fill rule 0" << endl;
    AngDistTarget << "@    s0 fill color 1" << endl;
    AngDistTarget << "@    s0 fill pattern 6" << endl;
    AngDistTarget << "@target G0.S0" << endl;
    AngDistTarget << "@type xy" << endl;
    double step = (AngleMax - AngleMin)/binAngle;
    while (AngleMin < AngleMax){
        long n = 0;
        for (long j = 0; j < angles.size(); ++j)
        {
            if (angles[j] >= AngleMin && angles[j] < AngleMin + step)
            {
                ++n;
            }
        }
        AngDistTarget << fixed << setprecision(5) << AngleMin + (step)/2 << "   " << n << endl;
        AngleMin += step;
    }
    AngDistTarget << "&" << endl;
}
