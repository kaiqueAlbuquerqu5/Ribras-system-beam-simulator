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

using namespace std;

#define PI 3.14159265

void Scatter(vector<vector<double>> BlackPoints, vector<vector<double>> RedPoints, string SecBeam, double Mean_Energy){
    ofstream scattergraf;
    scattergraf.open("Sec_Beam/OutputCharts/scattergraf.agr");
    scattergraf << "@version 50125" << endl;
    scattergraf << "@page size 792, 612" << endl;
    scattergraf << "@page scroll 5%" << endl;
    scattergraf << "@page inout 5%" << endl;
    scattergraf << "@link page off" << endl;
    scattergraf << "@default linewidth 1.0" << endl;
    scattergraf << "@default linestyle 1" << endl;
    scattergraf << "@map font 0 to \"Times-Roman\", \"Times-Roman\"" << endl;
    scattergraf << "@map color 0 to (255, 255, 255), \"white\"" << endl;
    scattergraf << "@map color 1 to (0, 0, 0), \"black\"" << endl;
    scattergraf << "@map color 2 to (255, 0, 0), \"red\"" << endl;
    scattergraf << "@map color 4 to (0, 0, 255), \"blue\"" << endl;
    scattergraf << "@map color 7 to (220, 220, 220), \"grey\"" << endl;
    scattergraf << "@map color 11 to (255, 165, 0), \"orange\"" << endl;
    scattergraf << "@g0 on" << endl;
    scattergraf << "@g0 hidden false" << endl;
    scattergraf << "@g0 type XY" << endl;
    scattergraf << "@with g0" << endl;
    scattergraf << "@    world -0.03, -0.03, 0.03, 0.03" << endl;
    scattergraf << "@    stack world 0, 0, 0, 0" << endl;
    scattergraf << "@    znorm 1" << endl;
    scattergraf << "@    view 0.115000, 0.080000, 0.925000, 0.890000" << endl;
    scattergraf << "@page background fill on" << endl;
    scattergraf << "@with line" << endl;
    scattergraf << "@    line on" << endl;
    scattergraf << "@    line loctype world" << endl;
    scattergraf << "@    line g0" << endl;
    scattergraf << "@    line 0.01, 0, 0, 0" << endl;
    scattergraf << "@    line linewidth 2.0" << endl;
    scattergraf << "@    line linestyle 4" << endl;
    scattergraf << "@    line color 4" << endl;
    scattergraf << "@    line arrow 0" << endl;
    scattergraf << "@    line arrow type 0" << endl;
    scattergraf << "@    line arrow length 1.000000" << endl;
    scattergraf << "@    line arrow layout 1.000000, 1.000000" << endl;
    scattergraf << "@line def" << endl;
    scattergraf << "@with string" << endl;
    scattergraf << "@    string on" << endl;
    scattergraf << "@    string loctype view" << endl;
    scattergraf << "@    string 0.333333333333, 0.94362745098" << endl;
    scattergraf << "@    string color 1" << endl;
    scattergraf << "@    string rot 0" << endl;
    scattergraf << "@    string font 0" << endl;
    scattergraf << "@    string just 0" << endl;
    scattergraf << "@    string char size 1.250000" << endl;
    scattergraf << "@    string def \"Dispersion of particles on target\"" << endl;
    scattergraf << "@with string" << endl;
    scattergraf << "@    string on" << endl;
    scattergraf << "@    string loctype world" << endl;
    scattergraf << "@    string g0" << endl;
    scattergraf << "@    string -0.011648511256, 0.031067538127" << endl;
    scattergraf << "@    string color 1" << endl;
    scattergraf << "@    string rot 0" << endl;
    scattergraf << "@    string font 0" << endl;
    scattergraf << "@    string just 0" << endl;
    scattergraf << "@    string char size 1.100000" << endl;
    scattergraf << "@    string def \"Energy On Target: " << Mean_Energy << " MeV\"" << endl;
    scattergraf << "@with string" << endl;
    scattergraf << "@    string on" << endl;
    scattergraf << "@    string loctype view" << endl;
    scattergraf << "@    string 1.04, 0.1" << endl;
    scattergraf << "@    string color 2" << endl;
    scattergraf << "@    string rot 0" << endl;
    scattergraf << "@    string font 0" << endl;
    scattergraf << "@    string just 0" << endl;
    scattergraf << "@    string char size 0.750000" << endl;
    scattergraf << "@    string def \"KAlbuquerque\"" << endl;
    scattergraf << "@with string" << endl;
    scattergraf << "@    string on" << endl;
    scattergraf << "@    string loctype view" << endl;
    scattergraf << "@    string 1.055, 0.069852941176" << endl;
    scattergraf << "@    string color 2" << endl;
    scattergraf << "@    string rot 0" << endl;
    scattergraf << "@    string font 0" << endl;
    scattergraf << "@    string just 0" << endl;
    scattergraf << "@    string char size 0.750000" << endl;
    scattergraf << "@    string def \"OCBSantos\"" << endl;
    scattergraf << "@with string" << endl;
    scattergraf << "@    string on" << endl;
    scattergraf << "@    string loctype world" << endl;
    scattergraf << "@    string g0" << endl;
    scattergraf << "@    string 0.0078, -0.001" << endl;
    scattergraf << "@    string color 4" << endl;
    scattergraf << "@    string rot 0" << endl;
    scattergraf << "@    string font 0" << endl;
    scattergraf << "@    string just 0" << endl;
    scattergraf << "@    string char size 0.600000" << endl;
    scattergraf << "@    string def \"1 cm\"" << endl;
    scattergraf << "@default sformat \"%.8g\"" << endl;
    scattergraf << "@background color 0" << endl;
    scattergraf << "@with g0" << endl;
    scattergraf << "@g0 hidden false" << endl;
    scattergraf << "@g0 type XY" << endl;
    scattergraf << "@    world -0.03, -0.03, 0.03, 0.03" << endl;
    scattergraf << "@    xaxes scale Normal" << endl;
    scattergraf << "@    yaxes scale Normal" << endl;
    scattergraf << "@    xaxis  on" << endl;
    scattergraf << "@    xaxis  label \"x (m)\"" << endl;
    scattergraf << "@    xaxis  label place auto" << endl;
    scattergraf << "@    xaxis  tick on" << endl;
    scattergraf << "@    xaxis  tick major 0.01" << endl;
    scattergraf << "@    xaxis  tick minor ticks 1" << endl;
    scattergraf << "@    xaxis  tick default 6" << endl;
    scattergraf << "@    yaxis  on" << endl;
    scattergraf << "@    yaxis  type zero false" << endl;
    scattergraf << "@    yaxis  label \"y (m)\"" << endl;
    scattergraf << "@    yaxis  label layout para" << endl;
    scattergraf << "@    yaxis  label place auto" << endl;
    scattergraf << "@    yaxis  label char size 1.000000" << endl;
    scattergraf << "@    yaxis  label font 0" << endl;
    scattergraf << "@    yaxis  label color 1" << endl;
    scattergraf << "@    yaxis  label place normal" << endl;
    scattergraf << "@    yaxis  tick on" << endl;
    scattergraf << "@    yaxis  tick major 0.01" << endl;
    scattergraf << "@    yaxis  tick minor ticks 1" << endl;
    scattergraf << "@    yaxis  tick default 6" << endl;
    scattergraf << "@    yaxis  tick place rounded true" << endl;
    scattergraf << "@    yaxis  tick in" << endl;
    scattergraf << "@    yaxis  tick major size 1.000000" << endl;
    scattergraf << "@    yaxis  tick major color 1" << endl;
    scattergraf << "@    yaxis  tick major linewidth 1.0" << endl;
    scattergraf << "@    yaxis  tick major linestyle 1" << endl;
    scattergraf << "@    yaxis  tick major grid off" << endl;
    scattergraf << "@    yaxis  ticklabel on" << endl;
    scattergraf << "@    yaxis  ticklabel format general" << endl;
    scattergraf << "@    yaxis  ticklabel prec 5" << endl;
    scattergraf << "@    yaxis  ticklabel place normal" << endl;
    scattergraf << "@    legend on" << endl;
    scattergraf << "@    legend loctype view" << endl;
    scattergraf << "@    legend 0.935, 0.89" << endl;
    scattergraf << "@    legend box color 1" << endl;
    scattergraf << "@    legend box pattern 1" << endl;
    scattergraf << "@    legend box linewidth 1.5" << endl;
    scattergraf << "@    legend box linestyle 1" << endl;
    scattergraf << "@    legend box fill color 0" << endl;
    scattergraf << "@    legend box fill pattern 1" << endl;
    scattergraf << "@    legend font 0" << endl;
    scattergraf << "@    legend char size 1.000000" << endl;
    scattergraf << "@    legend color 1" << endl;
    scattergraf << "@    legend length 4" << endl;
    scattergraf << "@    s0 hidden false" << endl;
    scattergraf << "@    s0 type xy" << endl;
    scattergraf << "@    s0 line type 2" << endl;
    scattergraf << "@    s0 line linestyle 1" << endl;
    scattergraf << "@    s0 line linewidth 2.5" << endl;
    scattergraf << "@    s0 line color 7" << endl;
    scattergraf << "@    s0 line pattern 1" << endl;
    scattergraf << "@    s0 fill type 1" << endl;
    scattergraf << "@    s0 fill rule 0" << endl;
    scattergraf << "@    s0 fill color 7" << endl;
    scattergraf << "@    s0 fill pattern 1" << endl;
    scattergraf << "@    s0 legend  \"Frame\"" << endl;
    scattergraf << "@    s1 hidden false" << endl;
    scattergraf << "@    s1 type xy" << endl;
    scattergraf << "@    s1 symbol 0" << endl;
    scattergraf << "@    s1 symbol size 1.000000" << endl;
    scattergraf << "@    s1 symbol color 7" << endl;
    scattergraf << "@    s1 symbol pattern 1" << endl;
    scattergraf << "@    s1 symbol fill color 11" << endl;
    scattergraf << "@    s1 symbol fill pattern 0" << endl;
    scattergraf << "@    s1 symbol linewidth 1.0" << endl;
    scattergraf << "@    s1 symbol linestyle 1" << endl;
    scattergraf << "@    s1 symbol char 65" << endl;
    scattergraf << "@    s1 line type 1" << endl;
    scattergraf << "@    s1 line linestyle 1" << endl;
    scattergraf << "@    s1 line linewidth 2.5" << endl;
    scattergraf << "@    s1 line color 7" << endl;
    scattergraf << "@    s1 fill type 1" << endl;
    scattergraf << "@    s1 fill rule 0" << endl;
    scattergraf << "@    s1 fill color 0" << endl;
    scattergraf << "@    s1 fill pattern 1" << endl;
    scattergraf << "@    s2 hidden false" << endl;
    scattergraf << "@    s2 type xy" << endl;
    scattergraf << "@    s2 symbol 1" << endl;
    scattergraf << "@    s2 symbol size 0.090000" << endl;
    scattergraf << "@    s2 symbol color 1" << endl;
    scattergraf << "@    s2 symbol pattern 1" << endl;
    scattergraf << "@    s2 symbol fill color 1" << endl;
    scattergraf << "@    s2 symbol fill pattern 14" << endl;
    scattergraf << "@    s2 symbol linewidth 2.5" << endl;
    scattergraf << "@    s2 symbol linestyle 1" << endl;
    scattergraf << "@    s2 symbol char 65" << endl;
    scattergraf << "@    s2 line type 0" << endl;
    scattergraf << "@    s2 legend  \"" << SecBeam << "\"" << endl;
    scattergraf << "@    s3 hidden false" << endl;
    scattergraf << "@    s3 type xy" << endl;
    scattergraf << "@    s3 symbol 1" << endl;
    scattergraf << "@    s3 symbol size 0.090000" << endl;
    scattergraf << "@    s3 symbol color 2" << endl;
    scattergraf << "@    s3 symbol pattern 1" << endl;
    scattergraf << "@    s3 symbol fill color 1" << endl;
    scattergraf << "@    s3 symbol fill pattern 14" << endl;
    scattergraf << "@    s3 symbol linewidth 2.5" << endl;
    scattergraf << "@    s3 symbol linestyle 1" << endl;
    scattergraf << "@    s3 symbol char 65" << endl;
    scattergraf << "@    s3 line type 0" << endl;
    scattergraf << "@    s3 legend  \"*" << SecBeam << "\"" << endl;
    scattergraf << "@target G0.S0" << endl;
    scattergraf << "@type xy" << endl;
    double twopi = 6.283;
    double contadorFrame;
    while(contadorFrame < twopi){
        double x = 0.021*cos(contadorFrame);
        double y = 0.021*sin(contadorFrame);
        scattergraf << x << " " << y << endl;
        contadorFrame += 0.001;
    }
    scattergraf << "&" << endl;
    scattergraf << "@target G0.S1" << endl;
    scattergraf << "@type xy" << endl;
    double contadorTarget = 0;
    while(contadorTarget < twopi){
        double x = 0.01*cos(contadorTarget);
        double y = 0.01*sin(contadorTarget);
        scattergraf << x << " " << y << endl;
        contadorTarget += 0.001;
    }
    scattergraf << "&" << endl;
    scattergraf << "@target G0.S2" << endl;
    scattergraf << "@type xy" << endl;
    for(long i = 0; i < BlackPoints.size(); ++i){
        for(long j = 0; j < 2; ++j){
            scattergraf << BlackPoints[i][j] << " ";
        }
        scattergraf << endl;
	}
    scattergraf << "&" << endl;
    scattergraf << "@target G0.S3" << endl;
    scattergraf << "@type xy" << endl;
    
    for(long i = 0; i < RedPoints.size(); ++i){
        for(long j = 0; j < 2; ++j){
            scattergraf << RedPoints[i][j] << " ";
        }
        scattergraf << endl;
	}
    
    scattergraf.close();
}
