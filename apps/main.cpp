#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>//used to convert stream
#include <unistd.h>//used to set directory
#include <math.h> 
#include <chrono>
#include <random>
#include <cstdlib>
#include "lib.hpp"
#define PI 3.14159265
#define RESET   "\033[0m"
#define BOLD    "\033[1m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"

using namespace std;



int main(){
    /*This code has the function of reading the input file written internally by the software and calling the subprograms (secbeam and biparametric)*/

    string Run_Sec_Beam;
    cin >> Run_Sec_Beam;
    long nVector;
    double primaryTargThick, massSecBeam, qSecBeam, Primary_Beam_Avarage_Energy, Qreac, ScatterAngMin, ScatterAngMax, FirstCollimatorRadius, FaradayCupRadius, LollipopRadius, SecondCollimatorRadius, SigmaPrimaryBeam, SigmaPrimaryEnergy, lolli2, Rloli2, currentFix;
    string primaryBeam, primaryTarget, secBeam, CrossSection, fixcurrent, straggling;
    cin >> nVector >> primaryBeam >> secBeam >> primaryTarget >> massSecBeam >> qSecBeam >> Qreac >> Primary_Beam_Avarage_Energy >> primaryTargThick >> ScatterAngMin >> ScatterAngMax >> SigmaPrimaryBeam >> SigmaPrimaryEnergy >> FirstCollimatorRadius >> FaradayCupRadius >> LollipopRadius >> SecondCollimatorRadius >> lolli2 >> Rloli2 >> CrossSection >> fixcurrent;
    if(primaryTargThick > 0){
        cout << endl << "********************************************************************"           << endl;
        cout <<BOLD<<YELLOW<< "  Warning: " << RESET << "The program will consider a cubic gaseous target with a  " << endl;
        cout <<         " 3.2 cm edge. In the case of a solid target, use the minus sign (-) " << endl;
        cout <<         "                   in front of the thickness in cm.                 " << endl;
        cout <<         "********************************************************************" << endl;
    }
    if(fixcurrent == "yes" || fixcurrent == "Yes" || fixcurrent == "YES" || fixcurrent == "y"){
        cin >> currentFix >> straggling;
    }else{
        currentFix = 0;
        cin >> straggling;
    }

    string Tag_Cont;
    double Mass_Cont, Charge_Cont;
    vector<string> Contamination_Tags;
    vector<double> Contamination_Mass, Contamination_Charge;
    while(true){
        cin >> Tag_Cont;
        char Break[10] = "END";
        if(strcmp(Tag_Cont.c_str(),Break)==0){
            break;
        }
        cin >> Mass_Cont >> Charge_Cont;
        Contamination_Tags.push_back(Tag_Cont);
        Contamination_Charge.push_back(Charge_Cont);
        Contamination_Mass.push_back(Mass_Cont);
    }
    
    if(Run_Sec_Beam == "yes" || Run_Sec_Beam == "Yes" || Run_Sec_Beam == "y"){
        Last_simulation('c',primaryBeam,secBeam,primaryTarget);
        secbeam(nVector, primaryTargThick, massSecBeam, qSecBeam, Primary_Beam_Avarage_Energy, Qreac, ScatterAngMin, ScatterAngMax, FirstCollimatorRadius, FaradayCupRadius, LollipopRadius, SecondCollimatorRadius, SigmaPrimaryBeam, SigmaPrimaryEnergy, lolli2, Rloli2, currentFix, primaryBeam, primaryTarget, secBeam, CrossSection, fixcurrent, straggling, Contamination_Tags, Contamination_Mass, Contamination_Charge);
    }else{
        Last_simulation('w',primaryBeam,secBeam,primaryTarget);
    }
    cout << "\n";

    /*
    Here the sec_beam ends and biparametric beguin.
    */
    double Detec_Angle, Detec_Size_x, Detec_Size_y, Detec_Size_z, Detec_Dist_Target, Detec_Radius, Sec_Targ_Thick;
    string Run_Biparametric, Detec_Tag, Detec_Shape;
    cin >> Run_Biparametric >> Detec_Tag;
    if(Detec_Tag == "END"){
        //If the user do not insert the first line of the biparametric input
        goto end;
    }
    cin >> Detec_Angle >> Detec_Shape;
    if(Detec_Shape == "circular" || Detec_Shape == "circ"){
        cin >> Detec_Radius;
        Detec_Size_x = 0, Detec_Size_y = 0;
    }else{
        cin >> Detec_Size_x >> Detec_Size_y;
        Detec_Radius = 0;
    }
    cin >> Detec_Size_z >> Sec_Targ_Thick >> Detec_Dist_Target;
    if(Run_Biparametric == "yes" || Run_Biparametric == "Yes" || Run_Biparametric == "y"){
        bipa(Detec_Angle, Detec_Radius, Detec_Size_x, Detec_Size_y, Detec_Size_z, Detec_Dist_Target, Detec_Tag, Sec_Targ_Thick);   
        bipa_contaminants(Detec_Angle, Detec_Radius, Detec_Size_x, Detec_Size_y, Detec_Size_z, Detec_Dist_Target, Detec_Tag, Sec_Targ_Thick);
    }
    
    end:
    cout << "********************************************************************"           << endl;
    cout << "          Program developed by OCBSantos and KAlbuquerque,          "           << endl;
    cout << "          using the twsp, cinema, kineq and stopx code, to          "           << endl;
    cout << "          simulate the secondary beam and the biparametrics         "           << endl;
    cout << "       in secondary scattering chamber of the RIBRAS facility.      "           << endl;
    cout << "********************************************************************"           << endl;
    cout << '\a';
    return 0;
}
