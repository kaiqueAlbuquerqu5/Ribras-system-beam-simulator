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

void write_input(const string filename, vector<string> sec_beam_list, vector<string> contaminant_list, vector<string> bipa_list, vector<string> geometric_list){
    /*This code has the function of writing the input file using the parameters provided in the input.dat file*/

    ofstream input_program;
    input_program.open(filename);
    //SecBeam input
    input_program << sec_beam_list[0] << endl;
    input_program << sec_beam_list[1] << "\t" << sec_beam_list[2] << "\t" << sec_beam_list[3] << "\t" << sec_beam_list[4] << "\t" << sec_beam_list[5] << "\t" << sec_beam_list[6] << "\t" << sec_beam_list[7] << "\t" << sec_beam_list[8] << "\t" << sec_beam_list[9] << "\t" << endl;
    input_program << geometric_list[0] << "\t" << geometric_list[1] << "\t" << geometric_list[2] << "\t" << geometric_list[3] << "\t" << geometric_list[4] << "\t" << geometric_list[5] << "\t" << geometric_list[6] << "\t" << geometric_list[7] << endl;
    input_program << geometric_list[8] << "\t" << geometric_list[9] << "\t" << "\n";
    input_program << sec_beam_list[10] << "\n";
    input_program << sec_beam_list[11] << "\n";
    input_program << sec_beam_list[12] << "\n" << "\n";

    //Contaminant input
    for(int i=0; i < contaminant_list.size(); ++i){
        input_program << contaminant_list[i] << endl;
    }
    input_program << "END" << endl << endl;
  

    //biparametric input
    input_program << bipa_list[0] << endl;
    input_program << bipa_list[1] << "\t" << bipa_list[4] << "\t" << bipa_list[5] << "\t" << bipa_list[6] << "\t" << bipa_list[2] << "\t" << bipa_list[8] << "\t" << bipa_list[3] << endl;

    vector<string> ejectiles, ejectiles_charges, recoils, Q_windows, Qgs_s;
    get_info_from_strig(bipa_list[9],ejectiles);
	for(int i=0; i< ejectiles.size(); ++i )
	{
		ejectiles_charges.push_back(sec_beam_list[6]);
	}
    get_info_from_strig(bipa_list[10],recoils);
    get_info_from_strig(bipa_list[11],Q_windows);
    get_info_from_strig(bipa_list[12],Qgs_s);
    for(int i = 0; i < ejectiles.size(); i++){
        input_program << sec_beam_list[3] << "\t" << bipa_list[7] << "\t" << ejectiles[i] << "\t" << recoils[i] << "\t" << Q_windows[i] << "\t" << Qgs_s[i] << "\t" << ejectiles_charges[i] << endl;
    }
    input_program << "END" << endl << endl; 


    //bipametric of contaminants
    for (int i = 0; i < contaminant_list.size(); i++){
        vector<string> Contaminant_info;
        get_info_from_strig(contaminant_list[i],Contaminant_info);
        input_program << Contaminant_info[0] << "\t" << bipa_list[7] << "\t" << Contaminant_info[0] << "\t" <<  bipa_list[7]  << "\t" << "0" << "\t" << "0" << "\t" << Contaminant_info[2] << endl;
    }
    input_program << "END";
}
    
    
int main(){
  vector<string> sec_beam_list, contaminant_list, bipa_list, geometric_list;
  get_parameters("input.dat",sec_beam_list, contaminant_list, bipa_list, geometric_list);
  write_input("Apps/in.dat", sec_beam_list, contaminant_list, bipa_list, geometric_list);
  system("./Apps/main < Apps/in.dat");
}
    
    
    
    
    
    
