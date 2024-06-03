#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>

#ifndef _H_TESTE
#define _H_TESTE
using namespace std;

string PointCheck(double x, double y, vector<vector<double>> data);
void ClearDirectorySB();
void get_info_from_strig(string text, vector<string> &informations);
void get_parameters(const string filename, vector<string> &sec_beam_list, vector<string> &contaminant_list, vector<string> &bipa_list, vector<string> &geometric_list);
void Last_simulation(char opt, string pri_beam, string sec_beam, string targ);
void Energy_Loss_stopx(double Thickness, double E_in, string Beam, string Target, string Straggling, double &E_out, double &sigE_dx1);
void Energy_Loss_Reaction_SB(double E_in, double Theta, string Primary_Beam, string Target, string Secundary_Beam, double Qreac, double &Energy_out);
void Energy_Loss_Reaction_Bipa(double E_in, double Theta, string Primary_Beam, string Target, string Sec_Beam, string Recoil, double Qreac, double &E_out_Sec_Beam, double &E_out_Recoil, double &Recoil_Angle);
void Current_Calculation(string secBeam, string fixcurrent, double currentFix, double massSecBeam, double Avarage_Energy, double qSecBeam, double &current);
void Read_data(vector<vector<double> > &Vdata, string data);
void maxValue(vector<vector<double>> Vdata, double &max);
void interpolation(vector<vector<double>> Vdata, double xint, double &fxint);
void Load_bar(double progress);
void Particle_Mass(string particle_tag, int &particle_mass);
void Primary_Target_Events(long nVector, double primaryTargThick, string CrossSection, double SigmaPrimaryBeam, double AngMin, double AngMax, vector<double> &Scatter_Angles, vector<double> &Phi_Angles, vector<vector<double>> &Reaction_Points, double &Nevents);
void run_stopx();
void run_kineq();
void run_cinema();
void run_twsp();
void run_strag();
void Input_kineq(string target, string beam, string eject, double energy, double angle, double Q);
void cinema_input(string Proj_Tag, string Targ_Tag, string Eject_Tag, string Recoil_Tag, double Energy, double Angle, double Qreac);
void read_cinema(double &En3, double &THLB4, double &En4);
void Input_stopx(string target, string beam, double energy, double thickness);
void Input_strag(string projectile, string target, double thickness, double energy);
void Input_twsp(double current, double phi1, vector<double> Reaction_Point, double primaryTargThick, double qSecBeam, double massSecBeam, double EnerSecBeam, double Scatter_Angle);
void plotGraf(long nInteractions, double primaryTargThick, double Qreac, double current, vector<vector<double>> Reaction_Points, string secBeam, string primaryBeam, string primaryTarget, double alpha, double energyPrimaryBeam, double qSecBeam, double massSecBeam, double FirstCollimatorRadius, double FardayCupRadius, double LollipopRadius, double SecondCollimatorRadius, vector<double> vecEnerSecBeam, vector<double> Scatter_Angles_Temp, vector<double> Phi_Angles_Temp, vector<double> &VectorEnergySecTarget, vector<double> &Scatter_Angles, vector<double> &Phi_Angles, vector<double> &Sec_Beam_Energies, double lolli2, double Rlolli2, vector<vector<double>> &BlackPoints, vector<vector<double>> &RedPoints, vector<double> &Arrival_Angles_Target);
void plotAngularDitribution(vector<double> angles, double AngleMin, double AngleMax, double binAngle, double angmaxy);
void MakeHistogram(double smaller, double larger, double maxy, double bin, double primaryTargThick, string primaryBeam, string primaryTarget, string secBeam);
void MakeHistogramTarget(double smallertarget, double largertarget, double maxytarget, double bintarget, double primaryTargThick, string primaryBeam, string primaryTarget, string secBeam);
void Scatter(vector<vector<double>> BlackPoints, vector<vector<double>> RedPoints, string secBeam, double Mean_Energy);
void Copy_file(string file_origin, string file_target);
void Program_Output(double Mean_Energy, double current, double Mean_energy_uncertainty, int NParticles_Target, int NParticles_Selected);
void Read_Vector(vector<double> &vec, string data);
void Read_Matrix(vector<vector<double>> &Vdata, string data);
void Clear_directory_Bipa();
void Save_Vector_Data(vector<double> V, string file_name);
void Reaction_Point(vector<double> &R_point, vector<double> U_vector, vector<double> I_point, double e_in_z);
void Secundary_Target_Events(vector<double> I_point, double Theta_In, vector<double> U_vector, double e, double Detec_Radius, double Detec_Size_x, double Detec_Size_y, double Detec_Dist_Target, double Detec_Angle, double &e_in, double &e_out, double &Theta_Out, double &Scatter_Angle, double &e_out_z, vector<double> &V_vector);
void Join_Biparametrics(vector<string> files);
void deleteDirectoryContents(const std::string& dir_path);
double e_broh(double a, double q, double brho);
double PDF_Normal(double Energy, double sigEnergy);
double PDF_RealNumber();
double read_stopx();
double read_kineq();
double read_strag();
double GetSmaller(vector<double> Values);
double GetLarger(vector<double> Values);
double MoreRepeats(vector<double> Values, double smaller, double larger, double bin);
double avarage(vector<double> vec);
double standard_deviation(vector<double> strags);
double rad(double deg);
double deg(double rad);
double KinematicAngle_RecoilParticle(vector<double> V, vector<double> U, int Min, int Mout, double Ein, double Eout);

void secbeam(long nVector, double primaryTargThick, double massSecBeam, double qSecBeam, double Primary_Beam_Avarage_Energy, double Qreac, double ScatterAngMin, double ScatterAngMax, double FirstCollimatorRadius, double FaradayCupRadius, double LollipopRadius, double SecondCollimatorRadius, double SigmaPrimaryBeam, double SigmaPrimaryEnergy, double lolli2, double Rloli2, double currentFix, string primaryBeam, string primaryTarget, string secBeam, string CrossSection, string fixcurrent, string straggling, vector<string> Contaminant_Tag, vector<double> Contaminant_Mass, vector<double> Contaminant_Charge);
void bipa_contaminants(double Detec_Angle, double Detec_Radius, double Detec_Size_x, double Detec_Size_y, double Detec_Size_z, double Detec_Dist_Target, string Detec_Tag, double Sec_Targ_Thick);
void bipa(double Detec_Angle, double Detec_Radius, double Detec_Size_x, double Detec_Size_y, double Detec_Size_z, double Detec_Dist_Target, string Detec_Tag, double Sec_Targ_Thick);
#endif
