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

void bipa(double Detec_Angle, double Detec_Radius, double Detec_Size_x, double Detec_Size_y, double Detec_Size_z, double Detec_Dist_Target, string Detec_Tag, double Sec_Targ_Thick){
	cout << endl << "********************************************************************" << endl;
    cout <<         "***************    Biparametric graphs of energies   ***************" << endl;
    //Opening the output ofstreams
    remove("Biparametric/Outputs/ReactionPoints.dat");
    remove("Biparametric/Outputs/DetectionPoints.dat");
    ofstream bipametric("Biparametric/Outputs/BiParametric.dat");
    ofstream EjectEnergies("Biparametric/Outputs/EjectileEnergies.dat");
    ofstream EjectAngles("Biparametric/Outputs/EjectileScatterAngles.dat");
    ofstream Eject_Z_Axis_Angles("Biparametric/Outputs/EjectileAngles_Z_Axis.dat");
    ofstream RecoilEnergies("Biparametric/Outputs/RecoilEnergies.dat");
    ofstream RecoilAngles("Biparametric/Outputs/RecoilScatterAngles.dat");
    ofstream Recoil_Z_Axis_Angles("Biparametric/Outputs/RecoilAngles_Z_Axis.dat");

    //Colecting the aspects of the detector
    Detec_Angle = rad(Detec_Angle);
    //Importing the incident (point, vector and theta_in) and Secundary beam energies
    vector<vector<double>> U_vectors, I_points;
    vector<double> Thetas_In, Sec_Beam_Energies;
    Read_Matrix(U_vectors,"Biparametric/Inputs/BeamDirectionOnTarget.dat");
    Read_Matrix(I_points,"Biparametric/Inputs/PointsOnTarget.dat");
    Read_Vector(Thetas_In, "Biparametric/Inputs/AnglesOnTarget.dat");
    Read_Vector(Sec_Beam_Energies,"Biparametric/Inputs/SecEnergySecTarget.dat");
    int Cont_times = 0;
    while(true){
        //Colecting the beams and target informations
        double Qwindow, Qgs;
        string Secundary_Beam, Secundary_Target, Ejectile_Beam, Recoil_Particle, Ejectile_Charge;
        cin >> Secundary_Beam;
        if(Secundary_Beam == "END"){
            if(Cont_times == 0){
                cout << endl << "       No reaction were entered" << endl;
            }
            break;
        }
        Cont_times += 1;
        cin >> Secundary_Target >> Ejectile_Beam >> Recoil_Particle >> Qwindow >> Qgs >> Ejectile_Charge;

        cout << "> Sec Target Reaction: " << Secundary_Target << "(" << Secundary_Beam << "," << Ejectile_Beam << ")" << Recoil_Particle << " and " << Ejectile_Beam << " detection on telescope." << endl;

        /*
        On this loop, the program simulate and compute the paths in
        the secundary target and the incident point on the detector, 
        thus obtaining the scattering angle.
        Furthermore, the program compute the losses of energy, in the
        two paths and in the reaction considering the straggling.
        */
        vector<double> Energy_Out_Sec_Target, Energy_Out_Recoil, Scatter_Angles, Scatter_Angles_Kinematics, Scatter_Angles_Recoil, Scatter_Angles_Recoil_Kinematics;
        cout << "--> Simulating reaction and detection point" << endl;
        if(Thetas_In.size() == 0){
            cerr << RED << "ERROR: " << RESET << "No particles reached the secondary target." << endl;
            exit(0);
        }
        for(int i = 0; i < Thetas_In.size(); ++i){
            Load_bar((1.0*i)/Thetas_In.size());
            double e_in, e_out, e_out_z, e_out_recoil, theta_out, scatter_angle, scatter_angle_recoil, scatter_angle_recoil_kinematics, E_out_dx1, E_out_kineq_Sec_Beam, E_out_dx2, E_Out_kineq_Recoil, E_out_dx2_Recoil, Q, strag;
            vector<double> V_vector;
            Secundary_Target_Events(I_points[i], rad(Thetas_In[i]), U_vectors[i], (-1)*pow(10,-2)*Sec_Targ_Thick, Detec_Radius, Detec_Size_x, Detec_Size_y, Detec_Dist_Target, Detec_Angle, e_in, e_out, theta_out, scatter_angle, e_out_z, V_vector);
            e_in = (-1)*pow(10,2)*e_in;
            e_out = (-1)*pow(10,2)*e_out;
            Q = Qgs - Qwindow*PDF_RealNumber();                

            //Stopping power to primary particle
            Energy_Loss_stopx(e_in, Sec_Beam_Energies[i],Secundary_Beam,Secundary_Target, "yes",E_out_dx1, strag);
            //Energies for secundary beam and recoil particle
            Energy_Loss_Reaction_Bipa(E_out_dx1, deg(scatter_angle), Secundary_Beam, Secundary_Target, Ejectile_Beam, Recoil_Particle, Q, E_out_kineq_Sec_Beam, E_Out_kineq_Recoil, scatter_angle_recoil);
            //Stopping power to secundary paticle
            Energy_Loss_stopx(e_out, E_out_kineq_Sec_Beam, Ejectile_Beam, Secundary_Target, "yes", E_out_dx2, strag);
            

            //Calculating the kinematic scattering angle of recoling particle
            int Min, Mout;
            Particle_Mass(Secundary_Beam, Min);
            Particle_Mass(Ejectile_Beam, Mout);
            scatter_angle_recoil_kinematics = KinematicAngle_RecoilParticle(U_vectors[i],V_vector,Min,Mout,E_out_dx1,E_out_kineq_Sec_Beam);

            //Stopping power to recoil particle
            e_out_z = (-1)*pow(10,2)*e_out_z;
            e_out_recoil = (e_out_z)/cos(scatter_angle_recoil_kinematics);
            Energy_Loss_stopx(e_out_recoil, E_Out_kineq_Recoil, Recoil_Particle, Secundary_Target, "yes", E_out_dx2_Recoil, strag);

            //Saving the informations
            Energy_Out_Sec_Target.push_back(E_out_dx2);
            Energy_Out_Recoil.push_back(E_out_dx2_Recoil);
            Scatter_Angles.push_back(deg(scatter_angle));
            Scatter_Angles_Kinematics.push_back(deg(theta_out));
            Scatter_Angles_Recoil.push_back(scatter_angle_recoil);
            Scatter_Angles_Recoil_Kinematics.push_back(deg(scatter_angle_recoil_kinematics));
        }


        /*
        On this loop, the program compute the energy
        loss in the deltaE part of the detector, considering 
        the straggling.
        */
        vector<double> Delta_E, E;
        cout << endl << "--> Calculanting energy losses on detector" << endl;
        for(int i = 0; i < Energy_Out_Sec_Target.size(); ++i){
            Load_bar((1.0*i)/Thetas_In.size());
            double E_out_Detec, strag;
            Energy_Loss_stopx(Detec_Size_z, Energy_Out_Sec_Target[i], Ejectile_Beam, Detec_Tag, "yes", E_out_Detec, strag);
            E.push_back(E_out_Detec);
            Delta_E.push_back(Energy_Out_Sec_Target[i]-E_out_Detec);
        }
        cout << "\n\n";

        //Biparametric output
        for(int i = 0; i < E.size(); ++i){
            bipametric << fixed << setprecision(5) << Energy_Out_Sec_Target[i] << "   " << Delta_E[i]  << endl;
        }

        //Sec beam output
        for(int i = 0; i < Energy_Out_Sec_Target.size(); ++i){
            EjectEnergies << Energy_Out_Sec_Target[i] << endl;
            EjectAngles << Scatter_Angles[i] << endl;
            Eject_Z_Axis_Angles << Scatter_Angles_Kinematics[i] << endl;
        }

        //Recoil output
        for(int i = 0; i < Energy_Out_Recoil.size(); ++i){
            RecoilEnergies << Energy_Out_Recoil[i] << endl;
            RecoilAngles << Scatter_Angles_Recoil[i] << endl;
            Recoil_Z_Axis_Angles << Scatter_Angles_Recoil_Kinematics[i] << endl;
        }

        //Cleaning the bin directory
        Clear_directory_Bipa();
    }
}

void bipa_contaminants(double Detec_Angle, double Detec_Radius, double Detec_Size_x, double Detec_Size_y, double Detec_Size_z, double Detec_Dist_Target, string Detec_Tag, double Sec_Targ_Thick){
	cout << endl << "********************************************************************" << endl;
    cout <<         "*****    Biparametric graphs of energies of the contaminants   *****" << endl << endl;

    //Colecting the aspects of the detector
    Detec_Angle = rad(Detec_Angle);
    //Importing the incident (point, vector and theta_in) and Secundary beam energies
    vector<vector<double>> U_vectors, I_points;
    vector<double> Thetas_In, Sec_Beam_Energies;
    Read_Matrix(U_vectors,"Biparametric/Inputs/BeamDirectionOnTarget.dat");
    Read_Matrix(I_points,"Biparametric/Inputs/PointsOnTarget.dat");
    Read_Vector(Thetas_In, "Biparametric/Inputs/AnglesOnTarget.dat");
    
    vector<string> file_names;
    deleteDirectoryContents("Biparametric/Outputs/Contaminants");
    int Cont_times = 0;
    while(true){
        beguin:
        //Colecting the beams and target informations
        double Qwindow, Qgs;
        string Secundary_Beam, Secundary_Target, Ejectile_Beam, Recoil_Particle, Ejectile_Charge;
        cin >> Secundary_Beam;
        if(Secundary_Beam == "END"){
            if(Cont_times == 0){
                cout << endl << "       No contaminants were entered" << endl;
            }
            break;
        }
        Cont_times += 1;
        cin >> Secundary_Target >> Ejectile_Beam >> Recoil_Particle >> Qwindow >> Qgs >> Ejectile_Charge;

        stringstream formatted_charge;
        formatted_charge << fixed << setprecision(0) << Ejectile_Charge;
        string file_name_output = "Biparametric/Outputs/Contaminants/" + Ejectile_Beam + "_Q" + formatted_charge.str() + "_Energy.dat";

        ofstream bipametric(file_name_output);
        file_names.push_back(file_name_output);

        string file_name_input = "Biparametric/Inputs/" + Ejectile_Beam + "_Q" + formatted_charge.str() + "_Energy.dat";
        Sec_Beam_Energies.clear();
        Read_Vector(Sec_Beam_Energies,file_name_input);

        if(Sec_Beam_Energies.size()==0){
            cout << YELLOW << "Warning: " << RESET << "Please, check whether the beam " << Ejectile_Beam << " was included in the entry card!\n" << endl;
            goto beguin;
        }

        cout << "> Sec Target Reaction: " << Secundary_Target << "(" << Secundary_Beam << "," << Ejectile_Beam << ")" << Recoil_Particle << " and " << Ejectile_Beam << " detection on telescope." << endl;
        /*
        On this loop, the program simulate and compute the paths in
        the secundary target and the incident point on the detector, 
        thus obtaining the scattering angle.
        Furthermore, the program compute the losses of energy, in the
        two paths and in the reaction considering the straggling.
        */
        vector<double> Energy_Out_Sec_Target, Energy_Out_Recoil, Scatter_Angles, Scatter_Angles_Kinematics, Scatter_Angles_Recoil, Scatter_Angles_Recoil_Kinematics;
        cout << "--> Simulating reaction and detection point" << endl;
        if(Thetas_In.size() == 0){
            cerr << RED << "ERROR: " << RESET << "No particles reached the secondary target." << endl;
            exit(0);
        }
        for(int i = 0; i < Thetas_In.size(); ++i){
            //Just the load bar
            Load_bar((1.0*i)/Thetas_In.size());


            double thick_in, thick_out, thick_out_z, thick_out_recoil, theta_out, scatter_angle, scatter_angle_recoil, scatter_angle_recoil_kinematics, E_out_dx1, E_out_kineq_Sec_Beam, E_out_dx2, E_Out_kineq_Recoil, E_out_dx2_Recoil, Q, strag;
            vector<double> V_vector;
            Secundary_Target_Events(I_points[i], rad(Thetas_In[i]), U_vectors[i], (-1)*pow(10,-2)*Sec_Targ_Thick, Detec_Radius, Detec_Size_x, Detec_Size_y, Detec_Dist_Target, Detec_Angle, thick_in, thick_out, theta_out, scatter_angle, thick_out_z, V_vector);
            thick_in = (-1)*pow(10,2)*thick_in;
            thick_out = (-1)*pow(10,2)*thick_out;
            Q = Qgs - Qwindow*PDF_RealNumber();                

            //Stopping power to primary particle
            Energy_Loss_stopx(thick_in, Sec_Beam_Energies[i],Secundary_Beam,Secundary_Target, "yes",E_out_dx1, strag);
            
            //Energies for secundary beam and recoil particle
            Energy_Loss_Reaction_Bipa(E_out_dx1, deg(scatter_angle), Secundary_Beam, Secundary_Target, Ejectile_Beam, Recoil_Particle, Q, E_out_kineq_Sec_Beam, E_Out_kineq_Recoil, scatter_angle_recoil);

            //Stopping power to secundary paticle
            Energy_Loss_stopx(thick_out, E_out_kineq_Sec_Beam, Ejectile_Beam, Secundary_Target, "yes", E_out_dx2, strag);

            //Calculating the kinematic scattering angle of recoling particle
            int Min, Mout;
            Particle_Mass(Secundary_Beam, Min);
            Particle_Mass(Ejectile_Beam, Mout);
            scatter_angle_recoil_kinematics = KinematicAngle_RecoilParticle(U_vectors[i],V_vector,Min,Mout,E_out_dx1,E_out_kineq_Sec_Beam);

            //Stopping power to recoil particle
            thick_out_z = (-1)*pow(10,2)*thick_out_z;
            thick_out_recoil = (thick_out_z)/cos(scatter_angle_recoil_kinematics);
            Energy_Loss_stopx(thick_out_recoil, E_Out_kineq_Recoil, Recoil_Particle, Secundary_Target, "yes", E_out_dx2_Recoil, strag);

            //Saving the informations
            Energy_Out_Sec_Target.push_back(E_out_dx2);
            Energy_Out_Recoil.push_back(E_out_dx2_Recoil);
            Scatter_Angles.push_back(deg(scatter_angle));
            Scatter_Angles_Kinematics.push_back(deg(theta_out));
            Scatter_Angles_Recoil.push_back(scatter_angle_recoil);
            Scatter_Angles_Recoil_Kinematics.push_back(deg(scatter_angle_recoil_kinematics));
        }


        /*
        On this loop, the program compute the energy
        loss in the deltaE part of the detector, considering 
        the straggling.
        */
        vector<double> Delta_E, E;
        cout << endl << "--> Calculanting energy losses on detector" << endl;
        for(int i = 0; i < Energy_Out_Sec_Target.size(); ++i){
            Load_bar((1.0*i)/Thetas_In.size());
            double E_out_Detec, strag;
            Energy_Loss_stopx(Detec_Size_z, Energy_Out_Sec_Target[i], Ejectile_Beam, Detec_Tag, "yes", E_out_Detec, strag);
            E.push_back(E_out_Detec);
            Delta_E.push_back(Energy_Out_Sec_Target[i]-E_out_Detec);
        }
        cout << "\n\n";

        //Biparametric output
        for(int i = 0; i < E.size(); ++i){
            bipametric << fixed << setprecision(5) << Energy_Out_Sec_Target[i] << "   " << Delta_E[i]  << endl;
        }

        Clear_directory_Bipa();
    }
    file_names.push_back("Biparametric/Outputs/BiParametric.dat");
    Join_Biparametrics(file_names);
}
