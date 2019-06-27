
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//																														//
// C++ implementation of the metapopulation model of mosquito population dynamics described in the article:				//
//																														//
// "Fine-scale modelling finds that breeding site fragmentation can reduce mosquito population persistence"				//
//  (Communications Biology, 2019)																						//				
//																														//
//																														//
//	Dr. Clare McCormack																									//
//	c.mccormack14@imperial.ac.uk																						//
//																														//	
//																														//
//																														//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<cmath>
#include<random>
#include<algorithm>
#include<time.h>
#include<fstream>
#include<string>
#include<sstream>
#include<direct.h>
#include<windows.h>
#include"randlib_par.h"
#include"Dispersal_Functions.h"
#include"Functions.h"
#include"Make_Landscape.h"

using namespace std;

/*Declare arrays used in code*/

	int***V;
	int***L;
	int***L_final;
	int***E;
	int**total_V;
	int**total_L;
	int**total_E;
	int**next_l;
	int**next_e;
	int**next_v;
	int**pop_v;
	int*disp_fun_v;
	int**neighbours_total;
	int**ext_sites_v;
	int***tau_v;
	int*max;
	int*extinctions_v;
	int**no_dist;
	int***seed_v;
	int***seed_l;
	int**seed_patch;
	int**cumulative_pop;
	int***time_occupancy;
	int*time_occupancy_all;
	int***rescue_v;
	int**count_extinctions;
	int**count_rescues;
	int***cumulative_deaths;

	double***K;
	double**X;
	double**avg_total_V;
	double**avg_total_L;
	double**avg_total_E;
	double**site_avg_v;
	double**site_sd_v;
	double**site_avg_l;
	double**site_sd_l;
	double**avg_rescues;
	double*run_avg_v;
	double*run_sd_v;
	double*run_avg_l;
	double*run_sd_l;
	double*avg_ext_sites_v;
	double*prop_extinct;
	double*ext_time_v;
	double**avg_no_occupied;
	double**avg_prop_occupied;
	double**avg_prob_established;
	double**no_occupied;
	double**prop_occupied;
	double**prob_established;
	double**count_occupancy;
	double*avg_count_occupancy;
	double*avg_lifetime_disp;
	double*lifetime_disp;

	int**offset;/*offset matrix describing location of possible neighbours*/
	double****neighbours;/*array of all neighbours of each patch*/
	int***no_neighbours;/*no. of neighbours a patch has at each distance*/
	double****dispersed;/*no. of adult mosquitoes dispersing to neighbouring patches at each distance*/

int main(int argc, char*argv[])

{
	/*Model Parameters*/

		int runs; /*no. of model runs*/
		int sites; /*total number of patches on grid*/
		int n_row; /*no of rows*/
		int n_col; /*no. of columns*/
		int dispersal; /*type of dispersal: 0-none, 1-with dispersal*/
		int capacity; /*type of carrying capacity (K): 1-homogeneous landscape 2-heterogeneous landscape with coefficient of variation 'sigma', and level of correlation 'corr' */
		int max_dist; /*max range of dispersal*/
		int kernel; /*type of dispersal kernel: 0-neg exp, 1-normal, 2- exp power with shape param=4, 3-inverse power law, 4- 2D t distribution */
		double disp_rate;/*rate of dispersal*/
		double mean_dist;/*mean distance travelled for dispersing mosquitoes*/
		double b; /*avg no. of eggs laid per day*/
		double mu_e; /*death rate of eggs*/
		double mu_l; /*death rate of larvae*/
		double mu_v; /*death rate of adult mosquitoes*/
		double gamma_e; /*rate of development from egg to larva*/
		double gamma_l; /*rate of development from larva to adult*/
		double Rm; /*mosquito reproduction number*/
		double sigma; /*co-efficient of variation for heterogeneous landscapes*/
		double corr; /*spatial extent of correlation (corresponds to 'r' in (Hancock,2012))*/
		double omega; /*power on density dependence*/
		double cap;/*specified carrying capacity (in terms of number of larvae)*/
		double amp; /*amplitude of seasonal forcing*/
		double phase; /*phase of seasonal forcing*/
		int init_l; /*initial larval population in each breeding site*/
		int init_v; /*initial adult mosquito population in each breeding site*/
		int seed_mos; /*adult population seed size*/
		int gon_length; /*length of gonotrophic cycle (days)*/
		int divisions; /*no. of levels of graularity we want to look at*/
		int m;/*factor by which no. of patches is reduced as granulariy is reduced */
		
	/*manually input parameter values*/

		runs = 100; n_row = 1; n_col = 1; capacity = 1; sigma = 0, corr = 0; amp = 0; phase = 0; omega = 1, cap = 0,
		dispersal = 0; kernel = 0; max_dist=12; disp_rate = 0.08, mean_dist = 5;
		Rm = 2.69; gon_length = 3; gamma_e = .25; gamma_l = 0.07;  mu_v = 0.1; mu_l = 0.025; mu_e = 0.01;
		init_v = 1; init_l = 1000; seed_mos = 0; m = 4, divisions = 1;

	/*read in parameter values from file*/

		//ifstream Params;

		//string ParamFileName = argv[1];
		//Params.open(ParamFileName);

		//string param_name, param_value_string;

		//while (!Params.eof())
		//{
		//	getline(Params, param_name, '\t');
		//	getline(Params, param_value_string, '\n');;

		//	if (param_name == "runs")				runs = stoi(param_value_string);
		//	if (param_name == "n_row")				n_row = stoi(param_value_string);
		//	if (param_name == "n_col")				n_col = stoi(param_value_string);
		//	if (param_name == "capacity")			capacity = stoi(param_value_string);
		//	if (param_name == "sigma")				sigma = stod(param_value_string);
		//	if (param_name == "corr")				corr = stod(param_value_string);
		//	if (param_name == "amp")				amp = stod(param_value_string);
		//	if (param_name == "phase")				phase = stod(param_value_string);
		//	if (param_name == "omega")				omega = stod(param_value_string);
		//	if (param_name == "cap")				cap = stod(param_value_string);
		//	if (param_name == "dispersal")			dispersal = stoi(param_value_string);
		//	if (param_name == "kernel")				kernel = stoi(param_value_string);
		//	if (param_name == "max_dist")			max_dist = stod(param_value_string);
		//	if (param_name == "disp_rate")			disp_rate = stod(param_value_string);
		//	if (param_name == "mean_dist")			mean_dist = stod(param_value_string);
		//	if (param_name == "Rm")					Rm = stod(param_value_string);
		//	if (param_name == "gon_length")			gon_length = stoi(param_value_string);
		//	if (param_name == "gamma_e")			gamma_e = stod(param_value_string);
		//	if (param_name == "gamma_l")			gamma_l = stod(param_value_string);
		//	if (param_name == "mu_e")				mu_e = stod(param_value_string);
		//	if (param_name == "mu_l")				mu_l = stod(param_value_string);
		//	if (param_name == "mu_v")				mu_v = stod(param_value_string);
		//	if (param_name == "seed_mos")			seed_mos = stoi(param_value_string);
		//	if (param_name == "init_l")				init_l = stoi(param_value_string);
		//	if (param_name == "init_v")				init_v = stoi(param_value_string);
		//	if (param_name == "m")					m = stoi(param_value_string);
		//	if (param_name == "divisions")			divisions = stoi(param_value_string);
		//}

		//Params.close();

	/*length of model runs*/

			int days = 365 * 7; 
			double dt = .125; /*size of time step*/
			int steps = (int)(((double)days) / dt) + 1;

	/*declare variables*/

			int current_e, current_l, current_v,ovi_v, n_ovi,
			new_e, new_l, new_v,
			deaths_e, deaths_l, deaths_v,
			out_e, out_l,
			dispersed_v;

			double equil_e, equil_l, equil_v, K_t, p_e, p_l, p_v, p_ovi, h_v, h_e, h_l;

			double pi = 4 * atan(1);
		
			init_l = init_l * n_row*n_col;/*total initial larval population across all sites*/

	/*seed random number generators*/

			initSeeds(unsigned(time(NULL)), unsigned(time(NULL)));

			srand(unsigned(time(NULL)));

			random_device rd;
			mt19937 gen(rd());

	///////////////
	/*Model Setup*/
	///////////////

	/*allocate memory*/

			site_avg_v = new double*[days];
			for (int t = 0; t < days; t++) site_avg_v[t] = new double[runs];/*average adult mosquito population in a patch at day t in run r*/

			site_sd_v = new double*[days];
			for (int t = 0; t < days; t++) site_sd_v[t] = new double[runs];/*sd of average adult mosquito population in a patch at dayt in run r*/

			site_avg_l = new double*[days];
			for (int t = 0; t < days; t++) site_avg_l[t] = new double[runs];/*average larval population in a patch at day t in run r*/

			site_sd_l = new double*[days];
			for (int t = 0; t < days; t++) site_sd_l[t] = new double[runs];/*sd of average larval population in a patch at day t in run r*/

			run_avg_v = new double[days];/*avg adult mosquito population of a patch at day t, averaged over all runs*/
			run_sd_v = new double[days];/*sd of avg adult mosquito population of a patch at day t, averaged over all runs*/

			run_avg_l = new double[days];/*avg larval population of a patch at day t, averaged over all runs*/
			run_sd_l = new double[days];/*sd of avg larval population of a patchat day t, averaged over all runs*/

			total_L = new int*[days];
			for (int t = 0; t < days; t++) total_L[t] = new int[runs];/*average total larval population at day t in run r*/

			total_V = new int*[days];
			for (int t = 0; t < days; t++) total_V[t] = new int[runs];/*average total adult mosquito population at day t in run r*/

			total_E = new int*[days];
			for (int t = 0; t < days; t++) total_E[t] = new int[runs];/*average total egg population at day t in run r*/

			avg_total_E = new double*[days];/*total egg population at day t, averaged over all runs*/
			for (int t = 0; t < days; t++) { avg_total_E[t] = new double[2]; }

			avg_total_L = new double*[days]; /*total larval population at day t, averaged over all runs*/
			for (int t = 0; t < days; t++) { avg_total_L[t] = new double[2]; }

			avg_total_V = new double*[days]; /*total adult mosquito population at dayt, averaged over all runs*/
			for (int t = 0; t < days; t++) { avg_total_V[t] = new double[2]; }

			ext_sites_v = new int*[days];
			for (int t = 0; t < days; t++) ext_sites_v[t] = new int[runs];/*average number of patches whose adult mosquito population goes extinct at day t in run r*/

			avg_ext_sites_v = new double[days];/*proportion of patches whose adult mosquito population goes extinct at time t, averaged over all runs*/

			extinctions_v = new int[runs];/*no of patches whose adult population is extinct at the end of run r*/
			ext_time_v = new double[runs];/*time at which the total adult population goes extinct in run r*/

			prop_extinct = new double[days]; /*proportion of patches extinct at day t, averaged over all runs*/

			seed_patch = new int*[runs];
			for (int r = 0; r < runs; r++) { seed_patch[r] = new int[2]; } /*location of seeded patch for each run*/

			no_occupied = new double*[days];
			for (int t = 0; t < days; t++) { no_occupied[t] = new double[runs]; }/*no. of patches occupied with adult mosquitoes at day t, for each run*/

			prob_established = new double*[days];
			for (int t = 0; t < days; t++) { prob_established[t] = new double[runs]; }/*count if total adult mosquito pop>0 at day t, for each run*/

			prop_occupied = new double*[days];
			for (int t = 0; t < days; t++) { prop_occupied[t] = new double[runs]; }/*proportion of patches occupied with adult mosquitoes at day t, for each run*/

			avg_no_occupied = new double*[days];/*no. of patches occupied at day t, averaged over all runs*/
			for (int t = 0; t < days; t++) { avg_no_occupied[t] = new double[2]; }

			avg_prop_occupied = new double*[days]; /*proportion of patches occupied at day t, averaged over all runs*/
			for (int t = 0; t < days; t++) { avg_prop_occupied[t] = new double[2]; }
	
			avg_prob_established = new double*[days]; /*proportion of runs with total adult pop >0 at day t, averaged over all runs*/
			for (int t = 0; t < days; t++) { avg_prob_established[t] = new double[2]; }

			count_occupancy = new double*[days];
			for (int t = 0; t < days; t++) { count_occupancy[t] = new double[runs]; } /*proportion of patches occupied by day t, for each run*/

			avg_count_occupancy = new double[days]; /*proportion of patches occupied by day t, averaged over all runs*/

			time_occupancy_all = new int[runs]; /*time by which all patches have been occupied at least once, for each run*/

			disp_fun_v = new int[2]; /*stores no. of adult mosquitoes dipsersing and no of adult mosquito deaths in a given timestep*/

			lifetime_disp = new double[runs]; /*lifetime dispersal distance for adult mosquitoes for each run*/

			avg_lifetime_disp = new double[2]; /*average lifetime dispersal distance for adult mosquitoes across all runs*/

			cumulative_deaths = new int**[n_row];
			for (int i = 0; i < n_row; i++) {cumulative_deaths[i] = new int*[n_col];
			for (int j = 0; j < n_col; j++) cumulative_deaths[i][j] = new int[runs];} /*cumulative adult mosquito deaths in each patch, for each run*/
		
	/*run model for each possible level of granularity*/

		int start_disp = dispersal; /*if start_disp=0, run model for landscape with and without dipsersal, else run for landscape with dispersal only*/
		int end_disp = 1;

		for (int alpha = 1; alpha < (divisions + 1); alpha++) {

			/*first keep track of carrying capacity and initial seeded pop of each patch at previous level of granularity*/

				double***old_K;
				old_K = new double**[runs];
				for (int r = 0; r < runs; r++) {
					old_K[r] = new double*[n_row];
					for (int i = 0; i < n_row; i++) { old_K[r][i] = new double[n_col]; }
				}

				int***old_seed_v;
				old_seed_v = new int**[runs];
				for (int r = 0; r < runs; r++) {
					old_seed_v[r] = new int*[n_row];
					for (int i = 0; i < n_row; i++) old_seed_v[r][i] = new int[n_col];
				}

				if (alpha > 1) {

					for (int r = 0; r < runs; r++) {

						for (int i = 0; i < n_row; i++) {

							for (int j = 0; j < n_col; j++) {

								old_K[r][i][j] = K[r][i][j];
								old_seed_v[r][i][j] = seed_v[r][i][j];

							}

							delete[] K[r][i];
							delete[] seed_v[r][i];

						}

						delete[] K[r];
						delete[] seed_v[r];

					}

					delete[] K;
					delete[] seed_v;
				}

				int rows = n_row;
				int cols = n_row;

			/*adjust mean and max distance travelled accordingly as granularity is reduced*/

				if (alpha > 1) {

					sites = (rows*cols) / m;

					max_dist = int(max_dist / double(2));
					mean_dist = mean_dist / double(2);
				}

				else sites = rows*cols;

				n_row = int(sqrt(sites));
				n_col = int(sqrt(sites));

			/*allocate memory for new grid size*/

				K = new double**[runs];
				for (int r = 0; r < runs; r++) {K[r] = new double*[n_row];
				for (int i = 0; i < n_row; i++) { K[r][i] = new double[n_col]; } /*carrying capacity of each patch for run r*/}

				V = new int**[steps];
				for (int t = 0; t < steps; t++) {V[t] = new int *[n_row];
				for (int i = 0; i < n_row; i++) V[t][i] = new int[n_col];} /*total adult mosquito population in patch (i,j) at time t*/

				L = new int**[steps];
				for (int t = 0; t < steps; t++) {L[t] = new int *[n_row];
				for (int i = 0; i < n_row; i++) L[t][i] = new int[n_col];} /*larval population in patch (i,j) at time t*/

				L_final = new int**[n_row];
				for (int i = 0; i < n_row; i++) {L_final[i] = new int *[n_col];
				for (int j = 0; j < n_col; j++) L_final[i][j] = new int[runs];}/*larval population in patch (i,j) at the final time step in run r*/

				E = new int**[steps];
				for (int t = 0; t < steps; t++) {E[t] = new int *[n_row];
				for (int i = 0; i < n_row; i++) E[t][i] = new int[n_col];} /*no. of eggs in patch (i,j) at time t*/

				pop_v = new int *[n_row];
				for (int i = 0; i < n_row; i++) pop_v[i] = new int[n_col]; /*stores adult mosquito population in patch (i,j) in a given timestep*/

				next_l = new int*[n_row];
				for (int i = 0; i < n_row; i++) next_l[i] = new int[n_col]; /*larval population in patch (i,j) at next time-step*/

				next_e = new int*[n_row];
				for (int i = 0; i < n_row; i++) next_e[i] = new int[n_col]; /*egg population in patch (i,j) at next time-step*/

				next_v = new int*[n_row];
				for (int i = 0; i < n_row; i++) next_v[i] = new int[n_col]; /*adult mosquito population in patch (i,j) at next time-step*/

				tau_v = new int**[n_row];
				for (int i = 0; i < n_row; i++) {tau_v[i] = new int*[n_col];
				for (int j = 0; j < n_col; j++) tau_v[i][j] = new int[runs];}/*time at which adult mosquito population in patch (i,j) first goes extinct in run r*/

				rescue_v = new int**[n_row];
				for (int i = 0; i < n_row; i++) {rescue_v[i] = new int*[n_col];
				for (int j = 0; j < n_col; j++) rescue_v[i][j] = new int[runs];}/*time at which adult mosquito population in patch (i,j) is first rescued in run r*/

				count_extinctions = new int*[n_row];
				for (int i = 0; i < n_row; i++) count_extinctions[i] = new int[n_col];/* no. of times patch (i,j) goes extinct (mosquito pop) in a single run*/

				count_rescues = new int*[n_row];
				for (int i = 0; i < n_row; i++) count_rescues[i] = new int[n_col];/* no. of times patch (i,j) is rescued (mosquito pop) in a single run*/

				avg_rescues = new double*[n_row];
				for (int i = 0; i < n_row; i++) avg_rescues[i] = new double[n_col];/* average no. of time patch (i,j) is rescued (mosquito pop)*/

				neighbours_total = new int*[n_row];
				for (int i = 0; i < n_row; i++) neighbours_total[i] = new int[n_col];/*total no. of neighbours patch (i,j) has in given dispersal range*/

				seed_v = new int**[runs]; 
				{for (int r = 0; r < runs; r++) {seed_v[r] = new int*[n_row];
				for (int i = 0; i < n_row; i++) { seed_v[r][i] = new int[n_col]; } /*seeded adult mosquito pop in each patch*/}}

				time_occupancy = new int**[n_row];
				for (int i = 0; i < n_row; i++) {
					time_occupancy[i] = new int*[n_col];
					for (int j = 0; j < n_col; j++) time_occupancy[i][j] = new int[runs];
				}/*time at which patch (i,j) becomes occupied with mosquitoes (for the first time) in run r*/
		
				cumulative_pop = new int*[n_row];
				for (int i = 0; i < n_row; i++) cumulative_pop[i] = new int[n_col]; /*cumulative susceptible adult mosquito pop in patch (i,j) at given timestep*/

	///////////////////
	////*RUN MODEL*////
    ///////////////////

		/*create directory for output files*/

			_mkdir("Results");

		/*run model with/without dipsersal for exact same landscape*/

			for (dispersal = start_disp; dispersal < (end_disp + 1); dispersal += 1) {

				/*string for unique filenames */

					stringstream ss_m; 

					ss_m << runs << ',' << sites << ',' << n_col << ',' << n_row << ',' << dispersal << ',' << kernel << ',' << max_dist << ',' << disp_rate << ',' << mean_dist << ','
					<< capacity << ',' << sigma << ',' << corr << ',' << amp << ',' << phase << ',' << omega << ',' << cap << ','
					<< Rm << ',' << gon_length << ',' << gamma_e << ',' << gamma_l << ',' << mu_l << ',' << mu_v << ',' << mu_e << ',' << init_l << ',' << init_v << ',' << seed_mos << ',' << m << ',' << divisions;

				/*initialise variables*/

					double avg_extinctions_v = 0; /*average proportion of sites which are extinct (zero adult pop) at end of run*/

					double avg_tau_v = 0; /*average time step at which an adult population first goes extinct*/

					double avg_length_v = 0; /*average no. of steps until an adult pop recovers from first extinction*/

					double avg_rescued_v = 0; /*average proportion of sites which go extinct and are subsequently rescued*/

					double avg_recovery_v = 0; /*average proportion of sites which go extinct and are recovered by end of time period*/

					double avg_zero_v = 0;/*average no. of steps between extinction and recovery*/
			
					double avg_ext_time_v = 0;/*average time until the total adult population goes extinct*/

					double avg_time_occupancy = 0; /*average time until an individual patch becomes occupied*/

					double avg_time_occupancy_all = 0;/*average time until all patches becomes occupied*/

					int failures_v = 0; /*no. of runs where total adult pop goes extinct*/

					int failures_l = 0; /*no. of runs where total larval pop goes extinct*/

					int runs_ext = 0; /*no. of runs in which at least one site goes extinct*/

					int runs_rescue = 0; /*no. of runs in which at least one site is rescued*/

					int runs_persist_v = 0;  /*no. of runs in which total final adult pop is greater than 0*/

					int runs_persist_l = 0;/*no. of runs in which total final larval pop is greater than 0*/

					int zero_v = 0; /*sum of number of time steps between extinction and recovery for each run*/

					int runs_occupancy = 0; /*number of runs where all patches become occupied*/

					double equil_time = 0; /*time at which adult population reaches equilibrium*/

			/*initialise arrays*/

					for (int t = 0; t < days; t++) {

						run_avg_v[t] = 0; run_sd_v[t] = 0; run_avg_l[t] = 0; run_sd_l[t] = 0;

						avg_ext_sites_v[t] = 0; prop_extinct[t] = 0; avg_count_occupancy[t] = 0;

						for (int i = 0; i < 2; i++) {

							avg_total_V[t][i] = 0;
							avg_total_L[t][i] = 0;  avg_total_E[t][i] = 0;

							avg_no_occupied[t][i] = 0; avg_prob_established[t][i] = 0; avg_prop_occupied[t][i] = 0;
						}
					}

					for (int t = 0; t < days; t++) { for (int r = 0; r < runs; r++) { ext_sites_v[t][r] = 0; site_avg_v[t][r] = 0; site_sd_v[t][r] = 0; site_avg_l[t][r] = 0; site_sd_l[t][r] = 0; total_V[t][r] = 0;  total_L[t][r] = 0;  total_E[t][r] = 0; count_occupancy[t][r] = 0; no_occupied[t][r] = 0; prop_occupied[t][r] = 0; prob_established[t][r] = 0; } }

					for (int i = 0; i < n_row; i++) { for (int j = 0; j < n_col; j++) { for (int r = 0; r < runs; r++) { L_final[i][j][r] = 0; tau_v[i][j][r] = 0;  time_occupancy[i][j][r] = 0; rescue_v[i][j][r] = 0; cumulative_deaths[i][j][r] = 0;} } }

				/*set equilibrium population of each patch*/

					if (dispersal == start_disp) {

							if (alpha == 1) {

								for (int i = 0; i < n_row; i++) { 

									for (int j = 0; j < n_col; j++) { 
									
									equil_l = int(init_l / double(sites)); 
									equil_v = equil_l*gamma_l / mu_v; 
									equil_e = b*gamma_l*equil_l / (mu_v*(mu_e + gamma_e));
								} 
								}

								equil_e = round(equil_e);
								equil_v = round(equil_v);
							}

							else {

								equil_l = equil_l*m;
								equil_v = equil_v*m;
								equil_e = equil_e*m;
							}
						}

				/*calculate oviposition rate from value of Rm*/

					b = Rm / ((gamma_e / (gamma_e + mu_e))*(gamma_l / (gamma_l + mu_l))*(1 / mu_v));

					if (b <= 0) { cout << "Warning:b<=0" << endl; system("pause"); }

				/*if we have dispersal, create matrix giving neighbours of each patch, distance to neighbours and probability of dispersing to those neighbours*/

					if (dispersal > 0) {

						int n_total, a;

						offset = Create_Offset_Matrix(max_dist);/*matrix of all possible locations of neighbours within dispersal range*/

						neighbours = new double***[n_row];
						dispersed = new double***[n_row];
						no_neighbours = new int**[n_row];
						no_dist = new int*[n_row];

						for (int i = 0; i < n_row; i++) {

							neighbours[i] = new double**[n_col];
							dispersed[i] = new double**[n_col];
							no_neighbours[i] = new int*[n_col];

							no_dist[i] = new int[n_col];

							for (int j = 0; j < n_col; j++) {

								n_total = Find_Number_of_Neighbours(i, j, n_row, n_col, max_dist);/*no. of neighbours of patch (i,j) within dispersal range*/

								neighbours_total[i][j] = n_total;

								neighbours[i][j] = new double*[n_total];

								for (int k = 0; k < n_total; k++) { neighbours[i][j][k] = new double[4]; }

								Create_Neighbours_Matrix(i, j, n_row, n_col, n_total, kernel, mean_dist, max_dist, neighbours, offset, disp_rate);/*fill in entries of neighbours matrix*/

								no_dist[i][j] = 1;

								for (int k = 1; k < n_total; k++) { if (neighbours[i][j][k][0]>neighbours[i][j][k - 1][0]) no_dist[i][j] += 1; }

								a = no_dist[i][j];

								dispersed[i][j] = new double*[a];
								no_neighbours[i][j] = new int[a];

								for (int d = 0; d < a; d++) { dispersed[i][j][d] = new double[3]; }
							}
						}

							/*find max flux out of a single patch*/

							double**flux;
							double max_f = 0;
							double i1, j1, p;

							flux = new double*[n_row];

							for (int i = 0; i < n_row; i++) { flux[i] = new double[n_col]; }

							for (int i = 0; i < n_row; i++) {

								for (int j = 0; j < n_col; j++)

								{
									flux[i][j] = 0;

									n_total = neighbours_total[i][j];

									for (int k = 0; k < n_total; k++) {

										flux[i][j] += neighbours[i][j][k][1];

									}

									if (flux[i][j]>max_f) max_f = flux[i][j];
								}
							}

							/*normalise values generated by kernel using max flux*/

							for (int i = 0; i < n_row; i++) {

								for (int j = 0; j < n_col; j++)

								{

									n_total = neighbours_total[i][j];

									for (int k = 0; k < n_total; k++) {

										p = neighbours[i][j][k][1];/*value generated by dispersal kernel for neighbour k of patch (i,j)*/
										i1 = neighbours[i][j][k][2];/*row number of the neighbouring patch*/
										j1 = neighbours[i][j][k][3];/*column number of the neighbouring patch*/

										if (max_f > 0) {

											if (i1 == i && j1 == j) {

												neighbours[i][j][k][1] = disp_rate*((1 - flux[i][j] / max_f) + p / max_f); /*take account of mosquitoes moving within the same patch*/
											}

											else
												neighbours[i][j][k][1] = disp_rate*(p / max_f);
										}

										else
											neighbours[i][j][k][1] = disp_rate;

									}
								}
							}

							//ofstream results_neighbours;
							//results_neighbours.open("Results/Neighbours " + ss_m.str() + ".csv");

							//for (int i = 0; i < n_row; i++) {

							//	for (int j = 0; j < n_col; j++)

							//	{
							//		n_total = neighbours_total[i][j];

							//		for (int k = 0; k < n_total; k++) {

							//			results_neighbours << i << ',' << j << ',' << neighbours[i][j][k][0] << ',' << neighbours[i][j][k][1] << ',' << neighbours[i][j][k][2] << ',' << neighbours[i][j][k][3] << endl;

							//		}
							//	}
							//}
							//
					}
	
		/*for each run*/

				for (int r = 0; r < runs; r++)

				{
					/*initialise variables*/
					int count = 0;
					extinctions_v[r] = 0;
					ext_time_v[r] = 0;
					time_occupancy_all[r] = 0;
					lifetime_disp[r] = 0;

					for (int i = 0; i < n_row; i++) { for (int j = 0; j < n_col; j++) { cumulative_pop[i][j] = 0; } }

					/*deterministic carrying capacity of each patch with initial larval population 'equil_l'*/
					
					double det_K = equil_l / pow(((gamma_l / mu_l)*((b*gamma_e) / (mu_v*(mu_e + gamma_e)) - 1) - 1), 1 / omega);

					if (det_K <= 0) { cout << "Warning:det_K<=0" << endl; system("pause"); }

					/*make landscape*/

					if (dispersal == start_disp) { Make_Landscape(capacity, K, old_K, r, n_row, n_col, alpha, det_K, sigma, corr); }

					ofstream results_capacity;
					results_capacity.open("Results/K " + ss_m.str() + ".csv");

					for (int i = 0; i < n_row; i++) {

						for (int j = 0; j < n_col; j++) {

							results_capacity << K[r][i][j] << ',';
						}
						
						results_capacity << endl;
					}

					/*set initial population of each patch*/

						if (seed_mos == 0) {

							/*same initial population in each patch*/

							for (int i = 0; i < n_row; i++) {

								for (int j = 0; j < n_col; j++) {

									V[0][i][j] = int(equil_v); L[0][i][j] = int(equil_l); E[0][i][j] = int(equil_e);

								}
							}
						}

						else {

							/*seed random patch with adult population*/

							if (alpha == 1) {

								/*keep seed the same in case of dispersal or no dispersal*/

								if (dispersal == start_disp) {

									seed_patch[r][0] = rand() % rows;
									seed_patch[r][1] = rand() % cols;

								}

								int row = seed_patch[r][0];
								int column = seed_patch[r][1];

								for (int i = 0; i < n_row; i++) {

									for (int j = 0; j < n_col; j++) {

										E[0][i][j] = 0;
										L[0][i][j] = 0;

										if (i == row && j == column) { seed_v[r][i][j] = seed_mos; }

										else
										{
											seed_v[r][i][j] = 0;
										}

										V[0][i][j] = seed_v[r][i][j];

									}
								}
							}

							/*for lower levels of resolution, sum value of inital pop at previous level of resolution to give inital pop values for current level of resolution*/

							else if (alpha > 1) {

								int x, y;
								int sum1;
								x = 0;
								y = 0;

								for (int i = 0; i < (2 * n_row); i += 2) {

									for (int j = 0; j < (2 * n_col); j += 2) {

										sum1 = 0;
										sum1 += old_seed_v[r][i][j] + old_seed_v[r][i][j + 1] + old_seed_v[r][i + 1][j] + old_seed_v[r][i + 1][j + 1];

										seed_v[r][x][y] = sum1;


										y = y + 1;
									}

									x = x + 1;
									y = 0;
								}

								for (int i = 0; i < n_row; i++) {

									for (int j = 0; j < n_col; j++) {

										E[0][i][j] = 0;
										L[0][i][j] = 0;
										V[0][i][j] = seed_v[r][i][j];

									}
								}

							}
						}

					/*for each time step*/
		
						for (int t = 0; t < steps - 1; t++)

						{
							
							int t1 = int(t / (1 / dt)); /*day number*/

							/*for each patch*/
				
							for (int i = 0; i < n_row; i++) {

								for (int j = 0; j < n_col; j++) {

									if (i == 0 && j == 0) { for (int x1 = 0; x1 < n_row; x1++) { for (int y1 = 0; y1 < n_col; y1++) { pop_v[x1][y1] = 0; } } }

									/*current population in (i,j) */
									current_v = V[t][i][j];
									current_l = L[t][i][j];
									current_e = E[t][i][j];

									/*Set thread number as code not multi-threaded*/
									int thread = 0;

									/*carrying capacity of patch (i,j) at time t*/
									if (cap == 0) { K_t = K[r][i][j] * (1 + amp*cos(2 * pi*(double(t) / (365 * (1 / dt)) + phase))); }/*carrying capacity of patch (i,j) at time t (with seasonal forcing)*/

									/*for single patch approximations-set equilibrium larval population to 'cap' (no seasonal variation)*/

									if (cap > 0) {K_t = cap / (pow(((gamma_l / mu_l)*((b*gamma_e) / (mu_v*(mu_e + gamma_e)) - 1) - 1), 1 / omega));}

									h_e = gamma_e + mu_e; /*total hazard of leaving egg compartment*/

									h_l = gamma_l + mu_l*(1 + pow(current_l / K_t, omega));/*total hazard of an individual larva leaving larval compartment*/

									h_v = mu_v; /*total hazard of leaving adult mosquito compartment*/

									p_e = 1 - exp(-h_e*dt);/*probability of leaving egg compartment*/
									p_l = 1 - exp(-h_l*dt);	/*probability of leaving larval compartment*/
									p_v = 1 - exp(-h_v*dt);/*probability of leaving adult mosquito compartment*/

									ovi_v = current_v; /*number of mosquitoes ovipositing at (i,j)*/
	
									/*determine number of adults mosquitoes leaving (i,j)*/	

									switch (dispersal)

									{

									case 0:

										deaths_v = ignbin_mt((long)current_v, p_v, thread); /*no. of adult mosquito deaths*/
										dispersed_v = 0; /*no. of adult mosquito dispersing from patch (i,j)*/
							
										break;

									case 1:

										Dispersal(i, j, n_row, n_col, runs, disp_fun_v, pop_v, current_v, neighbours, neighbours_total, dispersed, no_neighbours, no_dist, max_dist, mu_v, t, dt, disp_rate);/*no. of adult leaving site (i,j) through death or dispersal*/
										deaths_v = disp_fun_v[0];/*no. of adult mosquito deaths*/
										dispersed_v = disp_fun_v[1];/*no. of adult mosquito dispersing from patch (i,j)*/

										break;
									}

									p_ovi = 1 / double(gon_length); /*probability of laying eggs*/
									n_ovi = ignbin(ovi_v, p_ovi); /*number of adults laying eggs*/

									new_e = ignpoi_mt(b*gon_length*n_ovi*dt, thread);/*no. of eggs laid*/
									out_e = ignbin_mt((long)current_e, p_e, thread);/*total no. leaving egg compartment*/
									new_l = ignbin_mt((long)out_e, gamma_e / h_e, thread);/*no. of new larvae*/
									deaths_e = out_e - new_l; /*no. of egg deaths*/

									out_l = ignbin_mt((long)current_l, p_l, thread);/*total no. leaving larval compartment*/
									new_v = ignbin_mt((long)out_l, gamma_l / h_l, thread);/*no. of new adults among those leaving larval compartment*/
									deaths_l = out_l - new_v; /*no. of larval deaths*/

									next_e[i][j] = current_e + new_e - new_l - deaths_e;

									next_l[i][j] = current_l + new_l - new_v - deaths_l;

									next_v[i][j] = current_v + new_v - deaths_v - dispersed_v;

									/*population in (i,j) at time t+1, not accounting for adult mosquito dispersal*/

									E[t + 1][i][j] = next_e[i][j];
									L[t + 1][i][j] = next_l[i][j];
									V[t + 1][i][j] = next_v[i][j];
						
									//cumulative_deaths[i][j][r] += deaths_v;
								}
							}
					

							for (int i = 0; i < n_row; i++) {

								for (int j = 0; j < n_col; j++) {

									if (dispersal == 1) { V[t + 1][i][j] += pop_v[i][j]; } /*take account of adults moving to patch (i,j) at time t*/

									if (t % (int(1 / dt)) == 0) {

										int step = int(1 / dt);

										int prev = cumulative_pop[i][j];/*cumulative adult population in (i,j) at previous time step*/

										cumulative_pop[i][j] += V[t][i][j]; /*cumulative adult population in (i,j) at current time step*/

										/*total pop at time t1, accounting for adult mosquito dispersal*/

										total_V[t1][r] += V[t][i][j];
										total_L[t1][r] += L[t][i][j];
										total_E[t1][r] += E[t][i][j];

										/*count if the adult mosquito population goes extinct, and record time step at which this first occurs*/
										if (t1 > 0 && V[t - step][i][j] != 0 && V[t][i][j] == 0) { if (tau_v[i][j][r] == 0) { ext_sites_v[t1][r] += 1; tau_v[i][j][r] += (t1); } }

										/*count if the adult mosquito population is extinct at time t1*/
										if (V[t][i][j] == 0) { prop_extinct[t1] += 1; }

										/*record how many steps from initial extinction until the adult mosquito population recovers*/
										if (t1>0 && V[t - step][i][j] == 0 && V[t][i][j] != 0) { if (rescue_v[i][j][r] == 0) rescue_v[i][j][r] += (t1)-tau_v[i][j][r]; }

										/*at the final time step, count how many adult mosquito populations have gone extinct*/
										if (t1 == (days - 1) && V[t][i][j] == 0) { extinctions_v[r] += 1; }

										/*record larval population at final time step*/
										if (t1 == (days - 1)) { L_final[i][j][r] += L[t][i][j]; }

										/*count if patch (i,j) is occupied at time t1*/
										if (V[t][i][j] > 0) { no_occupied[t1][r] += 1; }

										/*count if patch (i,j) is newly occupied at time t1, and record time step at which this occurs*/
										if (prev == 0 && cumulative_pop[i][j] > 0) { time_occupancy[i][j][r] += t1; count += 1; }

										/*if all patches have been occupied at least once by time t1, record this time step*/
										if (time_occupancy_all[r] == 0 && count == sites) { time_occupancy_all[r] += t1; runs_occupancy += 1; }
									}
								}
							}

							if (t % (int(1 / dt)) == 0) {

								count_occupancy[t1][r] += count / double(sites); /*proportion of patches occupied at least once by time t1*/

								/*count if total pop is greater than zero at time t*/

								if (total_V[t1][r] > 0) {prob_established[t1][r] += 1;}
							}
				
						}

					/*count if total adult mosquito pop goes extinct in run, and record time at which this happens*/

						if (total_V[days - 1][r] == 0) {

								failures_v += 1;

								for (int t = 0; t < days - 2; t++) { if (total_V[t][r] != 0 && total_V[t + 1][r] == 0) { if (ext_time_v[r] == 0) ext_time_v[r] += (t + 1); } }
						}

						if (total_L[days - 1][r] == 0) {failures_l += 1;}

					/*For runs where larval population persists, calculate mean adult mosquito and larval population at each time step (averaged over all patches)*/

						if (total_L[days - 1][r] != 0) {

								for (int t = 0; t < days; t++)

								{
									site_avg_v[t][r] = double(total_V[t][r]) / sites;
									site_avg_l[t][r] = double(total_L[t][r]) / sites;

									run_avg_v[t] += site_avg_v[t][r];
									run_avg_l[t] += site_avg_l[t][r];

									avg_ext_sites_v[t] += double(ext_sites_v[t][r]) / (runs*sites);
								}
							}

					/*For runs where the adult population persists, calculate the mean time until a patch first becomes occupied*/

							int d, n1;

							if (total_V[days - 1][r] != 0) {

								d = 0;
								n1 = 0; /*no. of patches which become occupied (excluding seeded patch)*/

								for (int i = 0; i < n_row; i++) {

									for (int j = 0; j < n_col; j++) {

										d += time_occupancy[i][j][r];

										if (time_occupancy[i][j][r]> 0) n1 += 1;
									}
								}

								if (n1 > 0) { avg_time_occupancy += d / double(n1); }

								avg_time_occupancy_all += time_occupancy_all[r];
							}

					/*Calculate mean lifetime dispersal distance for run r*/

						//double dist;

						//for (int i = 0; i < n_row; i++) {

						//	for (int j = 0; j < n_col; j++) {

						//			for (int n = 0; n < neighbours_total[i][j];n++) {

						//				if (neighbours[i][j][n][2] == seed_patch[r][0] && neighbours[i][j][n][3] == seed_patch[r][1]) {

						//					dist = neighbours[i][j][n][0]; /*distance from seeded patch*/

						//					lifetime_disp[r] += (dist * cumulative_deaths[i][j][r]) / double(seed_mos);

						//					break;
						//				}

						//			}

						//		}
						//}

					/*Calculate mean proportion of runs for which total adult pop is non-zero at eact time step t*/ 

						 for (int t = 0; t < days; t++) { avg_prob_established[t][0] += prob_established[t][r] / double(runs); }

					/*Calculate average proportion of sites which are extinct at the end of a run*/

							avg_extinctions_v += extinctions_v[r] / (sites*double(runs));

					/*For each run, calculate mean time until adult pop in a patch first goes extinct, and mean no. of time steps until an adult pop is rescued*/

							n1 = 0; /*no. of patches which go extinct*/
							int n2 = 0; /*no. of patches which are rescued*/

							d = 0; int e = 0;
							double q = 0; double l = 0;

							for (int i = 0; i < n_row; i++) { for (int j = 0; j < n_col; j++) { e += (rescue_v[i][j][r]); d += tau_v[i][j][r]; if (tau_v[i][j][r] != 0) { n1 += 1; } if (rescue_v[i][j][r] != 0) { n2 += 1; } } }

							if (n1 > 0) {

								q = double(d) / n1;
								runs_ext += 1; /*count if at least one site goes extinct*/
							}

							else

								q = 0;

							if (n2 > 0)

							{
								l = double(e) / n2;
								runs_rescue += 1;/*count if at least one site is rescued*/
							}

							else

								l = 0;

							avg_tau_v += q;
							avg_length_v += l;

							/*For patches going extinct, calculate (1) average no of time steps between extinction and rescue, (2) the proportion which go extinct and are subsequently rescued and (3) the proportion which go extinct and are recovered at end of time period*/

							int t0 = 0; /*set t0>0 if we want to specify a burn in time before calculating summary stats*/

							if (r == 0) { for (int i = 0; i < n_row; i++) { for (int j = 0; j < n_col; j++) count_extinctions[i][j] = 0; } }

							for (int i = 0; i < n_row; i++) { for (int j = 0; j < n_col; j++) count_rescues[i][j] = 0; }

							int count_2;
							int sum_rescues = 0;
							zero_v = 0;

							for (int i = 0; i < n_row; i++) {

								for (int j = 0; j < n_col; j++) {

									for (int t = t0; t < steps - 2; t++) {

										if (V[t][i][j] != 0 && V[t + 1][i][j] == 0) count_extinctions[i][j] += 1;
										if (V[t][i][j] == 0 && V[t + 1][i][j] != 0) count_rescues[i][j] += 1; /*for each run, count the number of times the adult population is rescued*/
									}

								}
							}

							for (int i = 0; i < n_row; i++) { for (int j = 0; j < n_col; j++) { sum_rescues += count_rescues[i][j]; avg_rescues[i][j] += double(count_rescues[i][j]) / runs; } }

							for (int i = 0; i < n_row; i++) {

								for (int j = 0; j < n_col; j++) {

									count_2 = 0;

									if (count_rescues[i][j] > 0) {

										for (int t = t0; t < steps - 1; t++) {

											if (V[t][i][j] == 0 && rescue_v[i][j][r]>0) zero_v += 1;
											if (V[t][i][j] == 0 && V[t + 1][i][j] != 0) { count_2 += 1; }

											if (count_2 == count_rescues[i][j]) break; /*only want to consider cases where pop in rescued*/
										}
									}
								}
							}

							if (sum_rescues > 0) { avg_zero_v += zero_v / double(sum_rescues); }/*(1)*/

							for (int i = 0; i < n_row; i++) {

								for (int j = 0; j < n_col; j++) {

									if (tau_v[i][j][r]>0 && rescue_v[i][j][r]>0)   avg_rescued_v += 1 / double(n1);/*(2)*/
									if (tau_v[i][j][r] > 0 && V[steps - 2][i][j] > 0) avg_recovery_v += 1 / double(n1);/*(3)*/

								}
							}

							/*For runs where the adult population persists, calculate the mean proportion of sites occupied at each time step and mean proportion of sites which have been occupied at least once at each time step*/

							if (total_V[days - 1][r] != 0) {

								for (int t = 0; t < days; t++) {

									avg_no_occupied[t][0] += no_occupied[t][r];
									avg_prop_occupied[t][0] += no_occupied[t][r] / sites;
									avg_count_occupancy[t] += count_occupancy[t][r];

								}
							}

							cout << r << endl;

						}

	//////////////////////
	/*Summary Statistics*/
	//////////////////////

			runs_persist_v = runs - failures_v; /*no. of runs where adult population persists*/
			runs_persist_l = runs - failures_l; /*no. of runs where larval population persists*/

			/*calculate mean and variance of larval population across sites at final time step*/

				double var = 0;
				double mean = 0;
				double site_mean, site_var, diff, x, y;

				if (runs_persist_l > 0) {

					for (int i = 0; i < n_row; i++)

					{
						for (int j = 0; j < n_col; j++)

						{
							site_mean = 0;
							site_var = 0;
							y = 0;

							for (int r = 0; r < runs; r++)
							{

								if (total_L[days - 1][r] != 0) {

									site_mean += L_final[i][j][r] / double(runs_persist_l);
								}

							}

							for (int r = 0; r < runs; r++)

							{
								if (total_L[days - 1][r] != 0) {

									diff = L_final[i][j][r] - site_mean;
									x = pow(diff, 2);
									y += x;
								}
							}

							site_var = y / double(runs_persist_l);

							var += site_var / double(sites);
							mean += site_mean / sites;

						}
					}
				}

			/*calculate average proportion of patches which are extinct at each time step*/

				for (int t = 0; t < days; t++) { prop_extinct[t] = prop_extinct[t] / double(runs*sites); }

			/*calculate adult and larval population at time t, averaged over all sites and runs, and corresponding std dev*/

				for (int t = 0; t < days; t++) {

					run_avg_v[t] = run_avg_v[t] / (runs_persist_v);
					run_avg_l[t] = run_avg_l[t] / (runs_persist_l);

					run_sd_v[t] = Std_Dev_Doub(site_avg_v[t], total_V, days - 1, run_avg_v[t], runs, runs_persist_v);
					run_sd_l[t] = Std_Dev_Doub(site_avg_l[t], total_L, days - 1, run_avg_l[t], runs, runs_persist_l);

				}

			/*calculate total adult mosquito and larval population at time t, averaged over all runs where population persists, and corresponding std dev*/

				if (runs_persist_v > 0) {

					for (int t = 0; t < days; t++) {

						avg_total_E[t][0] = Mean_Value_Int(total_E[t], total_V, days - 1, runs, runs_persist_v);
						avg_total_L[t][0] = Mean_Value_Int(total_L[t], total_L, days - 1, runs, runs_persist_l);
						avg_total_V[t][0] = Mean_Value_Int(total_V[t], total_V, days - 1, runs, runs_persist_v);

						avg_total_E[t][1] = Std_Dev_Int(total_E[t], total_V, days - 1, avg_total_E[t][0], runs, runs_persist_v);
						avg_total_L[t][1] = Std_Dev_Int(total_L[t], total_V, days - 1, avg_total_L[t][0], runs, runs_persist_l);
						avg_total_V[t][1] = Std_Dev_Int(total_V[t], total_V, days - 1, avg_total_V[t][0], runs, runs_persist_v);

					}
				}

			/*calcute time at which adult mosquito population reaches equilibrium*/

				for (int t = 0; t < days; t++) {

					if ((fabs(avg_total_V[days - 1][0] - avg_total_V[t][0]) < .5) && equil_time == 0) { equil_time += t; }

				}

			/*Calculate the average time to extinction per run*/

				double k = 0;

				for (int r = 0; r < runs; r++)
				{
					if (ext_time_v[r] > 0) { k += ext_time_v[r]; }
				}

				if (failures_v > 0)

				{
					avg_ext_time_v = k / double(failures_v);
				}

				else

					avg_ext_time_v = 0;

				/*Calculate mean lifetime dispersal distance across runs*/

				//for (int r = 0; r < runs; r++) {

				//	avg_lifetime_disp[0] += lifetime_disp[r] / runs;
				//}

				//double sum = 0;

				//for (int r = 0; r < runs; r++) {

				//	diff = lifetime_disp[r] - avg_lifetime_disp[0];
				//	x = pow(diff, 2);
				//	sum += x;
				//}

				//avg_lifetime_disp[1] = sqrt(sum / double(runs));

			/*Calculate average time step at which each site first goes extinct*/

				ofstream results_km;
				results_km.open("Results/KM Sites " + ss_m.str() + ".csv");
				
				int s;

				for (int i = 0; i < n_row; i++) {

					for (int j = 0; j < n_col; j++) {

						double k = 0;
						double l = 0;
						double r1 = 0;

						if (i == 0) { s = j + 1; }

						else { s = i*n_row + j + 1; }

						for (int r = 0; r < runs; r++) {

							k += tau_v[i][j][r];

							if (tau_v[i][j][r]>0) r1 += 1;
						}

						if (r1 > 0) { l = k / r1; }

						else l = days / dt;

						results_km << s << ',' << l << endl;
					}
				}

			/*For runs where adult population persists, calculate average proportion of patches occupied and corresponding standard deviation*/

				if (runs_persist_v > 0) {

					for (int t = 0; t < days; t++) {

						avg_no_occupied[t][0] = avg_no_occupied[t][0] / runs_persist_v;
						avg_prop_occupied[t][0] = avg_prop_occupied[t][0] / runs_persist_v;
						avg_count_occupancy[t] = avg_count_occupancy[t] / runs_persist_v;

						for (int r = 0; r < runs; r++) {prop_occupied[t][r] = no_occupied[t][r] / sites;}
					}

					for (int t = 0; t < days; t++) {

						avg_no_occupied[t][1] = Std_Dev_Doub(no_occupied[t], total_V, days - 1, avg_no_occupied[t][0], runs, runs_persist_v);
						avg_prop_occupied[t][1] = Std_Dev_Doub(prop_occupied[t], total_V, days - 1, avg_prop_occupied[t][0], runs, runs_persist_v);
						avg_prob_established[t][1] = Std_Dev_Doub(prob_established[t], total_V, days - 1, avg_prob_established[t][0], runs, runs_persist_v);
					}

					avg_time_occupancy = avg_time_occupancy / runs_persist_v;
				}

			/*For runs where all patches become occupied, calculate mean time until all patches become occupied*/

				if (runs_occupancy>0) {avg_time_occupancy_all = avg_time_occupancy_all / runs_occupancy;}

			/*Calculate mean proportion of patches which are rescued and recovered, mean time to extinction of an individual patch*/

				if (runs_ext > 0) { avg_rescued_v = avg_rescued_v / runs_ext; avg_recovery_v = avg_recovery_v / runs_ext; avg_tau_v = avg_tau_v / runs_ext; }

			/*Calculate mean no. of days a patch is extinct for and mean time mean time between extinction and recovery */

				if (runs_rescue > 0) { avg_length_v = avg_length_v / runs_rescue; avg_zero_v = avg_zero_v / runs_rescue; }

	///////////////////////////
	/*Output Results to Files*/
	///////////////////////////

			/*stores adult mosquito population at each time step, averaged over all patches and all runs, and corresponding std dev between runs*/
			ofstream results;
			results.open("Results/Average Mosquito Population " + ss_m.str() + ".csv");

			/*stores larval population at each time step, averaged over all patches and all runs, and corresponding std dev between runs*/
			ofstream results_2;
			results_2.open("Results/Average Larval Populations " + ss_m.str() + ".csv");

			/*stores total adult mosquito population at time t (including all sub compartments), averaged over all runs, and corresponding std dev between runs*/
			ofstream results_3;
			results_3.open("Results/Total Adult Population " + ss_m.str() + ".csv");

			/*stores total larval population at time t, averaged over all runs, and corresponding std dev between runs*/
			ofstream results_4;
			results_4.open("Results/Total Larval Population " + ss_m.str() + ".csv");

			/*stores total egg population at time t, averaged over all runs, and corresponding std dev between runs*/
			ofstream results_5;
			results_5.open("Results/Total Egg Population" + ss_m.str() + ".csv");

			/*stores mean proportion of patches extinct at each time step*/
			ofstream results_6;
			results_6.open("Results/Extinctions " + ss_m.str() + ".csv");

			/*stores proportion of patches whose adult population first goes extinct at time t, averaged over all runs*/
			ofstream results_7;
			results_7.open("Results/Extinction Times " + ss_m.str() + ".csv");

			/*stores mean and variance of final larval pop across sites*/
			ofstream results_8;
			results_8.open("Results/Final Larval Pop " + ss_m.str() + ".csv");

			/*stores average number and proportion of patches occupied*/
			ofstream results_9;
			results_9.open("Results/Seeding " + ss_m.str() + ".csv");

			/*stores summary stats e.g. mean time to extinction*/
			ofstream summary;

			summary.open("Results/Summary " + ss_m.str() + ".csv");

			for (int t = 0; t < days; t++) {

				results << run_avg_v[t] << ',' << run_sd_v[t] << endl;
				results_2 << run_avg_l[t] << ',' << run_sd_l[t] << endl;
				results_3 << avg_total_V[t][0] << ',' << avg_total_V[t][1] << endl;
				results_4 << avg_total_L[t][0] << ',' << avg_total_L[t][1] << endl;
				results_5 << avg_total_E[t][0] << ',' << avg_total_E[t][1] << endl;
			
				results_7 << avg_ext_sites_v[t] << endl;
				results_6 << prop_extinct[t] << endl;
			}

			results_6 << endl;

			results_8<< mean << ',';
			results_8 << var << ',' << endl;

			for (int t = 0; t < days; t++) {

				results_9 << t << ',' << avg_no_occupied[t][0] << ',' << avg_no_occupied[t][1] << ',' << avg_prop_occupied[t][0] << ',' << avg_prop_occupied[t][1] << ',' << avg_prob_established[t][0] << ',' << avg_prob_established[t][1] << ',' << avg_count_occupancy[t] << endl;

			}

			summary << init_v << ',' << K[0][0][0] << ',' << avg_total_V[days - 1][0] << ',' << avg_extinctions_v << ',' << avg_ext_time_v << ',' << avg_tau_v << ',' << avg_length_v << ',' << avg_zero_v << ',' << avg_rescued_v << ',' << avg_recovery_v << ',' << avg_time_occupancy << ',' << avg_time_occupancy_all << ',' << runs_persist_v / double(runs) << ',' << runs_occupancy / double(runs) << ',' << equil_time << ',' << avg_lifetime_disp[0] << ',' << avg_lifetime_disp[1] << endl;

		}

	////////////////////////
	/*Detele Memory Stored*/
	////////////////////////

			if (dispersal>0) {

				for (int i = 0; i < n_row; i++) {

					for (int j = 0; j < n_col; j++) {

						int n2 = neighbours_total[i][j];

						for (int n1 = 0; n1 < n2; n1++) { delete[] neighbours[i][j][n1]; }

						delete[] neighbours[i][j];
					}

					delete[] neighbours[i];
				}

				for (int i = 0; i < n_row; i++) {
					for (int j = 0; j < n_col; j++) delete[] no_neighbours[i][j];
					delete[] no_neighbours[i];
				}

				for (int i = 0; i < n_row; i++) {
					for (int j = 0; j < n_col; j++) delete[] dispersed[i][j];
					delete[] dispersed[i];
				}

				for (int i = 0; i < n_row; i++) delete[] no_dist[i];

			}

			for (int t = 0; t < steps; t++) {
				for (int i = 0; i < n_row; i++) delete[] E[t][i];
				delete[] E[t];
			}
			delete[] E;

			for (int t = 0; t < steps; t++) {
				for (int i = 0; i < n_row; i++) delete[] L[t][i];
				delete[] L[t];
			}
			delete[] L;

			for (int t = 0; t < steps; t++) {
				for (int i = 0; i < n_row; i++) delete[] V[t][i];
				delete[] V[t];
			}
			delete[] V;

			for (int i = 0; i < n_row; i++) {
				for (int j = 0; j < n_col; j++) delete[] tau_v[i][j];
				delete[] tau_v[i];
			}
			delete[] tau_v;

			for (int i = 0; i < n_row; i++) {
				for (int j = 0; j < n_col; j++) delete[] rescue_v[i][j];
				delete[] rescue_v[i];
			}
			delete[] rescue_v;

			for (int i = 0; i < n_row; i++) delete[] pop_v[i];
			delete[] pop_v;

			for (int i = 0; i < n_row; i++) {
				for (int j = 0; j < n_col; j++) delete[] time_occupancy[i][j];
				delete[] time_occupancy[i];
			}
			delete[] time_occupancy;

			for (int i = 0; i < n_row; i++) delete[] neighbours_total[i];
			delete[] neighbours_total;

			for (int i = 0; i < n_row; i++) delete[] cumulative_pop[i];
			delete[] cumulative_pop;

			for (int i = 0; i < n_row; i++) delete[] next_v[i];
			delete[] next_v;

			for (int i = 0; i < n_row; i++) delete[] next_e[i];
			delete[] next_e;

			for (int i = 0; i < n_row; i++) delete[] next_l[i];
			delete[] next_l;

			if (alpha > 1) {

				for (int r = 0; r < runs; r++) {
					for (int i = 0; i < 2 * n_row; i++) delete[] old_K[r][i];
					delete old_K[r];
				}
				delete[] old_K;

				for (int r = 0; r < runs; r++) {
					for (int i = 0; i < 2 * n_row; i++) delete[] old_seed_v[r][i];
					delete[] old_seed_v[r];
				}
				delete[] old_seed_v;

			}

			if (alpha == divisions) {

				for (int r = 0; r < runs; r++) {
					for (int i = 0; i < n_row; i++) delete[] K[r][i];
					delete[]K[r];
				}
				delete[] K;

				for (int r = 0; r < runs; r++) {
					for (int i = 0; i < n_row; i++) delete[] seed_v[r][i];
					delete[] seed_v[r];
				}
				delete[] seed_v;
			}
		}

		for (int t = 0; t < days; t++) delete[] total_E[t];
		delete[] total_E;

		for (int t = 0; t < days; t++) delete[] total_L[t];
		delete[] total_L;

		for (int t = 0; t < days; t++) delete[] total_V[t];
		delete[] total_V;

		for (int t = 0; t < days; t++) delete[] site_avg_v[t];
		delete[] site_avg_v;

		for (int t = 0; t < days; t++) delete[] site_avg_l[t];
		delete[] site_avg_l;

		for (int t = 0; t < days; t++) delete[] site_sd_v[t];
		delete[] site_sd_v;

		for (int t = 0; t < days; t++) delete[] site_sd_l[t];
		delete[] site_sd_l;

		for (int t = 0; t < days; t++) delete[] ext_sites_v[t];
		delete[] ext_sites_v;

		for (int r = 0; r < runs; r++) delete[] seed_patch[r];
		delete[] seed_patch;

		for (int t = 0; t < days; t++) delete[] no_occupied[t];
		delete[] no_occupied;

		for (int t = 0; t < days; t++) delete[] prop_occupied[t];
		delete[] prop_occupied;

		for (int t = 0; t < days; t++) delete[] prob_established[t];
		delete[] prob_established;

		for (int t = 0; t < days; t++) delete[] count_occupancy[t];
		delete[] count_occupancy;

		for (int t = 0; t < days; t++) delete[] avg_total_V[t];
		delete[] avg_total_V;

		for (int t = 0; t < days; t++) delete[] avg_total_L[t];
		delete[] avg_total_L;

		for (int t = 0; t < days; t++) delete[] avg_total_E[t];
		delete[] avg_total_E;

		for (int t = 0; t < days; t++) delete[] avg_no_occupied[t];
		delete[] avg_no_occupied;

		for (int t = 0; t < days; t++) delete[] avg_prop_occupied[t];
		delete[] avg_prop_occupied;

		for (int t = 0; t < days; t++) delete[] avg_prob_established[t];
		delete[] avg_prob_established;

		delete[] run_avg_v;
		delete[] run_avg_l;
		delete[] run_sd_v;
		delete[] run_sd_l;
		delete[] avg_ext_sites_v;
		delete[] extinctions_v;
		delete[] ext_time_v;
		delete[] prop_extinct;
		delete[] time_occupancy_all;
		delete[] avg_count_occupancy;
		delete[] disp_fun_v;

	return 0;

}

