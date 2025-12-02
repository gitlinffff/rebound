/**
 * Modification log by Linfel
 * 1.add a 'ReadParticle' struct and a 'transform' function for processing initial
 * frame transformation of input dust particles.
 * 2. turn off force_radiation
 * 3. mass of dust particles are set to 0 for now
 */

/**
 * Dust evolution after DART impact
 *
 * This example shows how to integrate dust particles
 * using the IAS15 integrator with additional forces.
 * The example sets the function pointer `additional_forces`
 * to a function that describes the radiation forces and
 * non-gravitational perturbations of Dimorphos and Didymos.
 *
 * The output is custom too, as given in the function  `heartbeat`.
 * 
 * Written by Yun Zhang (2023 May)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <omp.h>
#include "rebound.h"

// define a ReadParticle struct for holding read-in particle data
typedef struct {
	double ID;
	double x, y, z;
	double vx, vy, vz;
	double mass, density;
	double a, e, i, Omega, omega, f;
} ReadParticle;

// math
typedef struct {
	double x;
	double y;
	double z;
} Vector3;

Vector3 crossProduct(Vector3 a, Vector3 b) {
	Vector3 axb;
	axb.x = a.y * b.z - a.z * b.y;
	axb.y = a.z * b.x - a.x * b.z;
	axb.z = a.x * b.y - a.y * b.x;
	return axb;
}

double vectorNorm(Vector3 v) {
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// Other function declarations
void transform(ReadParticle *p);
void force_radiation(struct reb_simulation* r);
void heartbeat(struct reb_simulation* r);

// Declare input variable (no const!)
int num_threads;
double r_dust;  // dust particle radius, m -> require to change SRP_coe as well!!!
double Q_pr;    // reflectivity coefficient of solar radiation pressure
double bs_eps;  // BS integrator relative and absolute tolerances
double tmax;    // time to end simulation, seconds
char fpath[256];// to store input data file path

// Define constants
const double G_const = 6.6743e-11; //  m^3 / kg s^2
const double AU = 1.495978707e11;
const double mass_star = 1.9884e30;
const double radius_star = 6.957e8;
const double mass_system = 5.5e11; // kg
const double sep_system = 1170.0; // seperation m
const double vol_didy = 0.2295409644951028;  // km^3
const double vol_dimor = 0.001830603200702610;
const double rho_dust = 3000; // dust particle density, kg/m^3
const double Rsq_didy = 850.0/2.0 * 850.0/2.0;
const double Rsq_dimor = 175.0/2.0 * 175.0/2.0;
//const double Rsq_long_dimor = 193.0/2.0 * 193.0/2.0;  // use its longest dimension
//const double Rsq_long_dimor = 177.0/2.0 * 177.0/2.0;  // use its longest dimension (for data_high)
const double Rsq_hill = 70500.0*70500.0;  // twice Hill radius of D-D system, m

// J2
const double J2_didy = 0.0956324486828653;  // J2 of Didymos  check and adjust
const double J2_dimor = 0.113929814552540;  // J2 of Dimorphos

// SRP
const double c = 2.99792458e8;         // speed of light.
const double Fsun = 1367.0;  /* integrated stellar flux at 1 au, W/m^2 */

/* 2022-Sep-26 23:14:24.1830 UTC (moment of the impact)
 * (https://ssd.jpl.nasa.gov/horizons/app.html#/) 
 * Coordinate Center: Sun (body center) [500@10] */
const Vector3 r_DSB_t0 = {1.556582267294774E+11, 1.349129152256939E+10, -8.638156994081538E+09}; // position of Didymos System Barycenter
const Vector3 v_DSB_t0 = {-7.322785629453647E+03, 3.319798419238497E+04, 9.918308846239352E+02}; // velocity of Didymos System Barycenter
const Vector3 r_Dimor_t0 = {1.556582259038835E+11, 1.349129068894670E+10, -8.638156976265389E+09};
const Vector3 v_Dimor_t0 = {-7.322906365387563E+03, 3.319810316934613E+04, 9.918029911906725E+02};

/* 2022-Sep-26 23:17:24.1830 UTC (180s after the impact)
 * (https://ssd.jpl.nasa.gov/horizons/app.html#/) 
 * Coordinate Center: Sun (body center) [500@10] */
const Vector3 r_DSB_t1 = {1.556569085406981E+11, 1.349726715213361E+10, -8.637978459636498E+09}; // position of Didymos System Barycenter
const Vector3 v_DSB_t1 = {-7.323756364209678E+03, 3.319789984814826E+04, 9.918851713181773E+02}; // velocity of Didymos System Barycenter
const Vector3 r_Dimor_t1 = {1.556569076939493E+11, 1.349726633984210E+10, -8.637978446769070E+09};
const Vector3 v_Dimor_t1 = {-7.323872270649070E+03, 3.319801993581914E+04, 9.918576504499956E+02};
const Vector3 r_Earth_t1 = {1.496913045881429E+11, 9.256618023871835E+09, -1.483832179474644E+06};
const Vector3 v_Earth_t1 = {-2.326419493617441E+03, 2.963124474205269E+04, -4.788949253207164E-01};

// matrix T that converts a vector from Sun body center frame to the Didymos system barycenter frame
const double T11 = -0.703595792353257;
const double T12 = -0.710438191316344;
const double T13 = 0.015183454875444;
const double T21 = 0.702824020343859;
const double T22 = -0.692584740738747;
const double T23 = 0.16237232930379;
const double T31 = -0.104839674791979;
const double T32 = 0.124915784491013;
const double T33 = 0.986612735258626;

// parameter tracking minimum dt within output interval
double dt_minimum = 1.e15;

// define output timing
static const double output_days[] = {0.0, 64.44, 78.65, 83.77, 92.66, 114.75, 131.29, 153.47, 155.31, 177.46, 198.9, 230.39};
//static const double output_days[] = {0.0, 0.01, 0.03, 0.05, 0.09, 0.1};

#define NUM_OUTPUTS (sizeof(output_days) / sizeof(output_days[0]))
static const int num_outputs = NUM_OUTPUTS;
static double output_sec[NUM_OUTPUTS];
static const int max_index = num_outputs - 1;
static int output_i = 0;
static double next_output_t = 0.;

int main(int argc, char* argv[]){
	
	// Parse command-line arguments
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
			num_threads = atoi(argv[++i]);
		} else if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) {
			r_dust = atof(argv[++i]);
		} else if (strcmp(argv[i], "-qpr") == 0 && i + 1 < argc) {
			Q_pr = atof(argv[++i]);
		} else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
			tmax = atof(argv[++i]);
		} else if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
			strncpy(fpath, argv[++i], sizeof(fpath));
			fpath[sizeof(fpath) - 1] = '\0'; // null-terminate safely
		} else {
			fprintf(stderr, "Usage: %s -n <num_threads> -r <r_dust> -qpr <Q_pr> -t <tmax> -f <rst_archive>\n", argv[0]);
			return 1;
		}
	}
	printf("Running with %d OpenMP threads\n", num_threads);
	printf("Running with r_dust = %.8e\n", r_dust);
	printf("Running with Q_pr = %e\n", Q_pr);
	printf("Running with tmax = %.2f\n", tmax);
	printf("Restarting with archive: %s\n", fpath);

	// Convert all output_days to seconds
	for (int i = 0; i <= max_index; i++) {
		output_sec[i] = output_days[i] * 86400.0;
	}
	next_output_t = output_sec[0];

	// Set the number of OpenMP threads to be the number of processors
	//int np = omp_get_num_procs();
	omp_set_num_threads(num_threads);
	
	// restart from a specified snapshot
	struct reb_simulationarchive* archive = reb_simulationarchive_create_from_file(fpath);// "archive.bin"
	struct reb_simulation* r = reb_simulation_create_from_simulationarchive(archive, 6); // -1 if the last snapshot
	reb_simulationarchive_free(archive);

	// print restarting information
	printf("===========================\nRestarting information:\n");
	printf("t: %f s = %f day\n", r->t, r->t/86400.);
	printf("dt: %f\n", r->dt);
	printf("N_active: %d\n", r->N_active);
	printf("G: %e\n", r->G);
	printf("BS integrator tolerance: %e, %e\n", r->ri_bs.eps_rel, r->ri_bs.eps_abs);
	printf("integrator = %d\n", r->integrator);
	struct reb_particle* particles = r->particles;
	const struct reb_particle p = particles[4];
	printf("p.r: %.8e\n", p.r);
	if (p.r != r_dust){
    char error_msg[70];
    sprintf(error_msg, "r_dust %.8e does not match with archive %.8e", r_dust, p.r);
    reb_simulation_error(r, error_msg);
		return 1;
	}

	// Reset function pointers:
	r->additional_forces   = force_radiation;
	r->heartbeat           = heartbeat;

	// --- LOGIC TO RESTART OUTPUT TIMING ---
	printf("Restarting from simulation time: %.2f days (%.2f seconds)\n", r->t/86400., r->t);

	// Find the next output time index (output_i)
	for (int i = 0; i <= max_index; i++) {
		if (output_sec[i] > r->t) {
			output_i = i;
			next_output_t = output_sec[i];
			printf("Resuming output from index: %d\n", output_i);
			printf("The next scheduled output time is: %.2f days (%.2f seconds)\n", 
						 next_output_t/86400., next_output_t);
			break; // Stop when the next future time is found
		}
	}
	// Handle the case where the simulation ran past the last defined output time  //
	// still some problem here
	if (r->t > output_sec[max_index]) {
		printf("All defined output times (up to %.2f seconds) have been passed.\n", output_sec[max_index]);
		// Set next_output_t to a large value to disable outputs
		next_output_t = tmax * 10.0; 
	}

	// start integration
	reb_simulation_save_to_file_interval(r, "archive1.bin", 864000.); // save for restart. 10 days between snapshots
	reb_simulation_integrate(r, tmax);
	fprintf(stdout, "\n");
}

void transform(ReadParticle *p) {
	// Units conversion from cgs to SI
	double new_x       = p->x / 100.0;
	double new_y       = p->y / 100.0;
	double new_z       = p->z / 100.0;
	double new_vx      = p->vx / 100.0;
	double new_vy      = p->vy / 100.0;
	double new_vz      = p->vz / 100.0;
	double new_mass    = p->mass / 1000.0;
	double new_density = p->density * 1000.0;

	// Rotate frame 180 degrees around y-axis
	new_x = -new_x;
	new_z = -new_z;
	new_vx = -new_vx;
	new_vz = -new_vz;

	// Update particle's info
	p->x = new_x;
	p->y = new_y;
	p->z = new_z;
	p->vx = new_vx;
	p->vy = new_vy;
	p->vz = new_vz;
	p->mass = new_mass;
	p->density = new_density;
}

void force_radiation(struct reb_simulation* r){
	double SRP_coe = Q_pr * Fsun/c * 3.0/4.0/rho_dust;

	struct reb_particle* particles = r->particles;
	const struct reb_particle Didymos = particles[0];
	const struct reb_particle Dimorphos = particles[1];
	const struct reb_particle star = particles[2];            // cache
	const int N = r->N;
	
	// Sun-Didymos vector
	double prx_didy_star = Didymos.x-star.x;
	double pry_didy_star = Didymos.y-star.y;
	double prz_didy_star = Didymos.z-star.z;
	double pr = 1.0/sqrt(prx_didy_star*prx_didy_star + pry_didy_star*pry_didy_star + prz_didy_star*prz_didy_star);
	prx_didy_star *= pr;
	pry_didy_star *= pr;
	prz_didy_star *= pr;
	
	// Sun-Dimorphos vector
	double prx_dimor_star = Dimorphos.x-star.x;
	double pry_dimor_star = Dimorphos.y-star.y;
	double prz_dimor_star = Dimorphos.z-star.z;
	pr = 1.0/sqrt(prx_dimor_star*prx_dimor_star + pry_dimor_star*pry_dimor_star + prz_dimor_star*prz_dimor_star);
	prx_dimor_star *= pr;
	pry_dimor_star *= pr;
	prz_dimor_star *= pr;
	
	
#pragma omp parallel for
	for (int i=0;i<N;i++){

		struct reb_particle p = particles[i];             // cache
		if ( p.m > 0. ) continue;                         // Only dust particles feel radiation forces
		
		// set up vector
		const double prx_star  = p.x-star.x;
		const double pry_star  = p.y-star.y;
		const double prz_star  = p.z-star.z;
		const double prx_didy = p.x-Didymos.x;
		const double pry_didy = p.y-Didymos.y;
		const double prz_didy = p.z-Didymos.z;
		const double prx_dimor = p.x-Dimorphos.x;
		const double pry_dimor = p.y-Dimorphos.y;
		const double prz_dimor = p.z-Dimorphos.z;
		double dfactor;
		unsigned int flag_noshadow = 1;
		
		/* radiation force */
		// check if is in the shadow of Didymos
		pr = prx_didy*prx_didy_star + pry_didy*pry_didy_star + prz_didy*prz_didy_star;
		if ( pr > 0.0 ) {
			double prx_body_rad = prx_didy - pr*prx_didy_star;  // radial vector of dust relative to the Sun-Didymos direction
			double pry_body_rad = pry_didy - pr*pry_didy_star;  // radial vector of dust relative to the Sun-Didymos direction
			double prz_body_rad = prz_didy - pr*prz_didy_star;  // radial vector of dust relative to the Sun-Didymos direction
			if ( prx_body_rad*prx_body_rad + pry_body_rad*pry_body_rad + prz_body_rad*prz_body_rad < Rsq_didy )
				flag_noshadow = 0;
		}
		
		// check if is in the shadow of Dimorphos
		pr = prx_dimor*prx_dimor_star + pry_dimor*pry_dimor_star + prz_dimor*prz_dimor_star;
		if ( pr > 0.0 ) {
			double prx_body_rad = prx_dimor - pr*prx_dimor_star;  // radial vector of dust relative to the Sun-Dimorphos direction
			double pry_body_rad = pry_dimor - pr*pry_dimor_star;  // radial vector of dust relative to the Sun-Dimorphos direction
			double prz_body_rad = prz_dimor - pr*prz_dimor_star;  // radial vector of dust relative to the Sun-Dimorphos direction
			if ( prx_body_rad*prx_body_rad + pry_body_rad*pry_body_rad + prz_body_rad*prz_body_rad < Rsq_dimor )
				flag_noshadow = 0;
		}
		
		// add radiation pressure if not in shadow
		if ( flag_noshadow ) {
			pr = sqrt(prx_star*prx_star + pry_star*pry_star + prz_star*prz_star);     // distance relative to star
			//const double prvx = p.vx-star.vx;
			//const double prvy = p.vy-star.vy;
			//const double prvz = p.vz-star.vz;
			//const double rdot = (prvx*prx_star + prvy*pry_star + prvz*prz_star)/pr;     // radial velocity relative to star
			dfactor = SRP_coe/r_dust * pow(AU/pr,2.0);

			// Equation (5) of Burns, Lamy, Soter (1979)
			//particles[i].ax += dfactor*((1.-rdot/c)*prx_star/pr - prvx/c);
			//particles[i].ay += dfactor*((1.-rdot/c)*pry_star/pr - prvy/c);
			//particles[i].az += dfactor*((1.-rdot/c)*prz_star/pr - prvz/c);
			particles[i].ax += dfactor*prx_star/pr;
			particles[i].ay += dfactor*pry_star/pr;
			particles[i].az += dfactor*prz_star/pr;
		}
		
		// J2 of Didymos
		pr   = prx_didy*prx_didy + pry_didy*pry_didy + prz_didy*prz_didy;
		dfactor  = 3.0*r->G*J2_didy*Didymos.m*Didymos.r*Didymos.r/2./pow(pr,3.5);
		particles[i].ax += dfactor*prx_didy*(prx_didy*prx_didy + pry_didy*pry_didy - 4.*prz_didy*prz_didy);
		particles[i].ay += dfactor*pry_didy*(prx_didy*prx_didy + pry_didy*pry_didy - 4.*prz_didy*prz_didy);
		particles[i].az += dfactor*prz_didy*(3.*(prx_didy*prx_didy + pry_didy*pry_didy) - 2.*prz_didy*prz_didy);
		
		// J2 of Dimorphos
		pr   = prx_dimor*prx_dimor + pry_dimor*pry_dimor + prz_dimor*prz_dimor;
		dfactor  = 3.0*r->G*J2_dimor*Dimorphos.m*Dimorphos.r*Dimorphos.r/2./pow(pr,3.5);
		particles[i].ax += dfactor*prx_dimor*(prx_dimor*prx_dimor + pry_dimor*pry_dimor - 4.*prz_dimor*prz_dimor);
		particles[i].ay += dfactor*pry_dimor*(prx_dimor*prx_dimor + pry_dimor*pry_dimor - 4.*prz_dimor*prz_dimor);
		particles[i].az += dfactor*prz_dimor*(3.*(prx_dimor*prx_dimor + pry_dimor*pry_dimor) - 2.*prz_dimor*prz_dimor);
	}
}


void reb_move_to_Didymos(struct reb_simulation* const r){
	const int N_real = r->N - r->N_var;
	if (N_real>0){
		struct reb_particle* restrict const particles = r->particles;
		struct reb_particle hel = r->particles[0];
		// Note: Variational particles will not be affected.
		for (int i=1;i<N_real;i++){
			particles[i].x  -= hel.x;
			particles[i].y  -= hel.y;
			particles[i].z  -= hel.z;
		}
		r->particles[0].x = 0.;
		r->particles[0].y = 0.;
		r->particles[0].z = 0.;
	}
}

void reb_simulation_move_to_DidyDimor_com(struct reb_simulation* const r){
	const int N_real = r->N - r->N_var;
	if (N_real>0){
		struct reb_particle* restrict const particles = r->particles;
		struct reb_particle Didy = r->particles[0];
		struct reb_particle Dimor = r->particles[1];
		// position and velocity of the center of mass of Didymos and Dimorphos
		double com_x = (Didy.m * Didy.x + Dimor.m * Dimor.x) / (Didy.m + Dimor.m);
		double com_y = (Didy.m * Didy.y + Dimor.m * Dimor.y) / (Didy.m + Dimor.m);
		double com_z = (Didy.m * Didy.z + Dimor.m * Dimor.z) / (Didy.m + Dimor.m);
		double com_vx = (Didy.m * Didy.vx + Dimor.m * Dimor.vx) / (Didy.m + Dimor.m);
		double com_vy = (Didy.m * Didy.vy + Dimor.m * Dimor.vy) / (Didy.m + Dimor.m);
		double com_vz = (Didy.m * Didy.vz + Dimor.m * Dimor.vz) / (Didy.m + Dimor.m);
		// Note: Variational particles will not be affected.
		for (int i=0;i<N_real;i++){
			particles[i].x  -= com_x;
			particles[i].y  -= com_y;
			particles[i].z  -= com_z;
			particles[i].vx  -= com_vx;
			particles[i].vy  -= com_vy;
			particles[i].vz  -= com_vz;
		}
	}
}

void heartbeat(struct reb_simulation* r){
//----------------track minimum dt--------------------
	if (r->dt < dt_minimum){
		dt_minimum = r->dt;
	}
	
//----------------output dt history-------------------
	if(reb_simulation_output_check(r, 60.0)){
		int N_tot = r->N;

		// Open file in append mode
		char filename[20] = "dt_history.csv";
		FILE* f_dt = fopen(filename, "a");
		if (f_dt == NULL) {
			char error_msg[50];
			sprintf(error_msg, "Could not open file: %s", filename);
			reb_simulation_error(r, error_msg);
			return;
		}

		// If file is empty, print header
		static int header_written = 0;
		if (!header_written){
			fprintf(f_dt, "N_tot, t, dt, dt_minimum, t/tmax%%\n");
			header_written = 1;
		}

		// Write values
		fprintf(f_dt, "%-10d %-15.6f %-15.6f %-15.6f %-8.4f\n", N_tot, r->t, r->dt, dt_minimum, r->t/tmax*100.0);

		fclose(f_dt);
		
		// reset dt_minimum
		dt_minimum = 1.e15;
	}
	
//----------------remove collided particles-----------------
	if(reb_simulation_output_check(r, 60.0)){  
		// In reality, dt is larger than 60 s. This chunk of code is executed every time steps

		struct reb_particle* particles = r->particles;
		const struct reb_particle Didymos = particles[0];
		const struct reb_particle Dimorphos = particles[1];
		int N = r->N;
		
		double dDisSQ_Didy, dDisSQ_Dimor;
		unsigned int N_remove = 0;
		unsigned int flag_remove;
		
		// delete and record collided particles
		FILE* f_c = fopen("collide.txt","ab+");
		if ( f_c == NULL){
			reb_simulation_error(r, "Can not open file: collide.txt.");
			return;
		}

		for ( int i=0;i<N;i++ ) {

			const struct reb_particle p = particles[i-N_remove];       // cache
			if ( p.m > 0. ) continue;                                  // Only delete dust particles
			
			dDisSQ_Didy  = pow(p.x-Didymos.x,2) + pow(p.y-Didymos.y,2) + pow(p.z-Didymos.z,2);
			dDisSQ_Dimor = pow(p.x-Dimorphos.x,2) + pow(p.y-Dimorphos.y,2) + pow(p.z-Dimorphos.z,2);
			
			flag_remove = 0;
			if (dDisSQ_Didy < Rsq_didy)
				flag_remove = 1; // collide with Didymos
			else if (dDisSQ_Dimor < Rsq_dimor)
				flag_remove = 2; // collide with Dimorphos
			else if ( dDisSQ_Didy > Rsq_hill )
				flag_remove = 3; // escaped particles
					
			if ((flag_remove == 1) || (flag_remove == 2)) {
				fwrite( &(flag_remove), sizeof(int), 1, f_c );
				fwrite( &(p.hash), sizeof(int), 1, f_c );
				fwrite( &(r->t), sizeof(double), 1, f_c );
				reb_simulation_remove_particle( r, i-N_remove, 1 );
				N_remove++;
			}
		}
		fclose(f_c);
		
		reb_simulation_move_to_DidyDimor_com(r);
		//reb_simulation_move_to_hel(r);
		//reb_move_to_Didymos(r);
	}
    
//----------------output all particles---------------------
	if(reb_simulation_output_check(r, next_output_t)){
		struct reb_particle* particles = r->particles;
		const int N = r->N;
		double di;

		reb_simulation_output_timing(r, tmax);
		printf("\n");

		// output particle position and velocity
		FILE* fp = fopen("particles.txt","ab+");
		if ( fp == NULL){
			reb_simulation_error(r, "Can not open file: particles.txt.");
			return;
		}
		
		fwrite( &(N), sizeof(int), 1, fp);
		fwrite( &(r->t), sizeof(double), 1, fp);
		fwrite( &(r_dust), sizeof(double), 1, fp);
		for ( int i=0; i<N; i++ ) {
			const struct reb_particle p = particles[i];
			di = (double)p.hash;
			fwrite( &(di), sizeof(double), 1, fp);
			fwrite( &(p.x), sizeof(double), 1, fp);
			fwrite( &(p.y), sizeof(double), 1, fp);
			fwrite( &(p.z), sizeof(double), 1, fp);
			fwrite( &(p.vx), sizeof(double), 1, fp);
			fwrite( &(p.vy), sizeof(double), 1, fp);
			fwrite( &(p.vz), sizeof(double), 1, fp);
		}
		fclose(fp);
		
		// Update next output time only when the r->t passes the next target time
		if (r->t >= next_output_t && output_i <= max_index) {
			output_i++;
			if (output_i <= max_index) {
				next_output_t = output_sec[output_i];
			} else{
				next_output_t = tmax * 10.;
			}
		}
 }
    
//----------------output orbital parameters--------------------
//particles orbits relative to barycenter of binary system
	if(reb_simulation_output_check(r, 4320000.0)){
		struct reb_particle* particles = r->particles;
		const struct reb_particle Didymos = particles[0];
		const struct reb_particle Dimorphos = particles[1];
		const int N = r->N;
		struct reb_orbit orbit;

		// open a file recording particles' orbital elements
		char filename[20];
		sprintf(filename, "a_e_t%d.csv", (int)r->t);
		FILE *f_ae = fopen(filename, "w");
		if (f_ae == NULL) {
			char error_msg[50];
			sprintf(error_msg, "Could not open file: %s", filename);
			reb_simulation_error(r, error_msg);
			return;
		}
		fprintf(f_ae, "ID,a_p,e_p,x,y,z,vx,vy,vz\n");

		// create a virtual body representing the com of Didymos and Dimorphos system
		struct reb_particle virtual_com;
		virtual_com.m = Didymos.m + Dimorphos.m;
		virtual_com.x = 0.;
		virtual_com.y = 0.;
		virtual_com.z = 0.;
		virtual_com.vx = 0.;
		virtual_com.vy = 0.;
		virtual_com.vz = 0.;

		for ( int i=0;i<N;i++ ) { 
			const struct reb_particle p = particles[i];
			orbit = reb_orbit_from_particle(r->G, p, virtual_com);
			fprintf(f_ae, "%d,%f,%f,%f,%f,%f,%f,%f,%f\n", p.hash, orbit.a, orbit.e, p.x, p.y, p.z, p.vx, p.vy, p.vz);
			//a and e may not be correct for Didy and Dimor
		}
			
		fclose(f_ae);
	}
}
