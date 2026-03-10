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
int num_threads = 8;
double r_dust = 1.e-3;  // dust radius, m -> require to change SRP_coe as well!!!
double Q_pr = 1.0;    // reflectivity coefficient of solar radiation pressure
int integrator_choice = 1; // Default to IAS15
double tmax; // time to end simulation, seconds
double bs_eps = 1.0e-5;    // Default BS tolerance
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

/* define output timing */
static const double output_days[] = {0.0, 64.44, 78.65, 83.77, 92.66, 114.75, 131.29, 153.47, 155.31, 177.46, 198.9, 230.39};
//static const double output_days[] = {0.0, 0.34, 0.74, 1.14, 1.74, 2.15, 3.72, 4.72, 5.70, 11.86, 14.91};

#define NUM_OUTPUTS (sizeof(output_days) / sizeof(output_days[0]))
static const int num_outputs = NUM_OUTPUTS;
static double output_sec[NUM_OUTPUTS];
static const int max_index = num_outputs - 1;
static int output_i = 0;
static double next_output_t = 0.;

// parameter tracking minimum dt within output interval
double dt_minimum = 1.e15;

int main(int argc, char* argv[]){
  /* Convert all output_days to seconds */
  for (int i = 0; i <= max_index; i++) {
    output_sec[i] = output_days[i] * 86400.0;
  }
	tmax = output_sec[max_index];
  next_output_t = output_sec[0];

	/* Parse command-line arguments */
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
			num_threads = atoi(argv[++i]);
		} else if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) {
			r_dust = atof(argv[++i]);
		} else if (strcmp(argv[i], "-qpr") == 0 && i + 1 < argc) {
			Q_pr = atof(argv[++i]);
		} else if (strcmp(argv[i], "-intg") == 0 && i + 1 < argc) {
			integrator_choice = atoi(argv[++i]);
		} else if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
			strncpy(fpath, argv[++i], sizeof(fpath));
			fpath[sizeof(fpath) - 1] = '\0';
		} else {
			fprintf(stderr, "Usage: %s -n <num_threads> -r <r_dust> -qpr <Q_pr> -intg <1|2> -f <input_file>\n", argv[0]);
			return 1;
		}
	}

	printf("Running with %d OpenMP threads\n", num_threads);
	printf("Running with r_dust = %e\n", r_dust);
	printf("Running with Q_pr = %e\n", Q_pr);
	printf("Running with tmax = %f (automatically set from output_days)\n", tmax);
	printf("Running with dust input file = %s\n", fpath);

	/* Set the number of OpenMP threads to be the number of processors */
	//int np = omp_get_num_procs();
	omp_set_num_threads(num_threads);
	
	// Setup simulation structure
	struct reb_simulation* r = reb_simulation_create();

	// Setup constants
	if (integrator_choice == 2) {
		r->integrator = REB_INTEGRATOR_BS;
		r->ri_bs.eps_rel       = bs_eps;
		r->ri_bs.eps_abs       = bs_eps;
		printf("Using Bulirsch-Stoer integrator (eps = %e)\n", r->ri_bs.eps_rel);
	} else {
		r->integrator = REB_INTEGRATOR_IAS15;
		printf("Using IAS15 integrator\n");
	}
	r->dt                  = 1e1;    // Initial timestep, s
	r->N_active            = 3;     // Only the Sun and the Didymos-Dimorphos system are massive, the dust particles are treated as test particles
	r->additional_forces   = force_radiation;
	r->heartbeat           = heartbeat;
	r->G                   = G_const;
	
	reb_simulation_configure_box(r,sep_system*5.,1,1,1);    

	// Didymos
	double mass_didy = mass_system*vol_didy/(vol_didy+vol_dimor);
	double mass_dimor = mass_system*vol_dimor/(vol_didy+vol_dimor);
	struct reb_particle Didymos = {0};
	double r_didy_com = -vol_dimor*sep_system/(vol_didy+vol_dimor); // distance of Didymos to center of mass of the system
	double v_didy_com = sqrt(-r->G*mass_dimor/pow(sep_system,2.0)*r_didy_com);
	Didymos.m    = mass_didy;
	Didymos.r    = 850.0/2.0;
	Didymos.x    = r_didy_com + T11 * (r_DSB_t1.x - r_DSB_t0.x) + T12 * (r_DSB_t1.y - r_DSB_t0.y) + T13 * (r_DSB_t1.z - r_DSB_t0.z);
	Didymos.y    = 0.0        + T21 * (r_DSB_t1.x - r_DSB_t0.x) + T22 * (r_DSB_t1.y - r_DSB_t0.y) + T23 * (r_DSB_t1.z - r_DSB_t0.z);
	Didymos.z    = 0.0        + T31 * (r_DSB_t1.x - r_DSB_t0.x) + T32 * (r_DSB_t1.y - r_DSB_t0.y) + T33 * (r_DSB_t1.z - r_DSB_t0.z);
	Didymos.vx   =              T11 * v_DSB_t1.x + T12 * v_DSB_t1.y + T13 * v_DSB_t1.z;
	Didymos.vy   = v_didy_com + T21 * v_DSB_t1.x + T22 * v_DSB_t1.y + T23 * v_DSB_t1.z;
	Didymos.vz   =              T31 * v_DSB_t1.x + T32 * v_DSB_t1.y + T33 * v_DSB_t1.z;
	Didymos.hash = 1;
	reb_simulation_add(r, Didymos);
    
	// Dimorphos
	struct reb_particle Dimorphos = {0};
	double r_dimor_com = vol_didy*sep_system/(vol_didy+vol_dimor);
	double v_dimor_com = -sqrt(r->G*mass_didy/pow(sep_system,2.0)*r_dimor_com);
	Dimorphos.m    = mass_dimor;
	Dimorphos.r    = 175.0/2.0;
	Dimorphos.x    = r_dimor_com + T11 * (r_DSB_t1.x - r_DSB_t0.x) + T12 * (r_DSB_t1.y - r_DSB_t0.y) + T13 * (r_DSB_t1.z - r_DSB_t0.z);
	Dimorphos.y    = 0.0         + T21 * (r_DSB_t1.x - r_DSB_t0.x) + T22 * (r_DSB_t1.y - r_DSB_t0.y) + T23 * (r_DSB_t1.z - r_DSB_t0.z);
	Dimorphos.z    = 0.0         + T31 * (r_DSB_t1.x - r_DSB_t0.x) + T32 * (r_DSB_t1.y - r_DSB_t0.y) + T33 * (r_DSB_t1.z - r_DSB_t0.z);
	Dimorphos.vx   =               T11 * v_DSB_t1.x + T12 * v_DSB_t1.y + T13 * v_DSB_t1.z;
	Dimorphos.vy   = v_dimor_com + T21 * v_DSB_t1.x + T22 * v_DSB_t1.y + T23 * v_DSB_t1.z;
	Dimorphos.vz   =               T31 * v_DSB_t1.x + T32 * v_DSB_t1.y + T33 * v_DSB_t1.z;
	Dimorphos.hash = 2;
	reb_simulation_add(r, Dimorphos);
	
	// Sun
	struct reb_particle star = {0};
	star.m  = mass_star;
	star.x = - T11 * r_DSB_t0.x - T12 * r_DSB_t0.y - T13 * r_DSB_t0.z;
	star.y = - T21 * r_DSB_t0.x - T22 * r_DSB_t0.y - T23 * r_DSB_t0.z;
	star.z = - T31 * r_DSB_t0.x - T32 * r_DSB_t0.y - T33 * r_DSB_t0.z;
	star.hash = 3;
	reb_simulation_add(r, star);

	// Earth (Hubble Space Telescope)
	struct reb_particle Earth = {0};
	Earth.m  = 0.;
	Earth.x = T11 * r_Earth_t1.x + T12 * r_Earth_t1.y + T13 * r_Earth_t1.z + star.x;
	Earth.y = T21 * r_Earth_t1.x + T22 * r_Earth_t1.y + T23 * r_Earth_t1.z + star.y;
	Earth.z = T31 * r_Earth_t1.x + T32 * r_Earth_t1.y + T33 * r_Earth_t1.z + star.z;
	Earth.vx = T11 * v_Earth_t1.x + T12 * v_Earth_t1.y + T13 * v_Earth_t1.z;
	Earth.vy = T21 * v_Earth_t1.x + T22 * v_Earth_t1.y + T23 * v_Earth_t1.z;
	Earth.vz = T31 * v_Earth_t1.x + T32 * v_Earth_t1.y + T33 * v_Earth_t1.z;
	Earth.hash = 4;
	reb_simulation_add(r, Earth);

	unsigned int N_particles = 4; // current number of particles (didy, dimor, sun, earth)
	unsigned int N_scanned = 0;   // record how many particles scanned in the input particle file
	unsigned int N_didy = 0;      // record initial # of particles within radius of Didymos
	unsigned int N_dimor = 0;     // record initial # of particles within radius of Dimorphos
	unsigned int N_hill = 0;      // record initial # of particles farther than Hill radius

	// Dust particles
	if (1){
		// open dust particles file
		FILE *dust_file = fopen(fpath, "r");
		if (dust_file == NULL) {
			fprintf(stderr, "Error: Could not open file %s\n", fpath);
			return 1;
		}

		// open a file examine the particles that are deleted
		FILE *f_dp = fopen("deleted_particles.csv", "w");
		if (f_dp == NULL) {
			reb_simulation_error(r, "Could not open file: deleted_particles.csv");
			return 1;
		}

		double disSQ_Didy, disSQ_Dimor;
		ReadParticle rp;
		while (fscanf(dust_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&rp.ID, &rp.x, &rp.y, &rp.z, &rp.vx, &rp.vy, &rp.vz, &rp.mass, &rp.density) == 9) {

			N_scanned++;
			// rotate the original coordinate system around its y-axis by 180 degree
			transform(&rp);

			struct reb_particle p = {0};
			p.m = 0.0;
			p.r = r_dust;
			p.x = rp.x + Dimorphos.x;
			p.y = rp.y + Dimorphos.y;
			p.z = rp.z + Dimorphos.z;
			p.vx = rp.vx + Dimorphos.vx;
			p.vy = rp.vy + Dimorphos.vy;
			p.vz = rp.vz + Dimorphos.vz;

			disSQ_Didy  = pow(p.x-Didymos.x,2) + pow(p.y-Didymos.y,2) + pow(p.z-Didymos.z,2);
			disSQ_Dimor = pow(p.x-Dimorphos.x,2) + pow(p.y-Dimorphos.y,2) + pow(p.z-Dimorphos.z,2);

			// skip particles that are farther than hill radius and that make up Didymos or Dimorphos
			if (disSQ_Didy < Rsq_didy){
				fprintf(f_dp, "%f,%f,%f,%d\n", rp.x, rp.y, rp.z, 1);
				N_didy++;
				continue;
			}
			if (disSQ_Dimor < Rsq_dimor){
				fprintf(f_dp, "%f,%f,%f,%d\n", rp.x, rp.y, rp.z, 2);
				N_dimor++;
				continue;
			}
			if (disSQ_Didy > Rsq_hill){
				fprintf(f_dp, "%f,%f,%f,%d\n", rp.x, rp.y, rp.z, 3);
				N_hill++;
				continue;
			}

			N_particles++;
			p.hash = N_particles;
			reb_simulation_add(r, p);
		}
		fclose(dust_file);
		fclose(f_dp);
	}

	fprintf(stdout, "Total # of particles scanned: %i\n", N_scanned);
	fprintf(stdout, "Total # of particles registered: %i\n", N_particles);
	fprintf(stdout, "# of particles within Didymos: %i\n", N_didy);
	fprintf(stdout, "# of particles within Dimorphos: %i\n", N_dimor);
	fprintf(stdout, "# of particles outside of Hill radius: %i\n", N_hill);

	system("rm -v particles.txt");
	system("rm -v collide.txt");

	reb_simulation_save_to_file_interval(r, "archive.bin", 864000.); // save for restart. 10 days between snapshots
//	reb_simulation_integrate(r, tmax);
	for (int i = 0; i < num_outputs; i++) {
		reb_simulation_integrate(r, output_sec[i]);
	}
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

void reb_simulation_move_to_DSB(struct reb_simulation* const r){
	const int N_real = r->N - r->N_var;
	if (N_real>0){
		struct reb_particle* restrict const particles = r->particles;
		struct reb_particle Didy = r->particles[0];
		struct reb_particle Dimor = r->particles[1];
		// position and velocity of the center of mass of Didymos and Dimorphos
		double com_x = (Didy.m * Didy.x + Dimor.m * Dimor.x) / (Didy.m + Dimor.m);
		double com_y = (Didy.m * Didy.y + Dimor.m * Dimor.y) / (Didy.m + Dimor.m);
		double com_z = (Didy.m * Didy.z + Dimor.m * Dimor.z) / (Didy.m + Dimor.m);
		// Note: Variational particles will not be affected.
		for (int i=0;i<N_real;i++){
			particles[i].x  -= com_x;
			particles[i].y  -= com_y;
			particles[i].z  -= com_z;
		}
	}
}

void heartbeat(struct reb_simulation* r){
//----------------track minimum dt--------------------
	if (r->dt < dt_minimum){
		dt_minimum = r->dt;
	}
	
//----------------output dt history------------------
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
		
		reb_simulation_move_to_DSB(r);
	}

//----------------output all particles---------------------
//	if(reb_simulation_output_check(r, next_output_t)){
	if (output_i <= max_index && r->t >= output_sec[output_i] - 1e-6) {
		reb_simulation_output_timing(r, tmax);
		printf("\n");

		struct reb_particle* particles = r->particles;
		const int N = r->N;
		double hash_val;

		// Calculate the CoM of Didymos and Dimorphos
		struct reb_particle Didy = particles[0];
		struct reb_particle Dimor = particles[1];
		double total_m = Didy.m + Dimor.m;

		double com_x  = (Didy.m * Didy.x  + Dimor.m * Dimor.x)  / total_m;
		double com_y  = (Didy.m * Didy.y  + Dimor.m * Dimor.y)  / total_m;
		double com_z  = (Didy.m * Didy.z  + Dimor.m * Dimor.z)  / total_m;
		double com_vx = (Didy.m * Didy.vx + Dimor.m * Dimor.vx) / total_m;
		double com_vy = (Didy.m * Didy.vy + Dimor.m * Dimor.vy) / total_m;
		double com_vz = (Didy.m * Didy.vz + Dimor.m * Dimor.vz) / total_m;

		// Open file for binary append
		FILE* fp = fopen("particles.txt","ab+");
		if ( fp == NULL){
			reb_simulation_error(r, "Can not open file: particles.txt.");
			return;
		}
		
		fwrite(&(N), sizeof(int), 1, fp);
		fwrite(&(r->t), sizeof(double), 1, fp);
		fwrite(&(r_dust), sizeof(double), 1, fp);
		for ( int i=0; i<N; i++ ) {
			const struct reb_particle p = particles[i];

			hash_val = (double)p.hash;
			double rx = p.x  - com_x;
			double ry = p.y  - com_y;
			double rz = p.z  - com_z;
			double rvx = p.vx - com_vx;
			double rvy = p.vy - com_vy;
			double rvz = p.vz - com_vz;

			fwrite( &hash_val, sizeof(double), 1, fp);
			fwrite( &rx,  sizeof(double), 1, fp);
			fwrite( &ry,  sizeof(double), 1, fp);
			fwrite( &rz,  sizeof(double), 1, fp);
			fwrite( &rvx, sizeof(double), 1, fp);
			fwrite( &rvy, sizeof(double), 1, fp);
			fwrite( &rvz, sizeof(double), 1, fp);
		}
		fclose(fp);

		// update next output timing
		output_i++;
		if (output_i <= max_index) {
			next_output_t = output_sec[output_i];
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
