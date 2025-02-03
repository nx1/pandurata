#include "../exe-dir/header.hh"
#include "../Utilities/SpectrumArray.hh"
#include "../PostProcess/PostProcess.hh"
#include "../PostProcess/MetricPP.hh"
#include "../PhysicsHelper/Mapper.hh"
#include <stdio.h>
#include <algorithm>

using namespace std;

int main(int argc, char *argv[]) {
  // *******************************************************
  // *******************************************************

  int nb = 0;
  int nx = 96;
  int ny = 64;
  int nz = 48;
  int px = atoi( argv[1] );
  int py = atoi( argv[2] );
  int pz = atoi( argv[3] );
  string fileStem = "ULX";	// for the master file

  MPIClass mpi;
  mpi.startMPI(argc, argv, px, py, pz);
  int rank   = mpi.getProcID();
  int nprocs = mpi.getNumProcs();

  int irefine = 0; 		// 0=off, 1=AMR, 2=geometric z-only, 3=geometric x-y-z
  int imet_pp = 2;      // original calculation ran on a spherical grid
  int nlevels = 1; 		// number of AMR levels

  int irestart = 0;
  int nStart=624, nStop=625, nSkip=1;
  int ndims    = 1;
  if(ny > 1) ndims = 2;
  if(nz > 1) ndims = 3;  

  // physical constants
  double pi      = Constants::PI;
  double tiny    = Constants::tinyNumber;
  double msolar  = Constants::solarMass;                 // ~2e+33 gm
  double mH      = Constants::protonMass;                // ~2e-24 gm
  double gee     = Constants::gravityConstant;           // ~7e-08 cm^3/(g s^2)
  double cee     = Constants::speedOfLight;              // ~3e+10 cm/s
  double kb      = Constants::boltzmanConstant;          // ~1e-16 erg/K
  double sig     = Constants::stefanBoltzman;
  double aR      = Constants::radiationConstant;         // ~8e-15 erg/(cc K^4)
  double sigmaT  = Constants::thomsonCrossSection;       // ~6.65e-25 (cm^2)
  double hp      = Constants::planckConstant;            // ~6.63e-27 (erg s)
  double evtoerg = Constants::evToErg;                   // ~1.6e-12 erg/eV
  double yrtosec = Constants::yearToSec;                 // ~3e+07   sec/yr
  double gtocm   = gee/(cee*cee);                        // cm/g

  // -------------------------------------------------------
  // problem (BH) parameters and units used in the cosmos runs
  // should not have to change
  // -------------------------------------------------------

  string baseDirName;	// the directory where the HDF files are stored
  baseDirName = "/globalscratch/fragilep/ULX/ULX_a9r50_48x32x24L2/";
   
  double bhmass    = 6.62; 	// BH mass in solar mass units
  double bhspin    = 0.9;   	// BH spin in BH mass units (0 < a < 1)
  double bhtilt    = 0.0;    // degrees (not radians)
  double mu        = 0.615; 	// mean molecular weight
  double gamma     = 5./3.; 	// adiabatic index
  double kappa     = 0.4;     // electron scattering (cm^2/g)

  double bhgm   = bhmass*msolar; 	// BH mass in grams
  double GM     = gee*bhgm;
  double rg     = GM/(cee*cee);
  double Ledd   = 4.*pi*gee*bhgm*mH*cee/sigmaT;

  // match units from the original calculation
  double lunit  = rg;                               // length (cm)
  double tunit  = lunit/cee;                        // time (s)
  double dunit  = 1.0;                              // density (g/cm^3)
  double munit  = dunit*lunit*lunit*lunit;          // mass (g)
  double tempunit = 1.;  // Dumped data is already in Kelvin (K)

  double vunit    = lunit/tunit;
  double eunit    = munit*vunit*vunit;
  double edunit   = eunit/pow(lunit,3);

  Constants::timeUnit        = tunit;
  Constants::massUnit        = munit;
  Constants::lengthUnit      = lunit;
  Constants::temperatureUnit = tempunit;

  // convert cgs to code units
  bhmass     = GM/(pow(lunit,3)/pow(tunit,2));
  bhspin    *= bhmass;
  bhtilt    *= pi/180.0;
  rg        /= lunit;

  double bhradius = bhmass + sqrt(bhmass*bhmass - bhspin*bhspin + tiny);

  double z1  = 1.0 + pow((1. - bhspin*bhspin/(bhmass*bhmass)),1./3.)*
                    (pow((1. + bhspin/bhmass),1./3.) + pow((1. - bhspin/bhmass),1./3.));
  double z2  = sqrt(3.*bhspin*bhspin/(bhmass*bhmass) + z1*z1);
  double rmb = 2.*bhmass - bhspin + 2.*sqrt(bhmass*(bhmass - bhspin+tiny));
  double rms = bhmass*(3. + z2 - sqrt((3. - z1)*(3. + z1 + 2.*z2)));
  if(bhspin < 0.) rms =  bhmass*(3. + z2 + sqrt((3. - z1)*(3. + z1 + 2.*z2)));

   //grid parameters
   double eta0      = 1.0;
   double r1        = 0.0;
   double rI        = 1.0;
   int ilogr = 1;
   
  // min and max coordinates of original calculation (in code units)
  double r_min =  0.9*bhradius;
  double r_max =  1000.;
   if(ilogr == 1) {
      r_min = eta0 + log(r_min/bhradius);
      r_max = eta0 + log(r_max/bhradius);
   }
  double t_min =  0.0;
  double t_max =  pi;
  double p_min =  -pi;
  double p_max =  pi;

  // -------------------------------------------------------
  // create the post-processor object
  // -------------------------------------------------------

  PPCLASS pp(mpi);

  // read master file BEFORE remapping (for the fileStem if it differs from default)
  pp.readMasterFile(baseDirName, fileStem);

  // read all data, consolidating domains into a single global array
  // necessary for remapping, but optional otherwise
  pp.readAllData();

  // note: time dump interval is ~92 Rg/c (the ISCO orbital period)
  int ndimsPP    = pp.getDimension();                   // grid dimension of stored dataset
  int totalSteps = pp.getNumTimeSteps();                // number of output dump files
  double timePPmin = pp.getCurrentTime(0);              // first dump time (cycle 0), PP (likely code) units
  double timePPmax = pp.getCurrentTime(totalSteps-1);   // final dump time, PP (likely code) units

  // startTime and stopTime are used only for time binning
  double stopTime, startTime;

  stopTime = pp.getCurrentTime(nStop);
  startTime = pp.getCurrentTime(nStart);

  if(rank == 0) {
     cout << " " << endl;
     cout << "HDF file info: " << endl;
     cout << "   working directory " << baseDirName << endl;
     cout << "   data contains " << totalSteps << " cycle dumps, from t = " << timePPmin << " to " << timePPmax << " (dumped units)" <<  endl;
     cout << " " << endl;
     cout << "start and stop indices = " << nStart << ", " << nStop << endl;
     cout << "start and stop times   = " << startTime << ", " << stopTime << endl;

     cout << " " << endl;
     cout << "BH radius = " << bhradius << endl;
     cout << "marginally bound radius  = " << rmb << endl;
     cout << "marginally stable radius = " << rms << endl;
     cout << " " << endl;
     cout << "mass   unit = " << munit    <<  " gm"  << endl;
     cout << "length unit = " << lunit    <<  " cm"  << endl;
     cout << "time   unit = " << tunit    <<  " sec" << endl;
     cout << "temp   unit = " << tempunit <<  " K "  << endl;
     cout << " " << endl;
  }

  // -------------------------------------------------------
  // create the grid and solver objects
  // -------------------------------------------------------
   
  // construct new grid
  Mapper mapper(&mpi);

  // build the KLM mesh
  KLMMesh *mesh;
  if(ndims == 1) {
     mesh = new KLMMesh(mpi,nlevels,nb,nx,r_min,r_max,0);
  } else if(ndims == 2) {
     mesh = new KLMMesh(mpi,nlevels,nb,nx,ny,r_min,r_max,t_min,t_max,0);
  } else {
     mesh = new KLMMesh(mpi,nlevels,nb,nx,ny,nz,r_min,r_max,t_min,t_max,p_min,p_max,0);
  }

  int numZones = (*mesh).getNumZones();
  int domZones = (*mesh).getFirstDomainZone();
  int bndZones = (*mesh).getFirstBoundaryZone();


  // instantiate field, metric, and physics objects
  if(rank == 0) cout << "setting physics objects..." << endl;

  vector<Physics *> physics;

  Field field(*mesh, mpi);

  Metric *metric;
  SphericalBL   metricbl (*mesh, field, mpi);
  metricbl.setBHSpin(bhspin);
  metricbl.setBHMass(bhmass);
  //metricbl.setBHTilt(bhtilt);  // Pandurata grid must be in untilted frame
  metricbl.setR0(bhradius);
  metricbl.setEta0(eta0);
  metricbl.setR1(r1);
  metricbl.setRindex(rI);
  if(ilogr == 1) metricbl.turnOnLogRCoordinate();
  metric = &metricbl;
  metric->setMetric();

  EOS *eos;
  eos = new IdealGas(*mesh, field, mpi);
  (*eos).setGamma(gamma); 
  (*eos).setMoleWeight(mu); 

  // opacity sets the frequency grids for transport, M1 and MC evolve on different grids
  Opacity opacity (*mesh, field, mpi);
  Opacity opacityM1(*mesh, field, mpi);

  OpacScatter* opacs = new OpacScatter(mpi);
  opacity.addScattering("scattering", opacs);

  OpacFF opacff;
  opacity.addAcceleration(&opacff);

  HydroHRSC *hydrohrsc;
  hydrohrsc = new HydroHRSC(*mesh, field, *metric, mpi, *eos);
  hydrohrsc->turnOnMHD();
  hydrohrsc->turnOnRelativity();
  hydrohrsc->turnOnFaceMetric();

  // not sure this is actually doing anything for us here.
  Primitive *primitive;
  primitive = new Primitive(*mesh, field, *metric, mpi);
  primitive->setEOS(eos);
  //physics.push_back(primitive);

   if(rank == 0) cout << "setting boundary conditions..." << endl;

  vector<BC *> bclist;
  SymmetricBC  bc_sym(*mesh, field, *metric, mpi);
  bclist.push_back( &bc_sym );

  for(int izone = 0; izone < numZones; izone++) {
      KLMZone *zone = (*mesh).getZone(izone);
      Vector3d pos  = (*mesh).getZone(izone)->getPosition();

      double x = pos.getX();
      double y = pos.getY();
      double z = pos.getZ();

      if(x > r_max) bc_sym.addZone(zone);
      if(y > t_max) bc_sym.addZone(zone);
      if(z > p_max) bc_sym.addZone(zone);

      if(x < r_min) bc_sym.addZone(zone);
      if(y < t_min) bc_sym.addZone(zone);
      if(z < p_min) bc_sym.addZone(zone);
  }

  RefineMesh   refine(*mesh, field, *metric, mpi);
  refine.setVerbose(0);

  if(irefine == 1) {
     refine.turnOnRefinement();
     refine.turnOnProperFieldTag();
     refine.turnOnConservativeInterpolation();   // default
  } else {
     refine.turnOffRefinement();
  }

  // set boundary conditions
  for(int i = 0; i < bclist.size(); i++) {
      refine.addBC( bclist[i] );
      for(int j = 0; j < physics.size(); j++) {
          physics[j]->addBC( bclist[i] );
      }
  }

  // -------------------------------------------------------
  // output dumps
  // -------------------------------------------------------
  if(rank == 0) cout << "setting dump objects..." << endl;

  // map of fields and associated multipliers needed to scale fields to code units
  // suffix identifies scalar fields ("_s") from vectors ("_v")
  map<string, double> fmap;
   fmap = {{"rho_s", 1.0}, {"temperature_s", 1.0}, {"velocity_v", 1.0}, {"magneticPressure_s", 1.0}};

  vector<int> izone_pp(domZones, -99999);
  vector<double> zoneCenter_pp(ndims*domZones, 0.);
  // loop over the cosmos data sets (different time dumps)
  for(int ncyc = nStart; ncyc < nStop; ncyc+=nSkip) {
     double time = pp.getCurrentTime(ncyc);   // use code units here
     pp.readData(ncyc);
     pp.getZoneCenter(zoneCenter_pp);

     if(rank == 0) cout << "\nstarting cycle " << ncyc << " of " << nStop << ", time = " << time << endl;

     // perform mapping
     if(rank == 0) cout << "starting remap cycle " << ncyc << endl;

     FILE *OUT;
     string fileName = "grid_map.dat";
     if(irestart == 0 && ncyc == nStart) {
        mapper.mapGrid(&pp, mesh, &field, &mpi, metric, ncyc, bhtilt, izone_pp);
        // save grid map to file
        OUT = fopen(fileName.data(), "w");
        for(int izone = 0; izone < domZones; izone++) fprintf(OUT,"%i \n", izone_pp[izone]);
        fclose(OUT);
     } else if(irestart == 1 && ncyc == nStart) {
        // read map from file
        OUT = fopen(fileName.data(), "r");
        for(int izone = 0; izone < domZones; izone++) fscanf(OUT,"%i", &izone_pp[izone]);
        fclose(OUT);
     }

     vector<double> magpres(numZones, 0.), temperature(numZones, 0.), rho(numZones, 0.), u0pKS(numZones, 1.);
     mapper.mapFields(&pp, mesh, &field, &mpi, metric, ncyc, imet_pp, fmap, izone_pp, u0pKS);

      if(rank == 0) cout << "completed remap cycle " << ncyc << endl;
     
     vector<double> sigtau_es(nx*nz, 0.), tau_es(nx*(ny+1)*nz, 0.);
     vector<double> diskbody(nx*nz, 0.), diskflux(nx*nz, 0.), Tdisk(nx*nz, 0.);
     vector<double> em_top(nx*nz, 0.5*pi), em_bot(nx*nz, 0.5*pi), ref_top(nx*nz, 0.5*pi), ref_bot(nx*nz, 0.5*pi);
     vector<Vector3d> velocity(numZones, Vector3d(0.,0.,0.));
     
     double corona_tau = 1.;
     double photo_tau  = corona_tau + 1.;

     fileName = "dump_times.dat";
     if(ncyc == nStart) OUT = fopen(fileName.data(), "w");
     else               OUT = fopen(fileName.data(), "a");
     if(OUT == NULL) {
        printf("problem opening output file in PostProcess\n");
        exit(0);
     }
     fprintf(OUT,"\t%12.5e \n", time);
     fclose(OUT);

     char str[5];
     snprintf (str, 5, "%04d", ncyc);
     fileName = "gr_"+string(str)+"_new.dat";
     OUT = fopen(fileName.data(), "w");
     fprintf(OUT,"%d %d %d \n", nx, ny, nz);
     for(int i = 0; i < nx; i++) {
        Vector3d pos = (*mesh).getZone(i)->getPosition();
        (*metric).getPhysicalPosition(pos);
        fprintf(OUT,"\t%12.5e", pos.getX());
     }
     fprintf(OUT,"\n");
     for(int j = 0; j < ny; j++) {
        Vector3d pos = (*mesh).getZone(j*nx)->getPosition();
        (*metric).getPhysicalPosition(pos);
        fprintf(OUT,"\t%12.5e", pos.getY());
     }
     fprintf(OUT,"\n");
     for(int k = 0; k < nz; k++) {
        Vector3d pos = (*mesh).getZone(k*nx*ny)->getPosition();
        (*metric).getPhysicalPosition(pos);
        fprintf(OUT,"\t%12.5e", pos.getZ()+pi);
     }
     fprintf(OUT,"\n");
     fclose(OUT);

     fileName = "rh_"+string(str)+"_new.dat";
     OUT = fopen(fileName.data(), "w");
     field.getScalarField("rho", rho);
     for(int k = 0; k < nz; k++) {
        for(int j = 0; j < ny; j++) {
           for(int i = 0; i < nx; i++) {
              fprintf(OUT,"\t%12.5e", rho[k*nx*ny+j*nx+i]*dunit);
           }
           fprintf(OUT,"\n");
        }
     }
     fclose(OUT);

     fileName = "te_"+string(str)+"_new.dat";
     OUT = fopen(fileName.data(), "w");
     field.getScalarField("temperature", temperature);
     for(int k = 0; k < nz; k++) {
        for(int j = 0; j < ny; j++) {
           for(int i = 0; i < nx; i++) {
              fprintf(OUT,"\t%12.5e", temperature[k*nx*ny+j*nx+i]*tempunit);
           }
           fprintf(OUT,"\n");
        }
     }
     fclose(OUT);

     // TODO: Check units
     fileName = "bb_"+string(str)+"_new.dat";
     OUT = fopen(fileName.data(), "w");
     field.getScalarField("magneticPressure", magpres);
     for(int k = 0; k < nz; k++) {
        for(int j = 0; j < ny; j++) {
           for(int i = 0; i < nx; i++) {
              fprintf(OUT,"\t%12.5e", magpres[k*nx*ny+j*nx+i]);
           }
           fprintf(OUT,"\n");
        }
     }
     fclose(OUT);

     FILE *OUT0, *OUT1, *OUT2, *OUT3;
     string fileName0 = "u0_"+string(str)+"_new.dat";
     string fileName1 = "u1_"+string(str)+"_new.dat";
     string fileName2 = "u2_"+string(str)+"_new.dat";
     string fileName3 = "u3_"+string(str)+"_new.dat";
     OUT0 = fopen(fileName0.data(), "w");
     OUT1 = fopen(fileName1.data(), "w");
     OUT2 = fopen(fileName2.data(), "w");
     OUT3 = fopen(fileName3.data(), "w");
     field.getVectorField("velocity", velocity);
     for(int k = 0; k < nz; k++) {
        for(int j = 0; j < ny; j++) {
           for(int i = 0; i < nx; i++) {
              int izone = k*nx*ny+j*nx+i;
              // Tilted KS coords
              double xpKS = zoneCenter_pp[3*izone_pp[izone]+0];
              double ypKS = zoneCenter_pp[3*izone_pp[izone]+1];
              double zpKS = zoneCenter_pp[3*izone_pp[izone]+2];
              double rpKS = sqrt(xpKS*xpKS + ypKS*ypKS + zpKS*zpKS);
              double sintp = sqrt(xpKS*xpKS + ypKS*ypKS)/rpKS;
              double costp = zpKS/rpKS;
              double sinpp = ypKS/sqrt(xpKS*xpKS + ypKS*ypKS);
              double cospp = xpKS/sqrt(xpKS*xpKS + ypKS*ypKS);
              double urpKS = u0pKS[izone]*velocity[izone].getX();
              double utpKS = u0pKS[izone]*velocity[izone].getY();
              double uppKS = u0pKS[izone]*velocity[izone].getZ();
              // Untilted KS coords
              double xKS = cos(bhtilt)*xpKS + sin(bhtilt)*zpKS;
              double yKS = ypKS;
              double zKS = -sin(bhtilt)*xpKS + cos(bhtilt)*zpKS;
              double rKS = sqrt(xKS*xKS + yKS*yKS + zKS*zKS);
              double sint = sqrt(xKS*xKS + yKS*yKS)/rKS;
              double cost = zKS/rKS;
              double cosp = xKS/sqrt(xKS*xKS + yKS*yKS);
              double dcostdtp = sin(bhtilt)*costp*cospp - cos(bhtilt)*sintp;
              double dcostdpp = -sin(bhtilt)*sintp*sinpp;
              double dtdtp = -1./sint*dcostdtp;
              double dtdpp = -1./sint*dcostdpp;
              double dsinpdtp = costp*sinpp/sint - sintp*sinpp*cost/(sint*sint)*dtdtp;
              double dsinpdpp = -sintp*cospp/sint - sintp*sinpp*cost/(sint*sint)*dtdpp;
              double dpdtp = 1./cosp*dsinpdtp;
              double dpdpp = 1./cosp*dsinpdpp;
              double urKS = urpKS;
              double utKS = dtdtp*utpKS + dtdpp*uppKS;
              double upKS = dpdtp*utpKS + dpdpp*uppKS;
              // Untilted BL coords
              Vector3d pos = (*mesh).getZone(izone)->getPosition();
              (*metric).getPhysicalPosition(pos);
              double r = pos.getX();
              double delta = r*r - 2.*bhmass*r + bhspin*bhspin;
              double upBL = upKS - bhspin/delta*urKS;
              Tensor4d met = (*mesh).getZone(izone)->getMetric();
              double g00 = met.getTT();
              double g01 = met.getTX();
              double g02 = met.getTY();
              double g03 = met.getTZ();
              double g11 = met.getXX();
              double g12 = met.getXY();
              double g13 = met.getXZ();
              double g22 = met.getYY();
              double g23 = met.getYZ();
              double g33 = met.getZZ();
              Tensor4d metinv = met.getInverse();
              double g00up = metinv.getTT();
              if(ilogr == 1) {
                 g01 /= r;
                 g11 /= (r*r);
                 g12 /= r;
                 g13 /= r;
              }
              double alpha = 1./sqrt(-g00up);
              double gam2 = (1. +g11*urKS*urKS + g22*utKS*utKS + g33*upBL*upBL +
                      2.*(g12*urKS*utKS + g13*urKS*upBL + g23*utKS*upBL));
              double gam  = sqrt(gam2);
              double u0BL = gam/alpha;

              fprintf(OUT0,"\t%12.5e", u0BL);
              fprintf(OUT1,"\t%12.5e", urKS);
              fprintf(OUT2,"\t%12.5e", utKS);
              fprintf(OUT3,"\t%12.5e", upBL);
              
              if(r > bhradius) {
                 if(!isfinite(u0BL)) {
                    cout << "u0 is nan" << endl;
                    cout << "   izone    = " << izone << "  izone pp = " << izone_pp[izone] << endl;
                    cout << "   position = " << pos << ", " <<
                    Vector3d(zoneCenter_pp[3*izone_pp[izone]+0],
                             zoneCenter_pp[3*izone_pp[izone]+1],
                             zoneCenter_pp[3*izone_pp[izone]+2]) << endl;
                 }
                 pos = (*mesh).getZone(izone)->getPosition();
                 (*metric).getCartesianPosition(pos);
                 Vector4d vtest = Vector4d(u0BL, urKS, utKS, upBL);
                 Tensor4d met   = Tensor4d(g00, g01, g02, g03, g11, g12, g13, g22, g23, g33);
                 double v2 = vtest.dot(met.dot(vtest));
                 
                 if(fabs(1.+v2) > 1.e-2) {
                    cout << "four-velocity is not properly normalized:" << endl;
                    cout << "   izone    = " << izone << "  izone pp = " << izone_pp[izone] << endl;
                    cout << "   position = " << pos << ", " <<
                    Vector3d(zoneCenter_pp[3*izone_pp[izone]+0],
                             zoneCenter_pp[3*izone_pp[izone]+1],
                             zoneCenter_pp[3*izone_pp[izone]+2]) << endl;
                    cout << "   v2 = " << v2 << endl;
                    cout << "   u0up = " << u0pKS[izone] << " " << u0BL << endl;
                    cout << "   velocity = " << velocity[izone] << endl;
                 }
              }
              
           }
           fprintf(OUT0,"\n");
           fprintf(OUT1,"\n");
           fprintf(OUT2,"\n");
           fprintf(OUT3,"\n");
        }
     }
     fclose(OUT0);
     fclose(OUT1);
     fclose(OUT2);
     fclose(OUT3);

     // Find photosphere
     fileName = "ta_"+string(str)+"_new.dat";
     OUT = fopen(fileName.data(), "w");
     for(int k = 0; k < nz; k++) {
        for(int i = 0; i < nx; i++) {
           vector<double> tauE_b(ny+1,0.0), tauE_t(ny+1,0.0);
           double thm1 = 0.;
           for(int j = 1; j <= ny; j++) {
              int izone = k*nx*ny+(j-1)*nx+i;
              Vector3d pos = (*mesh).getZone(izone)->getPosition();
              (*metric).getPhysicalPosition(pos);
              //double dth = pos.getY() - thm1;
              double dth = (*mesh).getZone(izone)->getVectorZoneLength()[1];
              tauE_t[j] = tauE_t[j-1] + kappa*rho[k*nx*ny+(j-1)*nx+i]*dunit*pos.getX()*lunit*dth;
              //if(i == 24 && k == 0) cout << j << " " << pos.getY() << " " << dth << " " << tauE_t[j] << endl;
              thm1 = pos.getY();
           }
           thm1 = pi;
           for(int j = ny-1; j >= 0; j--) {
              int izone = k*nx*ny+j*nx+i;
              Vector3d pos = (*mesh).getZone(izone)->getPosition();
              (*metric).getPhysicalPosition(pos);
              //double dth = pos.getY() - thm1;
              double dth = (*mesh).getZone(izone)->getVectorZoneLength()[1];
              tauE_b[j] = tauE_b[j+1] + kappa*rho[k*nx*ny+j*nx+i]*dunit*pos.getX()*lunit*dth;
              thm1 = pos.getY();
           }
           sigtau_es[k*nx+i] = tauE_t[ny];
           for(int j = 0; j <= ny; j++) {
              tau_es[k*nx*(ny+1)+j*nx+i] = min(tauE_b[j], tauE_t[j]);
           }
        }
     }
     for(int k = 0; k < nz; k++) {
        for(int j = 0; j <= ny; j++) {
           for(int i = 0; i < nx; i++) {
              fprintf(OUT,"\t%12.5e", tau_es[k*nx*(ny+1)+j*nx+i]);
           }
           fprintf(OUT,"\n");
        }
     }
     fclose(OUT);

     fileName = "ph_"+string(str)+"_new.dat";
     OUT = fopen(fileName.data(), "w");
     for(int k = 0; k < nz; k++) {
        for(int i = 0; i < nx; i++) {
           if(sigtau_es[k*nx+i] > 2.*photo_tau) {
              diskbody[k*nx+i] = 2.;
              for(int j = 1; j <= ny; j++) {
                 Vector3d pos = (*mesh).getZone(k*nx*ny+(j-1)*nx+i)->getPosition();
                 (*metric).getPhysicalPosition(pos);
                 //cout << j << " " << pos.getY() << " " << photo_tau << " " << tau_es[k*nx*ny+j*nx+i] << endl;
                 if(tau_es[k*nx*(ny+1)+j*nx+i] > photo_tau) {
                    ref_top[k*nx+i] = pos.getY();
                    break;
                 }
              }
              for(int j = ny-1; j >= 0; j--) {
                 Vector3d pos = (*mesh).getZone(k*nx*ny+j*nx+i)->getPosition();
                 (*metric).getPhysicalPosition(pos);
                 if(tau_es[k*nx*(ny+1)+j*nx+i] > photo_tau) {
                    ref_bot[k*nx+i] = pos.getY();
                    break;
                 }
              }
           }
           if(sigtau_es[k*nx+i] > 2.*corona_tau) {
              if(sigtau_es[k*nx+i] <= 2.*photo_tau) diskbody[k*nx+i] = 1.;
              //Tdisk[k*nx+i] = pow(0.5*diskflux[k*nx+i]/sig,0.25);
              double T_top = 1.e6, T_bot = 1.e6;
              for(int j = 1; j <= ny; j++) {
                 Vector3d pos = (*mesh).getZone(k*nx*ny+(j-1)*nx+i)->getPosition();
                 (*metric).getPhysicalPosition(pos);
                 //cout << "yes " << k << " " << i << " " << j << " " << pos.getY() << " " << tau_es[k*nx*(ny+1)+j*nx+i] << " " << corona_tau << endl;
                 if(tau_es[k*nx*(ny+1)+j*nx+i] > corona_tau) {
                    em_top[k*nx+i] = pos.getY();
                    // ion temperature
                    T_top = temperature[k*nx*ny+j*nx+i]*tempunit;
                    // temperature assuming radiation pressure dominated
                    double n_top = rho[k*nx*ny+j*nx+i]*dunit/mu/mH;
                    T_top = pow(3.*cee/4./sig*kb*n_top*T_top,0.25);
                    break;
                 }
              }
              for(int j = ny-1; j >= 0; j--) {
                 Vector3d pos = (*mesh).getZone(k*nx*ny+j*nx+i)->getPosition();
                 (*metric).getPhysicalPosition(pos);
                 if(tau_es[k*nx*(ny+1)+j*nx+i] > corona_tau) {
                    em_bot[k*nx+i] = pos.getY();
                    // ion temperature
                    T_bot = temperature[k*nx*ny+j*nx+i]*tempunit;
                    // temperature assuming radiation pressure dominated
                    double n_bot = rho[k*nx*ny+j*nx+i]*dunit/mu/mH;
                    T_bot = pow(3.*cee/4./sig*kb*n_bot*T_bot,0.25);
                    break;
                 }
              }
              Tdisk[k*nx+i] = 0.5*(T_top + T_bot);
           }
           fprintf(OUT,"\t%12.5e", diskbody[k*nx+i]);
        }
        fprintf(OUT,"\n");
     }
     for(int k = 0; k < nz; k++) {
        for(int i = 0; i < nx; i++) {
           fprintf(OUT,"\t%12.5e", sigtau_es[k*nx+i]);
        }
        fprintf(OUT,"\n");
     }
     for(int k = 0; k < nz; k++) {
        for(int i = 0; i < nx; i++) {
           fprintf(OUT,"\t%12.5e", Tdisk[k*nx+i]);
        }
        fprintf(OUT,"\n");
     }
     for(int k = 0; k < nz; k++) {
        for(int i = 0; i < nx; i++) {
           fprintf(OUT,"\t%12.5e", em_top[k*nx+i]);
        }
        fprintf(OUT,"\n");
     }
     for(int k = 0; k < nz; k++) {
        for(int i = 0; i < nx; i++) {
           fprintf(OUT,"\t%12.5e", em_bot[k*nx+i]);
        }
        fprintf(OUT,"\n");
     }
     for(int k = 0; k < nz; k++) {
        for(int i = 0; i < nx; i++) {
           fprintf(OUT,"\t%12.5e", ref_top[k*nx+i]);
        }
        fprintf(OUT,"\n");
     }
     for(int k = 0; k < nz; k++) {
        for(int i = 0; i < nx; i++) {
           fprintf(OUT,"\t%12.5e", ref_bot[k*nx+i]);
        }
        fprintf(OUT,"\n");
     }
     fclose(OUT);

  } // end loop over HDF data sets (n=nStart,nStop)


  // -------------------------------------------------------
  // wrap-up
  // -------------------------------------------------------
  mpi.Barrier(); // make sure the file writing finished
  TimeHydro.timerSummary(rank, nprocs);

  mpi.finalize();
  return 0;
}
