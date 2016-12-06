#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <assert.h>
#include <algorithm>
#include <fftw3.h>
#include <complex>
#include <omp.h>

using namespace std;

#define PI 3.14159

int grid = 1024; //Number of grids on one 
double length = 8000; //Units of h^-1 Mpc
double cell = length/grid; //Size of a cell in units of h^-1 Mpc, what's usually called $\Delta$
double cell3 = pow(cell,3.0);
int nbin = 30;

//Define structure of input galaxy and random files:
//(x,y,z) coordinates + weights + n(r)
struct sdss
{
 double x, y, z, w, n;
};

//Function which the initializes the FFTd k vector according to the usual standard, e.g. Fig. 12.2.2b of Numerical Recipes
double* initialize_k(double k0) //k0 = 2\pi/length is the basic wavenumber = 2\pi basic frequency
{
   int i25;
   double *myarray;
   
   myarray = new double[grid];
   
   for(i25 = 0; i25 < grid; i25++)
   {
      if(i25 < grid/2)
      {
         myarray[i25] = i25*k0;
      }
      else
      {
         myarray[i25] = (i25 - grid)*k0;
      }
   }

   return myarray;
}

//Function which returns i^2 if i < grid/2 and (grid - i)^2 if i < grid/2
double ksq(int k1)
{
   if(k1 < grid/2)
   {
      return k1*k1;
   }
   else
   {
      return (grid - k1)*(grid-k1);
   }
}

//Template used to initialize 1D arrays (vectors) of any type
template <class IV>
IV* init_1D_array(long size)
{
   long j1;
   
   IV *vec;
   vec = new IV[size];
   
   for(j1 = 0; j1 < size; j1++)
   {
      vec[j1] = 0.0; //Fill in the vector with zeros
   }
   
   return vec;
}

//template used to initialize 3D arrays of any type
template <class IA>
IA*** init_3D_array(long size1, long size2, long size3)
{
   long j2, j3, j4;

   IA ***arr;
   arr = new IA**[size1];

   for(j2 = 0; j2 < size1; j2++)
   {
      arr[j2] = new IA*[size2];
      
      for(j3 = 0; j3 < size2; j3++)
      {
         arr[j2][j3] = new IA[size3];
         
         for(j4 = 0; j4 < size3; j4++)
         {
            arr[j2][j3][j4] = 0.0; //Fill in the array with zeros
         }
      }
   }

   return arr;
}

//Template used to kill 3D arrays of any type, in order to free memory
template <class KA>
KA* kill_3D_array(KA ***karr, long sizex, long sizey)
{
   long j5, j6;
   
   for(j5 = 0; j5 < sizex; j5++)
   {
      for(j6 = 0; j6 < sizey; j6++)
      {
         delete[] karr[j5][j6];
      }
      delete[] karr[j5];
   }
   delete[] karr;
   
   return 0;
}

double sinc(double x)
{
   //If x << 1 sin(x)/x~1, want to avoid divergency for x = 0
   if(fabs(x) < 1e-8)
   {
      return 1;
   }
   else
   {
    return sin(x)/x;
   }
}

//Implements cloud-in-cell (CIC) mass assignment scheme (MAS) to galaxy and random catalogs
//"particles" is galaxy or random catalog, N is, length is length of box
int cic(sdss *particles, long N, double ***nr)
{
   long i1;

   for(i1 = 0; i1 < N; i1++)
   {
      //xp, yp, zp are coordinates of particles, divide by "cell" so to move to next cell just add or subtract 1
      double xp = particles[i1].x/cell;
      double yp = particles[i1].y/cell;
      double zp = particles[i1].z/cell;
      long xc = floor(xp); //Coordinates of central cell containing particle
      long yc = floor(yp);
      long zc = floor(zp);
      double dx = xp - xc;
      double dy = yp - yc;
      double dz = zp - zc;
      double tx = 1.0 - dx;
      double ty = 1.0 - dy;
      double tz = 1.0 - dz;

      //Impose periodic boundary conditions
      long xplus = xc + 1;
      long yplus = yc + 1;
      long zplus = zc + 1;
      if(xplus == grid)
      {
         xplus = 0;
      }
      if(yplus == grid)
      {
         yplus = 0;
      }
      if(zplus == grid)
      {
         zplus = 0;
      }
      
      //Fill eight cells around particle and divide by cell^3 to get number density
      nr[xc][yc][zc] = nr[xc][yc][zc] + particles[i1].w*tx*ty*tz/cell3;
      nr[xplus][yc][zc] = nr[xplus][yc][zc] + particles[i1].w*dx*ty*tz/cell3;
      nr[xc][yplus][zc] = nr[xc][yplus][zc] + particles[i1].w*tx*dy*tz/cell3;
      nr[xc][yc][zplus] = nr[xc][yc][zplus] + particles[i1].w*tx*ty*dz/cell3;
      nr[xplus][yplus][zc] = nr[xplus][yplus][zc] + particles[i1].w*dx*dy*tz/cell3;
      nr[xplus][yc][zplus] = nr[xplus][yc][zplus] + particles[i1].w*dx*ty*dz/cell3;
      nr[xc][yplus][zplus] = nr[xc][yplus][zplus] + particles[i1].w*tx*dy*dz/cell3;
      nr[xplus][yplus][zplus] = nr[xplus][yplus][zplus] + particles[i1].w*dx*dy*dz/cell3;
   }
 
  return 0;
}

//Function that returnes number of lines of a galaxy or random catalog, required since we will be allocating memory for them dynamically
long lines(string linesfile)
{  
   long numlines = 0;
   string linesstring;
   ifstream linesstream;

   linesstream.open( linesfile.c_str(), ios::in );
 
   if(linesstream.is_open())
   {
      while(getline(linesstream,linesstring), !linesstream.eof()) //Need to put getline(mystream, line)!! Otherwise with the eof part only it will read one more line before stopping the while loop
      {
         numlines++;
      }
     linesstream.close();
   }

   return numlines;
}

//Function which reads galaxy and random files
sdss* reader(string readername, long *num)
{
 ifstream readerstream;
 sdss *catalog;
 long i2;

 //Print what it is reading (galaxy or random) and the number of objects it contains
 cout << "Reading " << readername << endl;
 *num = lines(readername);
 cout << readername << " contains " << *num << " objects" << endl;

 catalog = new sdss[*num];
 readerstream.open( readername.c_str(), ios::in );
 
 for(i2 = 0; i2 < *num; i2++)
 {
    string str;
    
    getline(readerstream, str);
    stringstream stream(str);
    stream >> catalog[i2].x >> catalog[i2].y >> catalog[i2].z >> catalog[i2].w >> catalog[i2].n;
 }

 return catalog;
}

//Function to bin in k
int binner(double ***pkgrid, double *kavg, double *pavg, double normalization, double shotnoise)
{  
   //To decide in which bin does a certain k go, compare log10 of its magnitude with width defined below
   //Note we are working in units of the fundamental wavevector k0 = 2\pi/length
   double width = log10(grid/2.0)/nbin;
   double k2, logk;
   double *count = init_1D_array<double>(nbin) ;
   double *totlogk = init_1D_array<double>(nbin);
   double *totkpk = init_1D_array<double>(nbin);

   
   
   int i9, i10, i11, i12, bin;
   
   cout << "Binning in progress..." << endl;

   for(i9 = 0; i9 < grid; i9++)
   {  
      for(i10 = 0; i10 < grid; i10++)
      {
         for(i11 = 0; i11 < grid; i11++)
         {
            if(i9 == 0 && i10 == 0 && i11 == 0)
            {
               continue;
            }
                    
            k2 = ksq(i9) + ksq(i10) + ksq(i11); //Magnitude squared of k-vector in units of fundamental wavevector (2\pi/length)
            logk = 0.5*log10(k2); //log10 of the magnitude of k (in the units above), note the 0.5 to take the square root of k2
            bin = floor(logk/width); //Number of bin in which the k-vector will go, between 0 and nbin;
            
            if(bin >= 0 && bin < nbin)
            {
               count[bin] = count[bin] + 1; //Increase number of objects in the given bin
               totlogk[bin] = totlogk[bin] + logk; //Sum log(k) within the bin (need later to find average k and average P)
               totkpk[bin] = totkpk[bin] + (sqrt(k2)*pkgrid[i9][i10][i11]); //Sum k*P(k) within the bin (need later to find average k and average P)
            }
         }
      }
   }

 //For each bin find average k and average P if the number of objects in the bin is nonzero 
 for(i12 = 0; i12 < nbin; i12++)
 {
   if(count[i12] > 0)
   {
      kavg[i12] = pow(10.0,totlogk[i12]/count[i12]); //Average k is binned k
      pavg[i12] = totkpk[i12]/count[i12]; //Average product of k and P(k) in the bin
      pavg[i12] = pavg[i12]/kavg[i12]; //Divide by binned k to get average P(k) in the bin, like doing a spherical average
      pavg[i12] = pavg[i12]/normalization; //Divide by FKP normalization
      pavg[i12] = pavg[i12] - shotnoise; //Subtract shot noise, this is final estimate of galaxy power spectrum convolved with survey window function in the bin
      kavg[i12] = kavg[i12]*(2.0*PI/length); //Binned k in units of h Mpc^-1
   }
 }
 delete[] count;
 delete[] totlogk;
 delete[] totkpk;

 cout << "Binning completed successfully with " << nbin << " bins" << endl;

 return 0;
}
   

//Function which takes 3 indices and returns position row-major ordered FFT'd array
long rowmajor(int j1, int j2, int j3, long lastsize)
{
   long pos = (j1*grid + j2)*lastsize + j3;
   return pos;
}

//Function which computes FT of weighted overdensity and corrects for MAS
int ps(double ***Fr, double normalization, double shotnoise, double *kavg, double *pavg)
{     
   int numthreads = 1;
   int i3, i4, i5, i6, i7, i8, counter;
   long complexz = grid/2 + 1; //Size of FFT'd array in the last variable (z) because of how r2c transform works
   double *fkp1d = init_1D_array<double>(grid*grid*(grid+2)); //Contains weighted overdensity organized in a 1D vector for FFT, will be overwritten by FT of weighted overdensity in 1D
   double ***power;
   double *kgrid;
   double k0 = 2.0*PI/length;

   kgrid = initialize_k(k0);
   
   counter = 0;
   //Fill real array, requires padding (adding 2 additional spots for each value of the 3rd dimension) for in-place transforms
   for(i3 = 0; i3 < grid; i3++)
   {    
      for(i4 = 0; i4 < grid; i4++)
      {       
         for(i5 = 0; i5 < grid; i5++)
         {
            fkp1d[counter] = Fr[i3][i4][i5];
            counter = counter + 1;
         }
         counter = counter + 2; //For padding due to FFTW geometry
      }
   }
    
   //Call fftw_init_threads to perform one-time initialization to use threads, returns 0 if error
   int fftwinit;
   fftwinit = fftw_init_threads();
    
   if(fftwinit == 0)
   {    
      cout << "FFT initialized unsuccessfully" << endl;
   }
   else
   {
      cout << "FFT initialized successfully" << endl;
   }

   //Call fftw_plan_with_nthreads before creating parallelizable plan with nthreads threads
   fftw_plan_with_nthreads(1);

   //Define FFTW plan: dft (discrete FT), r2c (real to complex), 3d, with real data of dimension grid^3
   //The flag FFTW_ESTIMATE just builds what it thinks is a reasonable plan which could be suboptimal
   //If initialization time is not important, use the flag FFTW_MEASURE instead
   fftw_plan plan;
   plan = fftw_plan_dft_r2c_3d(grid, grid, grid, fkp1d, (fftw_complex *)fkp1d, FFTW_ESTIMATE); //FFT is done in-place, input is overwritten by output

   //Execute FFTW plan defined before
   fftw_execute(plan);

   //Deallocate the plan once we are done with it
   fftw_destroy_plan(plan);

   //Get rid of memory and resources internally allocated by FFTW
   //Also gets rid of threads data
   //DO NOT execute any previously created plan after calling this function
   fftw_cleanup_threads();

   //Done FFT, now correct for MAS dividing FT of overdensity by MAS window function
   //CIC MAS window function is sinc[pi*k/(2*k_Ny)] \equiv sinc[k/k_CIC]
   //k_CIC \equiv 2*k_Ny/pi, where k_Ny = \pi/cell is the Nyquist frequency, so k_CIC = 2/cell
   double kcic = 2.0/cell;

   power = init_3D_array<double>(grid, grid, grid);

   double sincx, sincy, sincz; //sinc[k/k_CIC,x_i]
   double sinc2x, sinc2y, sinc2z, Wcic; //sinc^2[k/k_CIC,x_i] and W_CIC = \Pi_i=1^3 sinc^2[k/k_CIC,x_i]

//   #pragma omp parallel for schedule(dynamic)
   for(i6 = 0; i6 < grid; i6++)
   {
      sincx = sinc(kgrid[i6]/kcic);
      sinc2x = sincx*sincx;
      for(i7 = 0; i7 < grid; i7++)
      {
         sincy = sinc(kgrid[i7]/kcic);
         sinc2y = sincy*sincy;
         for(i8 = 0; i8 < complexz; i8++)
         {   
            sincz = sinc(kgrid[i8]/kcic);
            sinc2z = sincz*sincz;
             
            Wcic = sinc2x*sinc2y*sinc2z;
             
            //Retrieve position of element in the FFT'd array, stored in row-major order
            long fftpos = rowmajor(i6,i7,i8,complexz);
            long refftpos = 2*fftpos; //Position of real part;
            long imfftpos = refftpos + 1; //Position of imaginary part;

            //Divide FT of overdensity field by volume of cell (due to normalization conventions) and by MAS window function to correct for gridding
            fkp1d[refftpos] = fkp1d[refftpos]*cell3/Wcic;
            fkp1d[imfftpos] = fkp1d[imfftpos]*cell3/Wcic;

            //Take modulus squared of FT of overdensity field
            power[i6][i7][i8] = fkp1d[refftpos]*fkp1d[refftpos];
            power[i6][i7][i8] = power[i6][i7][i8] + fkp1d[imfftpos]*fkp1d[imfftpos];

            //Because of symmetry fkp1d[grid-i8] = (fkp1d[i8])*, special cases are i8 = 0 and i8 = grid/2, but note that grid - grid/2 = grid/2
            if(i8!=0)
            {
               power[i6][i7][grid-i8] = power[i6][i7][i8];
            }
         }   
      }
   }   

   delete[] fkp1d;
   delete[] kgrid;

   cout << "FFT completed successfully" << endl;
   binner(power, kavg, pavg, normalization, shotnoise);

   kill_3D_array(power, grid, grid);

   return 0;
}

int main(int argc, char* argv[])
{
   string galiname, raniname, galext, ranext, myout, myoutext;
   sdss *gal, *ran;
   long numgal, numran, i13, i14, i15, i16, i17, i18, i19, i20, i21, i22, i23, i24;
   double xming, yming, zming, xminr, yminr, zminr, xmin, ymin, zmin;
   double dngwg, dnrwr, alpha, di19, dnw2, di20, di20sq, dn2w2;
   double ***ngal, ***nran, ***fkp, *kfinal, *pfinal;
   int imain = 1;

   double normalization = 1.0;
   double shotnoise = 1.0;

//   galiname = "galaxy"; //Name of galaxy catalog without extension
//   raniname = "random"; //Name of random catalog without extension

   while(imain < argc)
   {
      string arg = argv[imain];
      if(arg == "-gal")
      {
        imain++;
        galiname = argv[imain];
        imain++;
      }
      if(arg == "-ran")
      {
        imain++;
        raniname = argv[imain];
        imain++;
      }
      if(arg == "-galext")
      {
        imain++;
        galext = argv[imain];
        imain++;
      }
      if(arg == "-ranext")
      {
        imain++;
        ranext = argv[imain];
        imain++;
      }
      if(arg == "-out")
      {
        imain++;
        myout = argv[imain];
        imain++;
      }
      if(arg == "-outext")
      {
        imain++;
        myoutext = argv[imain];
        imain++;
      }
  }

   string galifile = "../Outputs/" + galiname + "." + galext;
   string ranifile = "../Outputs/" + raniname + "." + ranext;

   gal = reader(galifile, &numgal);
   ran = reader(ranifile, &numran);

   xming = 1.0e33;
   yming = 1.0e33;
   zming = 1.0e33;
   xminr = 1.0e33;
   yminr = 1.0e33;
   zminr = 1.0e33;

   //Find minimum coordinates of galaxy catalog
   for(i13 = 0; i13 < numgal; i13++)
   {
      if(gal[i13].x < xming)
      {
         xming = gal[i13].x;
      }
      if(gal[i13].y < yming)
      {
         yming = gal[i13].y;
      }
      if(gal[i13].z < zming)
      {
         zming = gal[i13].z;
      }
   }

   //Find minimum coordinates of random catalog
   for(i14 = 0; i14 < numran; i14++)
   {
      if(ran[i14].x < xminr)
      {
         xminr = ran[i14].x;
      }
      if(ran[i14].y < yminr)
      {
         yminr = ran[i14].y;
      }
      if(ran[i14].z < zminr)
      {
         zminr = ran[i14].z;
      }
   }

   if(xming <= xminr)
   {
      xmin = xming;
   }
   else
   {
      xmin = xminr;
   }
   
   if(yming <= yminr)
   {
      ymin = yming;
   }
   else
   {
      ymin = yminr;
   }
   
   if(zming <= zminr)
   {
      zmin = zming;
   }
   else
   {
      zmin = zminr;
   }

   xmin = xmin - 200.0;
   ymin = ymin - 200.0;
   zmin = zmin - 200.0;

   //Subtract minimum coordinates to galaxy and random catalogues
   for(i15 = 0; i15 < numgal; i15++)
   {
      gal[i15].x = gal[i15].x - xmin;
      gal[i15].y = gal[i15].y - ymin;
      gal[i15].z = gal[i15].z - zmin;
   }
   
   for(i16 = 0; i16 < numran; i16++)
   {
      ran[i16].x = ran[i16].x - xmin;
      ran[i16].y = ran[i16].y - ymin;
      ran[i16].z = ran[i16].z - zmin;
   }
   
   //Required to calculate alpha = \int n_g(r)w_g(r)/\int n_r(r)w_r(r)
   double ngwg = 0.0;
   double nrwr = 0.0;
   
   //Required to calculate shot noise and FKP normalization
   double nw2 = 0.0;
   double n2w2 = 0.0;

   //Calculate ngwg = \int n_g(r)w_g(r)
   for(i17 = 0; i17 < numgal; i17++)
   {
      dngwg = gal[i17].w;
      ngwg = ngwg + dngwg;
   }

   //Calculate nrwr = \int n_r(r)w_r(r)
   for(i18 = 0; i18 < numran; i18++)
   {
      dnrwr = ran[i18].w;
      nrwr = nrwr + dnrwr;
   }

   alpha = ngwg/nrwr;
 
   //Calculate nw2 = \int \bar{n}(r)w_r^2(r)
   for(i19 = 0; i19 < numran; i19++)
   {
      di19 = ran[i19].w;
      dnw2 = di19*di19;
      nw2 = nw2 + dnw2;
   }   
   //Note additional factor of alpha and dropping of a factor of \bar{n} in the integrand in converting from integral to sum. This is explained in FKP before Eq.(2.4.1)
   nw2 = alpha*nw2;

   //Calculate n2w2 = \int \bar{n}^2w_r^2(r)
   for(i20 = 0; i20 < numran; i20++)
   {
      di20 = ran[i20].w;
      di20sq = di20*di20;
      dn2w2 = di20sq*ran[i20].n;
      n2w2 = n2w2 + dn2w2;
   }
   //Same as per above for additional factor of alpha missing factor of \bar{n} in the integrand
   n2w2 = alpha*n2w2;
   
   //FKP normalization
   normalization = n2w2;

   //FKP shot noise
   shotnoise = (1.0 + alpha)*nw2/(normalization);
   
   //Initialize array for galaxy and random catalogs, and array for weighted difference
   ngal = init_3D_array<double>(grid, grid, grid);
   nran = init_3D_array<double>(grid, grid, grid);
   fkp = init_3D_array<double>(grid,grid,grid);

   //Put galaxy and random catalogs on grid according to the CIC MAS
   cout << "CIC gridding beginning" << endl;
   cic(gal, numgal, ngal);
   cic(ran, numran, nran);
   cout << "CIC gridding completed" << endl;
   //Now ngal and nran contain galaxy and random catalogs weighted and placed on grid

   for(i21 = 0; i21 < grid; i21++)
   {
      for(i22 = 0; i22 < grid; i22++)
      {
         for(i23 = 0; i23 < grid; i23++)
         {
            fkp[i21][i22][i23] = ngal[i21][i22][i23] - alpha*nran[i21][i22][i23];
         }
      }
   }
   
   //Don't need ngal and nran anymore since we have their FKP weighted difference, free up memory
   kill_3D_array(ngal, grid, grid);
   kill_3D_array(nran, grid, grid);

   //Initialize arrays that will contain final binned values of k and P(k)
   kfinal = init_1D_array<double>(nbin);
   pfinal = init_1D_array<double>(nbin);

   //Compute power spectrum and bin in k (since ps function calls binner function)
   ps(fkp, normalization, shotnoise, kfinal, pfinal);

   //Print binned [k,P(k)] to file
   ofstream pkbinned;
   stringstream kpkout;
   kpkout << "../Outputs/" << myout << "." << myoutext;
   pkbinned.open( kpkout.str().c_str(), ios::out );
   pkbinned.precision(8);
   
   for(i24 = 0; i24 < nbin; i24++)
   {
      pkbinned << kfinal[i24] << ' ' << pfinal[i24] << endl;
   }
   
   pkbinned.close();

   cout << "Your power spectrum is awaiting you in Outputs/" << myout << "." << myoutext << endl;
   
   return 0;
}
