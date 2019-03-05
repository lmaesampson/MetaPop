
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <random>
#include <vector>
#include <array>
#include <functional>
#include <numeric>
#include <cassert>
#include <algorithm>
// #include <armadillo>

//compile with g++ -std=c++14 -O2 metapop.cpp lodepng.cpp  -o metapop

#include "lodepng.h"
#include "basiccolors.hpp"

using namespace std;

// Number of age classes:
constexpr int nac = 1;

////Number of columns, rows, total patches:
constexpr int ncl = 20;
constexpr int nrw = ncl;
constexpr int nptch = ncl * nrw;  //square grid.
//constexpr int mbs = nptch * nptch;
constexpr int lnumtarg = 1212;  //patch where infecteds are introduced
constexpr float pi = 3.14159265359;
constexpr int targnum = 10;     //numer of infecteds introduced

#include "mixingMatrices/BetaMc1.5_norm.hpp"
#define loadbmat 0


template <class T>
void printAsLine(T& c) {
    for( typename T::iterator i = c.begin(); i != c.end(); i++ ) {
	std::cout << *i << " " ;
    }
    // cout << endl;
}

std::ostream& nl(std::ostream& out){
  return out << "\n";
}


void writeFile(string fn, string content) {
    std::ofstream fb(fn);
    if (fb.is_open())
    {
	fb << content;
    }
}

vector<float> computeMeanAndStandardDev(vector<float> datap) {
    float max = 0.0;

    vector<float> data;
    vector<float> meanAndStdev;

    meanAndStdev.reserve(4);
    data.reserve(4);

    for(auto x : datap) {
	if (!(isnan(x))) {
	    data.push_back(x);
	}
    }
    float min = data[0];

    for(auto i : data) {
	if (i>max) {
	    max = i;
	}
	if (i<min) {
	    min = i;
	}
    }

    // vector<float> data( a,a + sizeof( a ) / sizeof( a[0] ) );
    float mean = accumulate( data.begin(), data.end(), 0.0f )/ data.size();
    vector<float> zero_mean( data );
    transform( zero_mean.begin(), zero_mean.end(), zero_mean.begin(),bind2nd( minus<float>(), mean ) );
    float deviation = inner_product( zero_mean.begin(),zero_mean.end(), zero_mean.begin(), 0.0f );
    deviation = sqrt( deviation / ( data.size() - 1 ) );
    meanAndStdev.push_back(min);
    meanAndStdev.push_back(max);
    meanAndStdev.push_back(mean);
    if(isnan(deviation)) {
	meanAndStdev.push_back(0);
    } else {
	meanAndStdev.push_back(deviation);
    }
    return meanAndStdev;
}

float computeStandardDev(vector<float> data) {
    // vector<float> data( a,a + sizeof( a ) / sizeof( a[0] ) );
    float mean = accumulate( data.begin(), data.end(), 0.0f )/ data.size();
    vector<float> zero_mean( data );
    transform( zero_mean.begin(), zero_mean.end(), zero_mean.begin(),bind2nd( minus<float>(), mean ) );
    float deviation = inner_product( zero_mean.begin(),zero_mean.end(), zero_mean.begin(), 0.0f );
    deviation = sqrt( deviation / ( data.size() - 1 ) );
    return deviation;
}


// Parameter set:
struct parset {
  float popvar;
  float patchpop;
  float ttsi; // total count of Infecteds
  float tscs; // total count of timesteps with cases
  string fn; // filename stem
  int images_on; // image switch
  int logs_on; // log switch
  int iters;
  int stoch;
  float betamean;
  int permiter;
  float beta;
  float conprob;
  float constren;
  float maxcasecount;
  float maxsusccount;
  vector<int> patchtargets;
  parset() {};
};

// Simple coordinate class, holding an index into 1D patch arrays:
struct crd {
    int x;
    int y;
    int n;
    crd() { };
    crd(int vx, int vy, int ndx) : x(vx), y(vy), n(ndx) { };
};

// Generate a sequence:
template<class T>
vector<T> genseq(T start, T stop, T inc) {
    vector<T> req;
    float cv = start;
    while(cv < stop) {
	req.push_back(cv);
	cv += inc;
    }
    return(req);
}

// Just sum over a nax by nptch 2D array:
float sumNpxNa( std::array<std::array<float, nac>, nptch>  m) {
    float r = 0.0;
    for(int i=0; i<nac; i++) {
	for(int j=0; j<nptch; j++) {
	    r += m[i][j];
	}
    }
    return r;
}

// return vector with each patch sum:
vector<float> sumEach( std::array<std::array<float, nac>, nptch> &m) {
    float r = 0.0;
    vector<float> rv;
    for(int i=0; i<nptch; i++) {
	r = 0.0;
	for(int j=0; j<nac; j++) {
	    r += m[j][i];
	}
	rv.push_back(r);
    }
    return rv;
}

// return vector with each patch sum:
void spp( std::array<std::array<float, nac>, nptch> &m) {
    for(int i=0; i<nac; i++) {
	for(int j=0; j<nptch; j++) {
	    // p i j m[i][j]
	    std::cout <<"["<<__FILE__<<"|"<<__FUNCTION__<<"|"<<__LINE__<<"]"<< "nac i: " << i <<"npt j: " << j <<" m[i][j]: " << m[i][j] << " " << std::endl;


	}
    }
}

// Sum a col by row array
float sumNrxNc(std::array<std::array<float, ncl>, nrw> &m) {
    float r = 0.0;
    for(int i=0; i<ncl; i++) {
	for(int j=0; j<nrw; j++) {
	    r += m[i][j];
	}
    }
    return r;
}


// Just sum by nac and reshape nptch 2D array:
std::array<std::array<float, ncl>, nrw>
sumReshape(std::array<std::array<float, nac>, nptch>  &m) {
    std::array<std::array<float, ncl>, nrw> rt;
    float r = 0.0;
    int p=0;
    int q=0;
    int nclc = 0;
    for(int j=0; j<nptch; j++) {
	r = 0.0;
	for(int i=0; i<nac; i++) {
	    r += m[i][j];
	}
	rt[p][q] = r;
	nclc += 1;
	if(nclc>ncl) {
	    p += 1;
	    q = 0;
	}
    }
    return rt;
}


float sumX(std::array<float, nptch> &m) {
    float r = 0.0;
    for(int i=0; i<ncl; i++) {
	r += m[i];
    }
    return r;
}



void mat2csv(const std::array<std::array<float, ncl>, nrw> &m,
	     std::string n
	    ) {

    stringstream a;

    for(int u=0; u<nrw; u++) {
	for(int j=0; j<ncl; j++) {
	    a << floor(m[u][j]) << " ";
	}
	a << endl;
    }

    std::stringstream a_fn_ss;
    a_fn_ss << n<< ".csv";
    std::string a_fn_str = a_fn_ss.str();
    std::ofstream a_ofst(a_fn_str.c_str());
    if (a_ofst.is_open())
    {
	a_ofst << a.str();
    }
    else
    {
	std::cerr << "Could not open " << n << ".csv" << std::endl;
    }
}

//The image argument has width * height RGBA pixels or width * height * 4 bytes
void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height)
{
    //Encode the image
    unsigned error = lodepng::encode(filename, image, width, height);
    //if there's an error, display it
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

//void arpng(const arma::mat& pnar, const vector<int>& cm, std::string fs){
void arpng(const std::array<std::array<float, nac>, nptch> &pnar, const vector<int>& cm,float maxcount, std::string fs){
  //  mat ar = normalise(pnar);
  unsigned xwidth  = nrw ;
  unsigned xheight = ncl;
  int i = 0;
  // for ffmpeg
  // while(xwidth % 2 != 0){ xwidth--;};
  // while(xheight % 2 != 0){ xheight--;};
  unsigned int imscalar = 4;
  static vector<unsigned char> ximage;
  ximage.resize(xwidth*xheight*imscalar);
  for(unsigned y = 0; y < xheight; y++){
      for(unsigned x = 0; x < xwidth;  x++) {
	// double w = 255.0 * ar(y,x);
	double w = 100.0 * pnar[0][i]/maxcount;
// p w pnar[0][i] maxcount
	unsigned uw = static_cast<unsigned>(w);
	int ci = uw*3;
	ximage[imscalar * xwidth * y + imscalar * x + 0] = cm[ci];
	ximage[imscalar * xwidth * y + imscalar * x + 1] = cm[ci+1];
	ximage[imscalar * xwidth * y + imscalar * x + 2] = cm[ci+2];
	ximage[imscalar * xwidth * y + imscalar * x + 3] = 255;
	i++;
      }
    }
    std::stringstream bmgss;
    bmgss << fs << ".png";
    std::string bmg_str = bmgss.str();
    unsigned error = lodepng::encode(bmg_str.c_str(), ximage, xwidth, xheight);
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}



//void arpng(const arma::mat& pnar, const vector<int>& cm, std::string fs){
void arpng2(const std::array<std::array<float, nac>, nptch> &pnar,
	    const std::array<std::array<float, nac>, nptch> &pnar2,
	    const vector<int>& cm,
	    float maxcount,
	    float maxcount2,
	    std::string fs){
  //  mat ar = normalise(pnar);
  unsigned xwidth  = nrw * 2;
  unsigned xheight = ncl;
  int i = 0;
  int i2 = 0;
  // for ffmpeg
   while(xwidth % 2 != 0){ xwidth--;};
   while(xheight % 2 != 0){ xheight--;};
  unsigned int imscalar = 4;
  static vector<unsigned char> ximage;
  ximage.resize(xwidth*xheight*imscalar);
  for(unsigned y = 0; y < xheight; y++){
      for(unsigned x = 0; x < xwidth;  x++) {

	double w = 0;
	if(x<xwidth/2){
	   w = 100.0 * pnar[0][i]/maxcount;
	   i++;
	   //	   cout << "o";
	}else{
	   w = 100.0 * pnar2[0][i2]/maxcount2;
	   i2++;
	   //	   cout << "x";
	}

	unsigned uw = static_cast<unsigned>(w);
	int ci = uw*3;
	int hp = imscalar * xwidth * y + imscalar * x + 0;

    //cout << imscalar << endl;
    if(imscalar>=1){
	ximage[imscalar * xwidth * y + imscalar * x + 0] = cm[ci];
	ximage[imscalar * xwidth * y + imscalar * x + 1] = cm[ci+1];
	ximage[imscalar * xwidth * y + imscalar * x + 2] = cm[ci+2];
	ximage[imscalar * xwidth * y + imscalar * x + 3] = 255;
	//	i++;
      }
          if(imscalar<1){
              ximage[imscalar * xwidth * y + imscalar * x + 0] = cm[ci];
              ximage[imscalar * xwidth * y + imscalar * x + 1] = cm[ci+1];
              ximage[imscalar * xwidth * y + imscalar * x + 2] = cm[ci+2];
              ximage[imscalar * xwidth * y + imscalar * x + 3] = 255;
              //	i++;
          }
      }
    }

    std::stringstream bmgss;
    bmgss << fs << ".png";
    std::string bmg_str = bmgss.str();
    unsigned error = lodepng::encode(bmg_str.c_str(), ximage, xwidth, xheight);
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}


// Modified for betaMatrix
////////////////////////////////////////////////////////////////////////////////////////////////////
//This function takes a vector of susceptible individuals in each patch
//(S), a vector of infected individuals in each patch (Infe)
//	a single value for transmission (Beta), the alpha parameter in a discrete-time TSIR model (alpha)
//	a number of births (births), and a vector of vaccination rates in ecah patch (vax)
//This function steps through one time step of a TSIR model
//This function returns a vector of susceptible individuals in each
//patch (S), and a vector of Infected individuals in each patch (I)
vector< std::array<std::array<float, nac>, nptch> >
epistepTSIRVmetaAge(
    const std::array<std::array<float, nac>, nptch> &S,
    const std::array<std::array<float, nac>, nptch> &I,
    const std::array<std::array<float, nac>, nptch> &R,
    const std::array<std::array<float, nac>, nptch> &V,
    const float Beta,
    //std::array<float, nptch> Betas,
    const std::array<std::array<float, nptch>, nptch> &betaMatrix,
    const float alpha,
    const float births,
    //const std::array<const std::array<float, nac>, nptch> vax,
    const std::array<float, nptch> vax,
    parset &pars,
    const vector<int> &cpt
)
{
    /* It <- array( 0, dim=c(nptch,age.classes)); */
    /* St <- array( 0, dim=c(nptch,age.classes)); */
    std::array<std::array<float, nac>, nptch> It, St, Rt, Vt;

    int ii = 0;

    std::binomial_distribution<int> bindist(9,0.5);

    float tsi = 0.0;

    float loc_S = 0.0f;
    float loc_I = 0.0f;
    float loc_R = 0.0f;
    float loc_V = 0.0f;


    //iterate through patches
    for(int ix=0; ix<nptch; ix++) {

	ii = cpt[ix];

	float possibleIs = 0.0;
	float possibleSs = 0.0;
	float possibleRs = 0.0;
	float possibleVs = 0.0;

	for(int j=0; j<nac; j++) {

	    float fullBetaI = 0.0;
	    float bi = 0.0;
        float mRate, vRate;
	    possibleIs = 0.0;
	    possibleSs = 0.0;

	    possibleRs = 0.0;
	    possibleVs = 0.0;

	    loc_S = S[ii][j];  //S,I,R,V in this patch
	    loc_I = I[ii][j];
	    loc_R = R[ii][j];
	    loc_V = V[ii][j];

	    float tpop = loc_S + loc_I + loc_R + loc_V;

        float BetaN = (Beta*pars.betamean)/tpop; //  "Beta" here is merely the seasonal component
        //float BetaN = pars.betamean/tpop;
        
        

	    // S′ = −βSI,
	    possibleSs = loc_S + (-BetaN * loc_S*loc_I);
	    // R′ = αI
	    possibleRs = loc_R + (alpha * loc_I);

        /*float v1 = rand();
        if(v1>0.99){
            possibleIs += 1;
        }*/
        //Stochastic introduction: draw an integer from a Poisson distribution with rate equal to the sum of virgin introduction and migratory introduction
        mRate = 0.;
        for(int nx=0; nx<nptch; nx++) {
            if(betaMatrix[ii][nx] > 0.0) { //adds up all beta*I for all patches that have non-zero connectivity with this patch
                if(nx!=j)
                {mRate += betaMatrix[ii][nx] * I[nx][j];}  // beta*I}
                if(nx==j)
                {mRate += 0.;}
            }
        }
        
        mRate /= 5.;
        vRate = 0.0001; //If the virgin introduction rate is always set to zero (which it has been for many simulations), do you lose the 'white noise' effect that changes ecc?
        
        if((vRate+mRate)<0.0001)
        {
            vRate = 0.0001;
        }
        
        // std::default_random_engine generator(time(0)*(ii+j));
        std::random_device r;
        std::seed_seq seed_seq{r()};
        std::mt19937 generator{seed_seq};
        
        std::poisson_distribution<int> pDist(vRate+mRate);
        int intros = pDist(generator);
        
        // I′ =βSI−αI,
        possibleIs = loc_I + (BetaN * loc_S*(loc_I+intros)) - (alpha * loc_I);
        
        if(possibleIs < 0 ) {
            //assert(possibleIs > 0 );
            possibleIs = 0;
        }
        if(possibleIs>loc_S) {
            //assert(possibleIs<loc_S);
            possibleIs=loc_S;
        }
        
        tsi += possibleIs;
       
        
        //std::cout << "generator:" << generator << std::endl;
        //std::cout <<"rate: " << mRate+vRate << " intros " << intros << std::endl;
        
        possibleIs += intros;

	    It[ii][j] = possibleIs;
	    Rt[ii][j] = possibleRs;

	    // Vaccination:
	    float tsbirths = (tpop / 1000.0) * (births/12.0); // 12 because br is measured at half year. births per thousand per six months (half a year is 12 time steps)
	    float ftsbirths = tsbirths;
	    float vaccedbirths = (vax[ii])*(ftsbirths);       //vaccination is implemented in the birth rate

	    Vt[ii][j] = V[ii][j] + vaccedbirths;
	    St[ii][j] = possibleSs + (ftsbirths - vaccedbirths);
	    if(St[ii][j]<0) {
		St[ii][j]=0;
	    }

	    // DEATH
	    float deathscalar = tsbirths/(tpop+tsbirths);  //keeps population constant

	    float Vtd = (Vt[ii][j]*deathscalar);
	    float Itd = (It[ii][j]*deathscalar);
	    float Rtd = (Rt[ii][j]*deathscalar);
	    float Std = (St[ii][j]*deathscalar);

	    Vt[ii][j] -= Vtd;
	    It[ii][j] -= Itd;
	    Rt[ii][j] -= Rtd;
	    St[ii][j] -= Std;

	    //track max number of cases:
	    if(pars.maxcasecount<It[ii][j]){    pars.maxcasecount=It[ii][j];	    }
	    if(pars.maxsusccount<St[ii][j]){    pars.maxsusccount=It[ii][j];	    }

	}
    }

    pars.ttsi += tsi;
    if(tsi > 0.0) pars.tscs += 1;

    vector<std::array<std::array<float, nac>, nptch> > retStor;
    retStor.push_back(St);
    retStor.push_back(It);
    retStor.push_back(Rt);
    retStor.push_back(Vt);
    return(retStor);

}



////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  epiTSIRmeta
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//This function takes a vector of the initial number of susceptible
//individuals in each patch (S),
//	a vector of the initial number of infected individuals in each
//	patch (Infe), a matrix of the vaccination rate at each
//	timestep (row) and patch (column) a vector of transmission at
//	each timestep (Beta), the alpha parameter in a discrete-time
//	TSIR model (alpha) a number of births (births), and a number
//	of timesteps (T)
//This function runs a TSIR model for a range of patches This function
//returns a matrix of susceptible individuals in each patch (column)
//at each timestep (row) (S),
//	and a matrix of Infected individuals in each patch (column) at
//	each timestep (row) (I)
const vector< vector< std::array<std::array<float, nac>, nptch> > >   //time stepping function.
epiTSIRVmeta(
    const std::array<std::array<float, nac>, nptch> &S,
    const std::array<std::array<float, nac>, nptch> &I,
    const std::array<std::array<float, nac>, nptch> &R,
    const std::array<std::array<float, nac>, nptch> &V,
    const std::array<float, nptch> vax,
    const vector<float> bet,
    const float alpha,
    const float births,
    const int T,
    parset &pars
)
{

    //I.stor<-array(0, dim = c((T+1), nptch, age.classes))
    //S.stor<-array(0, dim = c((T+1), nptch, age.classes))

    vector< std::array<std::array<float, nac>, nptch> > Sstor, Istor, Rstor, Vstor,out;

    Sstor.reserve(T);
    Istor.reserve(T);
    Rstor.reserve(T);
    Vstor.reserve(T);

    Istor.push_back(I);
    Sstor.push_back(S);
    Rstor.push_back(R);
    Vstor.push_back(V);

    std::array<std::array<crd, nrw>, ncl> lup;
    std::array<crd, nptch> lnpt;
    int cx = 0;
    //stringstream a;
    for(int  i=0; i<nrw; i++) {
	for(int j=0; j<ncl; j++) {     //lup is a lookup table - gives you mapping between i,j and the patch number, and gives the neighbors. Creating lookup table here.
	    lup[i][j] = crd(i,j,cx);
	    // a << "(" << i << "," << j <<","<<cx<<") " ;
	    lnpt[cx] = crd(i,j,cx); //coordinate structure.
	    cx++;
	}
	//a << endl;
    }
    // // // This will print the mapping of the flat array to the grid:
    // std::ofstream lup_ofst("lup.txt");
    // if (lup_ofst.is_open()){
    //    lup_ofst << a.str();
    //  }

    // Create betaMatrix, load with zeros
    static std::array<std::array<float, nptch>, nptch> betaMatrix;
    if(loadbmat==0){
        static std::array<std::array<float, nptch>, nptch> sbetaMatrix;
    }
    auto rand_unif = bind(uniform_real_distribution<> {0.0,1.0},default_random_engine{});   //random mixing matrix, or can load betaMatrix from a file.
    if(loadbmat==1) {
      betaMatrix = sbetaMatrix;
    }
    else
    {
	float obmdiv = 0.1;
	float bmdiv = 1.0 - (4.0 * obmdiv);  //setting up standard deviations.
	for(int  i=0; i<nptch; i++) {
	    for(int j=0; j<nptch; j++) {
		betaMatrix[i][j] = 0.0;
		if(rand_unif()<pars.conprob)
        {
            betaMatrix[i][j] = pars.constren;
        }//questionable, apparently. based on connection probability and connection strength. Lets things be completely unconnected.
		if(i==j) {
		    betaMatrix[i][j] = bmdiv;
		}
            //cout << betaMatrix[i][j] << pars.conprob << pars.constren <<endl;

	    }
	}

	int ix, ixn,ixe,ixw,ixs;
	for(int  i=0; i<ncl; i++) {
	    for(int j=0; j<nrw; j++) {

		ix = lup[i][j].n;

		if(i==j) {
		    betaMatrix[ix][ix] = bmdiv;
		}

		if(i>0) {
		    ixe = lup[i-1][j].n;
		    betaMatrix[ixe][ix] = obmdiv;
		}
		if(i<(ncl-1)) {
		    ixw = lup[i+1][j].n;
		    betaMatrix[ixw][ix] = obmdiv;
		}
		if(i<=(ncl-2)) {
		    ixw = lup[i-nrw][j].n;
		    if(ixw>=0&&ixw<=nptch){
// p ixw ix npatch
		    betaMatrix[ixw][ix] = obmdiv;
		    }
		}
		if(j>0) {
		    ixn = lup[i][j-1].n;
		    betaMatrix[ixn][ix] = obmdiv;
		}
		if(i<=(nrw-2)) {
		    ixs = lup[i][j+1].n;
		    betaMatrix[ixs][ix] = obmdiv;
		}
		// p i j ix
	    }
	}
    }



    vector<int> cpt;  //cpt is a vector of patch numbers. we start at a different patch number at every time step.
    cpt.reserve(nptch);
    for(int i=0; i<nptch; i++) {
	cpt.push_back(i);
    }

    int imc = 0;
    //iterate through timesteps
    for(int i=0; i<T; i++) {
	//cout << "T:" << i << endl;

	random_shuffle(cpt.begin(), cpt.end()); //shuffle the patch numbers to start at a different patch each timestep

	out = epistepTSIRVmetaAge( //this function does one timestep
		  Sstor[i],
		  Istor[i],
		  Rstor[i],
		  Vstor[i],
		  bet[i],
		  betaMatrix,
		  alpha,
		  births,
		  vax,
		  pars,
		  cpt
	      );
	// #store output in matrices
	// S.stor[(i+1),,]<-out$S I.stor[(i+1),,]<-out$I
	// Sstor[(i+1)] = out[0];  Istor[(i+1)] = out[1];

	Sstor.push_back(out[0]);
	Istor.push_back(out[1]);
	Rstor.push_back(out[2]);
	Vstor.push_back(out[3]);

	// if images:
	if(pars.images_on) {
	    // if(i>=0) {

	    //	std::array<std::array<float, nac>, nptch> Iu = out[1];
	    //	vector<float> gdi; // img
	    //	float bgv = 0.0;
	    //	float mbgv = 0.0;
	    //	gdi.reserve(nptch);
	    //	for(int ip=0; ip<nptch; ip++) {
	    //	    float asum = 0.0;
	    //	    for(int ag=0; ag<nac; ag++) {
	    //		asum += Iu[ip][ag];
	    //	    }
	    //	    if(asum>bgv) {
	    //		bgv = asum;
	    //	    }
	    //	    if(mbgv>asum || i==0) {
	    //		mbgv = asum;
	    //	    }
	    //	    gdi.push_back(asum);
	    //	}

	    //	std::array<std::array<float, nrw>, ncl> dmit;

	    //	std::stringstream a;
	    //	int gdr = 0;
	    //	for(int u=0; u<nrw; u++) {
	    //	    for(int m=0; m<ncl; m++) {
	    //		a << floor(gdi[gdr]) << " ";
	    //		gdr++;
	    //		dmit[u][m] = gdi[gdr];
	    //	    }
	    //	    a << endl;
	    //	}

	    //	std::array<std::array<float, nac>, nptch> Su = out[0];
	    //	vector<float> sgdi; // img
	    //	float sbgv = 0.0;
	    //	float msbgv = 0.0;
	    //	sgdi.reserve(nptch);
	    //	for(int ip=0; ip<nptch; ip++) {
	    //	    float asum = 0.0;
	    //	    for(int ag=0; ag<nac; ag++) {
	    //		asum += Su[ip][ag];
	    //	    }
	    //	    if(asum>sbgv) {
	    //		sbgv = asum;
	    //	    }
	    //	    if(msbgv>asum || i==0) {
	    //		msbgv = asum;
	    //	    }

	    //	    sgdi.push_back(asum);
	    //	}

	    //	std::array<std::array<float, nrw>, ncl> sdmit;

	    //	gdr = 0;
	    //	for(int u=0; u<nrw; u++) {
	    //	    for(int m=0; m<ncl; m++) {
	    //		gdr++;
	    //		sdmit[u][m] = sgdi[gdr];
	    //	    }
	    //	}

	    //	vector<float> xgdi; // img
	    //	float xsbgv = 0.0;
	    //	for(int ip=0; ip<nptch; ip++) {
	    //	    float asum = 0.0;
	    //	    float bsum = 0.0;
	    //	    for(int ag=0; ag<nac; ag++) {
	    //		asum += Su[ip][ag];
	    //		bsum += Iu[ip][ag];
	    //	    }
	    //	    if(asum>xsbgv) {
	    //		xsbgv = asum;
	    //	    }
	    //	}

	    //	std::array<std::array<float, nrw>, ncl> xdmit;
	    //	gdr = 0;
	    //	for(int u=0; u<nrw; u++) {
	    //	    for(int m=0; m<ncl; m++) {
	    //		gdr++;
	    //		xdmit[u][m] = sgdi[gdr];
	    //	    }
	    //	}


	    //	float csc = 255.0 / bgv;
	    //	float ssc = 255.0 / sbgv - msbgv;

	    //	unsigned width  = nrw;
	    //	unsigned height = ncl;
	    //	unsigned int imscalar = 4;

	    //	std::vector<unsigned char> image;
	    //	image.resize(width * height * imscalar);

	    //	for(unsigned y = 0; y < height; y++) {
	    //	    for(unsigned x = 0; x < width;  x++) {
	    //		image[imscalar * width * y + imscalar * x + 0] = 255;
	    //		image[imscalar * width * y + imscalar * x + 1] = 255 - (csc * static_cast<unsigned>(floor(  dmit[y][x])));
	    //		image[imscalar * width * y + imscalar * x + 2] = 255 - (ssc * static_cast<unsigned>(floor( sdmit[y][x])));
	    //		image[imscalar * width * y + imscalar * x + 3] = 255;
	    //	    }
	    //	}

	    //	string ist = to_string(imc);
	    //	int gen00len = 5 - ist.length();
	    //	string ztring;
	    //	string zr("0");
	    //	for(int z=0; z<gen00len; z++) {
	    //	    ztring.append(zr);
	    //	}
	    //	std::stringstream p_fn_ss;
	    //	p_fn_ss << "out_" << pars.fn  << "_" << ztring  << imc << ".png";
	    //	std::string p_fn_str = p_fn_ss.str();
	    //	encodeOneStep(p_fn_str.c_str(), image, width, height);

	    //	imc++;
	    // }
	}

    }

    vector< vector<std::array<std::array<float, nac>, nptch> > > retStor;
    retStor.reserve(4);
    retStor.push_back(Sstor);
    retStor.push_back(Istor);
    retStor.push_back(Rstor);
    retStor.push_back(Vstor);

    return(retStor);
}

// This function takes a number of patches (n), a base vaccination rate
// (prob), a number of timesteps to run the TSIR model for (timesteps) --
// plots/output will only be for the second half a patch population
// (patchpop, 2e5 by default), and a mean beta (beta.mean - 20 by
// default (beta is basically R0)) This function generates epi
// parameters and plots output This function returns a matrix of actual
// vaccination rates (vax), and the distribution of phase angles
// (phase)
void epimeta(int iter,
	     const float vmaxprob,
	     const float vminprob,
	     const int timesteps,
	     string filename,
	     const float patchpop,
	     const float births,
	     parset &pars
	    )
{
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();

    //float betameanPP	= betamean/patchpop;
    float alpha		= 0.986;
    float seas		= 0.0;
    int T		= timesteps;

    vector<float> vax, sax;
    vax.reserve(nptch);
    sax.reserve(nptch);

    // float minvax = vminvax;
    // float maxvax = prob;

    vax = genseq<float>(vmaxprob,vminprob,((vmaxprob-vminprob)*((1.0/static_cast<float>(nptch))))); //vaccination set in each patch. Right now it's a smooth gradient. If you want a Gaussian, that's how pop is set.
    sax = genseq<float>(0.0,1.0,((1.0-0.0)*((1.0/static_cast<float>(nptch))))); //scale across space.

    std::array<float, nptch> matvax;

    myclock::duration d = myclock::now() - beginning;
    unsigned seed2 =  d.count();
    std::default_random_engine generator(seed2); //random number generator with random seed.

    if(pars.stoch==0) {
	std::default_random_engine generatorns; // no seed
	generator = generatorns;  //random number generator - with no seed.
    }

    std::normal_distribution<float> vdist(vmaxprob,vminprob);  //This line sets the vacc distribution as a normal dist centered at vmaxprob, with variance vminprob
    for(int i=0; i<nptch; i++)
    {
	float lvac = vdist(generator);
        if(lvac>=1.){lvac=0.99;}
        if(lvac>vmaxprob){lvac=vmaxprob;}
        matvax[i] = lvac; // * sax[i]; //for now, no spatial gradient.
        //cout  << lvac << endl;

    }

    vector<float> agepopscaler;
    agepopscaler = genseq<float>(5.0,60.0,2.0);

    //S <- array(0,dim=c(n,nac))# nac = age classes
    std::array<std::array<float, nac>, nptch> S,I,R,V;   //Create array (nagexnpatch) for S,I,R,V

    std::normal_distribution<float> ndistribution(patchpop,(patchpop * pars.popvar)); //popvar is the percentage standard deviation

    if(pars.popvar>0.0) {
	for(int i=0; i<nptch; i++)
	{
	    // S[i,] = round(patchpop*(1-vax[i])/seq(5,29,1) )
	    float gpatchpop  = ndistribution(generator);
	    if(gpatchpop <= 0) {
		gpatchpop = 1000;
	    }
	    for(int j=0; j<nac; j++) {
		// S[i][j] = gpatchpop*(1-matvax[i]); // application of vaccination
		S[i][j] = gpatchpop;
		I[i][j] = 0.0;
		// R[i][j] = 0.0;
		//V[i][j] = 0.0;
	    }
	}
    }
    else
    {
	for(int i=0; i<nptch; i++)
	{
	    for(int j=0; j<nac; j++) {
		S[i][j] = pars.patchpop;
		I[i][j] = 0.0;
		// R[i][j] = 0.0;
		// V[i][j] = 0.0;
	    }
	}
    }

    R = I;
    V = I;

    // Seed initial I population:
    //someIs <- round(rnorm(n=(nptch*nac),mean=10,sd=2))
    //I <- array(someIs,dim=c(nptch,nac))

    //int lnum = lnumtarg;     //this is the patch into which we introduce cases

    vector<int> patchtargets = pars.patchtargets;
    for(int i=0; i<nptch; i++)
    {
	for(int j=0; j<nac; j++)
	{
	    I[i][j] = 0.0;
	    if(std::find(patchtargets.begin(), patchtargets.end(), i) != patchtargets.end()) {
		float number = targnum;  //add targnum to the infecteds
		if(number<0) {
		    number = 0;
		}
		I[i][j] = number;
		S[i][j] -= number;
           // I[i][j] = 0.0;
	    }
	}
    }

    //##transform beta to a vector with seasonal forcing
    //beta.seas<-beta.mean*(1+seas*cos(2*pi*(1:24)/24))
    vector<float> singlebetaseas; // 1..24
    for(int i=0; i<24; i++)
    {
	//betaseas[bi]=betameanPP*(1.0+seas*math::cos(2.0*math::pi*(bi)/24.0))
	float bi = static_cast<float>(i);
	// Remove betameanPP:
	//float nv = betameanPP * (1.0 + seas * cos(2.0*pi*((bi)/24.0)));
	// Just save the seasonal component:
    float nv =  (1.0 + seas * cos(2.0*pi*((bi)/24.0)));///(1.+seas);  //vector that's really just a cosine. 24 timesteps per year.
	singlebetaseas.push_back(nv);
    }

    vector<float> betaseas;
    for(int i=0; i<timesteps; i++)
    {
	betaseas.insert(std::end(betaseas), std::begin(singlebetaseas), std::end(singlebetaseas));
    }

    //out = epiTSIRVmeta(S,I,vax,BETA,alpha,births,T,nptch,nr,iter)
    vector< vector< std::array<std::array<float, nac>, nptch> > > out = epiTSIRVmeta(S,I,R,V,
	    matvax,betaseas,alpha,
	    births,T, pars); //This is where the model starts

    // --> The returned 'out' stacks:
    // out[0][i] S
    // out[1][i] I
    // out[2][i] R
    // out[3][i] V

    if(pars.logs_on) {
	unsigned ol = out[0].size();    //integer is the length of the vector in the first position of out (which happens to be S). Should be number of timesteps.

	vector<float> sumS, sumI,sumR, sumV;

	sumS.reserve(ol); //allocates the memory
	sumI.reserve(ol);
	sumR.reserve(ol);
	sumV.reserve(ol);


	string ist;
	int gen00len;
	string ztring;
	string zr("0");
        


	/*for(unsigned int i=0; i<ol; i++) {
	    float sums = sumNpxNa(out[0][i]); //function defined to sum through patches and ages
	    float sumi = sumNpxNa(out[1][i]);
	    float sumr = sumNpxNa(out[2][i]);
	    float sumv = sumNpxNa(out[3][i]);
	    sumS.push_back(sums); //adds this into the pre-defined vectors
	    sumI.push_back(sumi);
	    sumR.push_back(sumr);
	    sumV.push_back(sumv);
	    // write cases image:
        

	     ist = to_string(i);
	     gen00len = 5 - ist.length();
	     ztring.clear();
	    for(int z=0; z<gen00len; z++) {
	      ztring.append(zr);
	    }

	    stringstream png_bs;
	    png_bs << "out_" << filename  << "_" << ztring << i;
	    //	    arpng2(out[1][i],out[0][i],moreland,pars.maxcasecount,pars.maxsusccount,png_bs.str());
	    arpng2(out[1][i],out[0][i],cubehelix,pars.maxcasecount,pars.maxsusccount,png_bs.str());
	}*/

	std::stringstream sumSs, bss, sbss, ibss, vbss, rbss, pbss; //defines a string you can quickly append to

	vector<float> irv,srv,vrv,rrv;

	irv.reserve(T);
	srv.reserve(T);
	vrv.reserve(T);
	rrv.reserve(T);

	sumSs << "pop, sumS, sumI, sumR, sumV" << nl; //header is the first line of the stringstream
        

	int npm1 = nptch - 1 ;

	for(unsigned int i=0; i<ol; i++) {
	    float psz = (sumS[i]) + (sumI[i]) + (sumR[i]) + (sumV[i]) ;
	    sumSs << floor(psz) << ", " << floor(sumS[i]) << ", "
		  << floor(sumI[i]) << ", "  << floor(sumR[i]) << ", " << floor(sumV[i]) << nl; //same syntax as writing to cout

	    irv = sumEach(out[1][i]); //another summing function
	    srv = sumEach(out[0][i]);
	    rrv = sumEach(out[2][i]);
	    vrv = sumEach(out[3][i]);

	    ibss << " " << floor(irv[0]) ;
	    sbss << " " << floor(srv[0]) ;
	    rbss << " " << floor(rrv[0]) ;
	    vbss << " " << floor(vrv[0]) ;
	    pbss << " " << floor(irv[0]+srv[0]+rrv[0]+vrv[0]) ;
	    for(int p=1; p<npm1; p++) {
		ibss << " " << floor(irv[p]);
		sbss << " " << floor(srv[p]);
		rbss << " " << floor(rrv[p]);
		vbss << " " << floor(vrv[p]);
		pbss << " " << floor(irv[p]+srv[p]+rrv[p]+vrv[p]);
	    }
	    ibss << " " << floor(irv[npm1])  << nl;
	    sbss << " " << floor(srv[npm1])  << nl;
	    rbss << " " << floor(rrv[npm1])  << nl;
	    vbss << " " << floor(vrv[npm1])  << nl;
	    pbss << " " << floor(irv[npm1]+srv[npm1]+rrv[npm1]+vrv[npm1])  << nl;

	}


	//stringstream bx; //creating file names
	//bx << "sums_" << filename  << "_" << iter << ".csv"; //filename is passed in
	//std::ofstream sof(bx.str());
	//if (sof.is_open()) {
	  //  sof << sumSs.str();
	//}
	stringstream fn_bs;
	fn_bs << "cases_" << filename  << "_" << iter << ".csv";
	writeFile(fn_bs.str(),ibss.str()); //Brian has his own writeFile function
	stringstream fn_sbs;
	fn_sbs << "susc_" << filename  << "_" << iter << ".csv";
	writeFile(fn_sbs.str(),sbss.str());
	stringstream fn_rbs;
	fn_rbs << "recov_" << filename  << "_" << iter << ".csv";
	writeFile(fn_rbs.str(),rbss.str());
	stringstream fn_vbs;
	fn_vbs << "vacc_" << filename  << "_" << iter << ".csv";
	writeFile(fn_vbs.str(),vbss.str());
	stringstream fn_pbs;
	//fn_pbs << "pop_" << filename  << "_" << iter << ".csv";
	//writeFile(fn_pbs.str(),pbss.str());

    }
}

int main(int argc, char** argv)
{
   /* if(loadbmat==1)
    {
    #include "betamtest.hpp"
    }*/

    int timesteps = 420;

    float vaccprob	= 0.99; // Max vacc rate, presently
    float vaccprob2	= 0.9;// Min vacc rate, presently
    float betamean	= 8.0;
    float patchpop	= 10000;
    float births	= 40.0; // Births is now count per thousand

    // Birth rate compares the average annual number of births during a
    // year per 1,000 persons in the population at midyear; also known as
    // crude birth rate. [6-44]

    parset pars;
    pars.popvar = 0;
    pars.images_on = 1;

    string arg;

    for (int i = 1; i < argc; i+=2) {
	if (i + 1 != argc) {
	    if (strncmp(argv[i],"beta",4) == 0) {
		arg = string(argv[i+1]);
		betamean = stof(arg);
		pars.betamean = betamean;
	    } else if (strncmp(argv[i],"maxvaccprob",8) == 0) {
		arg = string(argv[i+1]);
		vaccprob = stof(arg);
	    } else if (strncmp(argv[i],"minvaccprob",9) == 0) {
		arg = string(argv[i+1]);
		vaccprob2 = stof(arg);
	    } else if (strncmp(argv[i],"timesteps",9) == 0) {
		arg = string(argv[i+1]);
		timesteps = stoi(arg);
	    } else if (strncmp(argv[i],"patchpop",8) == 0) {
		arg = string(argv[i+1]);
		patchpop = stof(arg);
		pars.patchpop = patchpop;
	    }
	    else if (strncmp(argv[i],"birthrate",9) == 0) {
		arg = string(argv[i+1]);
		births = stof(arg);
	    }
	    else if (strncmp(argv[i],"popstddev",9) == 0) {
		arg = string(argv[i+1]);
		pars.popvar = stof(arg);
	    }
	    else if (strncmp(argv[i],"conprob",7) == 0) {
		arg = string(argv[i+1]);
		pars.conprob = stof(arg);
	    }
	    else if (strncmp(argv[i],"constren",8) == 0) {
		arg = string(argv[i+1]);
		pars.constren = stof(arg);
	    }
	    else if (strncmp(argv[i],"images",6) == 0) {
		arg = string(argv[i+1]);
		pars.images_on = stoi(arg);
	    }
	    else if (strncmp(argv[i],"logs",4) == 0) {
		arg = string(argv[i+1]);
		pars.logs_on = stoi(arg);
	    }
	    else if (strncmp(argv[i],"fn",2) == 0) {
		arg = string(argv[i+1]);
		pars.fn = arg;
	    }
	    else if (strncmp(argv[i],"stochastic",10) == 0) {
		arg = string(argv[i+1]);
		pars.stoch = stoi(arg);
	    }
	    else if (strncmp(argv[i],"iter",4) == 0) {
		arg = string(argv[i+1]);
		pars.iters = stoi(arg);
	    }
	    else if (strncmp(argv[i],"patchtargets",11) == 0) {
	      break;
	    } else {
		cout << i << endl;
		cout << argv[i] << endl;
		cout << "Invalid arguments\n";
		exit(0);
	    }
	}
    }


    vector<int> patchtargets = {0,5};
    // get patch targets to seed infection:
    // for (int i = 1; i < argc; i++) {
    //   if (i + 1 != argc) {
    //	if (strncmp(argv[i],"patchtargets",11) == 0) {
    //	  while(i + 1 != argc) {
    //	    arg = string(argv[i+1]);
    //	    int ttar = stoi(arg);
    //	    patchtargets.push_back(ttar);
    //	    i++;
    //	  }
    //	}
    //   }
    // }

	pars.patchtargets = patchtargets;
    int numberofiterates = pars.iters;

    pars.permiter = 0;

    vector<float> res;
    vector<float> fres;
    res.reserve(numberofiterates);
    fres.reserve(numberofiterates);

    for(int i=0; i<numberofiterates; i++) {
	pars.ttsi = 0.0;
	pars.tscs = 0.0;
	epimeta(i, vaccprob, vaccprob2, timesteps, pars.fn, patchpop, births, pars); //This is Amalie's function
	res.push_back(pars.ttsi);
	fres.push_back(pars.tscs);
    }

    cout << " "<< vaccprob << " " << vaccprob2;
    cout <<" " << timesteps << " " << pars.fn << " " << patchpop << " " << betamean << " ";
    cout  << births << " " << pars.popvar << " " << pars.conprob << " " << pars.constren << " ";

//   To print to screen: cout  << variable << endl;

    if(pars.iters>1) {
	vector<float> fsdv;
	fsdv.reserve(4);
	fsdv = computeMeanAndStandardDev(fres);
	cout << fsdv[0] << " " << fsdv[1] << " " << fsdv[2] << " " << fsdv[3] << " " ;

	vector<float> sdv;
	sdv.reserve(4);
	sdv = computeMeanAndStandardDev(res);
	cout << sdv[0] << " " << sdv[1] << " " << sdv[2] << " " << sdv[3] << " " ;

	//reward
	// cout << (( (( sdv[1] - sdv[0] ) / sdv[1] ) + (( fsdv[1] - fsdv[0] ) / fsdv[1] ) ) / 2.0) - 1.0;
	// cout << (( sdv[1] - sdv[0] ) / sdv[1] )  - 1.0;
	cout << sdv[3];

    }
    cout << " " << endl;
    res.clear();
    fres.clear();

    // 0 replicate number
    //1 vaccprob
    //2 vaccprob2
    //3 timesteps
    //4 fn
    //5 patchpop
    //6 betamean
    //7 births
    //8 div
    //9 popvar
    //10 cases_count 1
    // ...
    //19 cases_count 10
    //20 min
    //21 max
    //22 mean
    //23 sd


    exit(EXIT_SUCCESS);
}

// to debug:
// lldb -- metapop arg1 arg2 etc...
// run

// lldb -- ./metapop beta 8.0 maxvaccprob 0.4 minvaccprob 0.1 timesteps 20 patchpop 100000.0 birthrate 0.02 popstddev 0.1 images 1 iter 1 stochastic 1

// ./metapop beta 2.0 maxvaccprob 0.2 minvaccprob 0.1 timesteps 4 patchpop 10000 birthrate 40.0 popstddev 0.0 images 0 logs 0 fn afn9 iter 1 stochastic 1  > thefirst.txt
