#pragma once

#include "leptoproto.h"

#include "libjsondata.h"

#include <iostream>
#include <string>
#include <fstream>
#include <thread>
#include <chrono>
#include <utility>

#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Dense>
#include <algorithm>
#include <map>

using namespace std;
using Eigen::MatrixXd;
using namespace Eigen;

using namespace mty::lib;


namespace leptoproto{

#define mpl 1.22e19
#define hc2 3.89379338e-28
#define hc2b 0.3893793656e9
#define hbar 6.582119514e-25
#define c 29979245800
#define sbyrho 2.743855e08

#define MAXSTEPS 3000

#define EP 0.01


const int ki = 12;
const int li = 13;
const double xki[ki] = {
0.2049481473588502838963231E-04,
0.3116434971088132598181122E-02,
0.3567129378808611226048581E-01,
0.1637038106438552243872225E+00,
0.4810226336846232187336123E+00,
0.1093749789758291561354613E+01,
0.2113272551481096111677349E+01,
0.3658647547020350981572865E+01,
0.5870258102208375884709533E+01,
0.8943499910675444211687219E+01,
0.1321831883068720803236893E+02,
0.1952704436176789300899542E+02};

const double aki[li] = {
0.1664210817022919684179579E-03,
0.1002537753936819105781108E-01,
0.6402700468555379941006606E-01,
0.1745721822360148306828806E+00,
0.2762012799800637007904487E+00,
0.2670162411605854677879026E+00,
0.1524436358317874768111974E+00,
0.4772451419171662163940487E-01,
0.7344423584282093647907940E-02,
0.4696806037099412766549128E-03,
0.9212653029824790631065364E-05,
0.2645218576013667721828686E-07};

const double xs[li] ={
    -1.,
    -0.8333333333333334,
    -0.6666666666666666,
    -0.5,
    -0.3333333333333333,
    -0.16666666666666666,
    0,
    0.16666666666666666,
    0.3333333333333333,
    0.5,
    0.6666666666666666,
    0.8333333333333334,
    1.};

const double delta = 2./36.;
const double as[li] = {
    1.*delta,
    4.*delta,
    2.*delta,
    4.*delta,
    2.*delta,
    4.*delta,
    2.*delta,
    4.*delta,
    2.*delta,
    4.*delta,
    2.*delta,
    4.*delta,
    1.*delta
};

const double delta1 = 1./36.;
const double as1[li] = {
    1.*delta1,
    4.*delta1,
    2.*delta1,
    4.*delta1,
    2.*delta1,
    4.*delta1,
    2.*delta1,
    4.*delta1,
    2.*delta1,
    4.*delta1,
    2.*delta1,
    4.*delta1,
    1.*delta1
};

const double xs1[li] = {
    0,
    0.08333333333333333,
    0.16666666666666666,
    0.25,
    0.3333333333333333,
    0.4166666666666667,
    0.5,
    0.5833333333333334,
    0.6666666666666666,
    0.75,
    0.8333333333333334,
    0.9166666666666666,
    1.};

vector<double> const temp {0.001,0.00107855,0.00124135,0.00146572,0.00173064,0.00209632,0.00257195,0.00299836,0.00387152,0.00499896,0.00613318,0.00781887,0.00923186,0.0117695,0.0138967,0.0168333,0.0196224,0.0243828,0.0295331,0.0357713,0.0401299,0.0461851,0.05594,0.0660461,0.0760109,0.0932506,0.100678,0.111505,0.123499,0.135044,0.142107,0.145769,0.157366,0.163491,0.165567,0.178733,0.185666,0.185618,0.195309,0.200319,0.205471,0.218983,0.227483,0.230332,0.239298,0.255056,0.271857,0.289765,0.31283,0.333425,0.355377,0.383645,0.408908,0.430309,0.452826,0.476503,0.495087,0.534513,0.555319,0.576962,0.599464,0.663911,0.698601,0.754189,0.7937,0.83525,0.890261,0.998665,1.06445,1.14919,1.22492,1.39183,1.48357,1.66428,1.84335,2.14871,2.67002,3.19304,3.58214,4.22949,4.86791,5.74775,7.42174,8.76275,10.4791,14.4232,16.8126,19.1031,23.7358,27.666,31.0359,33.5058,36.6385,40.0631,43.2509,47.2936,52.3792,58.0118,62.6271,68.4816,74.8861,82.9409,90.6965,101.741,115.604,133.046,141.818,305.292,202.821,258.559,409.624,496.176,578.424,709.656,827.28,928.114,1000};

vector<double> const geff {10.5372,10.5372,10.5384,10.5398,10.5412,10.6599,10.6617,10.3117,10.548,10.7844,10.7861,10.7881,11.0237,10.7916,10.793,10.7946,11.2642,11.8514,12.5555,13.2596,13.7288,14.1983,15.0195,15.7233,16.3099,17.014,17.483,18.4204,19.2408,20.1781,21.5834,22.7544,24.0428,25.5651,27.2042,28.8439,31.654,34.2296,36.4544,38.7961,40.4353,42.5432,45.1191,48.3973,49.9195,51.2079,52.2621,53.3163,54.2535,55.6589,57.0644,58.4699,59.7582,60.8123,61.9835,63.5059,64.3257,64.9117,66.434,67.488,68.3078,69.5965,71.4701,72.6414,73.2272,74.1643,75.3355,76.5072,77.5614,78.3816,79.2016,79.7881,80.4911,81.3116,81.7807,82.6015,83.0717,83.4244,83.7766,84.0121,84.0133,84.0147,84.0169,84.4866,84.9564,85.5445,86.2482,87.0688,88.3585,89.7646,90.5851,91.6394,92.4597,93.5141,94.6855,95.7399,96.7944,97.849,99.1374,100.075,100.661,101.481,102.184,103.239,103.943,104.529,104.998,105.941,105.703,105.94,105.709,105.828,105.712,105.831,105.832,105.95,105.95};

struct plist{
    std::vector<int> procid;
};

void initialize(param_t &pp);
void findAmp(std::string pname,int &id);
double geffT(double x);
double Yeq(std::string name,double T);
double Massp(std::string name,param_t &pp);
double jac(complex_t sqss,Eigen::VectorXd &mm,vector<double> const &rr,int const &intp,int const &outp);
std::vector<std::vector<complex_t>> pmom_std(
        complex_t sqss,std::vector<double> &mm,vector<double> const &cthz,vector<double> const &rr,int const &inp, int const &outp);
Eigen::MatrixXcd pmom(complex_t sqss,Eigen::VectorXd &mm,vector<double> const &cthz,vector<double> const &rr,int const &inp, int const &outp);
Eigen::MatrixXcd Mass(Eigen::VectorXcd &th,Eigen::VectorXcd &mm);
Eigen::MatrixXcd expM(Eigen::MatrixXd M,int &dim);
complex_t L(complex_t xx,complex_t yy,complex_t zz);

complex_t Decay(int &pid,param_t &pp,Process prc,int loop);
//double gamma(int &pid,double const &T,param_t &data,Process prc);
double gamma(int &pid,double const &T,param_t &data,Process prc,int loop);
double epsilon1(std::string ll,std::string ptl,param_t &pp);
void initialize(std::vector<Process> &proc);
void BESolver(std::ofstream &ofile,param_t &pp);
//double bisection(double a, double b,double z_cut,param_t &pp);
double bisection(double a, double b,double z_cut,param_t &pp,std::string name);
double washout(double zz,param_t &pp,std::string name);
//double washout(double zz,param_t &pp);
double freezeout(param_t &pp);

}
