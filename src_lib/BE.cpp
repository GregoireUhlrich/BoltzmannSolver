#include "BE.h"
#include <array>



using namespace std;

using Eigen::MatrixXd;
using namespace Eigen;


using namespace mty::lib;


namespace leptoproto{

std::map<std::string, double> ptlmap;

vector<std::string> Xptl;

vector<std::string> SMptl;

std::map<std::string, double> gptl;

std::vector<Process> procL;
std::vector<Process> procQ;

Eigen::MatrixXi GamQa, GamQb, GamLa, GamLb;
Eigen::VectorXi DecayQ, DecayL;

double Mscale;
//=====================================================================

complex_t Complex(int const &x,double const &y)
{
    return x + y*1i;
}

double geffT(double x)
{
    int size = temp.size();
    
    int i = 0;                  // find left end of interval for interpolation
    if ( x >= temp[size - 2] )  // special case: beyond right end
    {
        i = size - 2;
    }
    else
    {
        while ( x > temp[i+1] ) i++;
    }
    double xL = temp[i], yL = geff[i], xR = geff[i+1], yR = geff[i+1];      // points on either side (unless beyond ends)
    if ( x < xL ) yR = yL;
    if ( x > xR ) yL = yR;
    
    
    double dydx = ( yR - yL ) / ( xR - xL );       // gradient
    
    return yL + dydx * ( x - xL );                 // linear interpolation
}

double Massp(std::string name,param_t &pp)
{
    double mi;
    std::string ptl_mass = "m_";
    
    if (pp.realParams.find(ptl_mass + name) == pp.realParams.end())
    { // no mass,
        mi = 0.;
    }
    else{
        auto pt_mass = pp.realParams[ptl_mass + name];
        mi = *pt_mass;
        
    }
    
    return mi;
}

double Yeq(std::string name,double T)
{
    double yy = 0.;
    double mi = ptlmap[name];//Massp(name,pp);
    double gi = 1.;
    gi *= gptl[name];
    if(mi>0.)
    {
        yy = 45.*gi/geffT(T)*pow(mi/(2.*M_PI*T),2.)*cyl_bessel_k(2, mi/T);
    }
    if(mi==0.)
    {
        yy = 45./2.*gi/geffT(T)*pow(M_PI,-4.); //cout << "h1!!" << endl;
    }
    return yy;
}

Eigen::MatrixXcd expM(Eigen::MatrixXd M,int &dim)
{
    Eigen::MatrixXcd res,ev;
    Eigen::EigenSolver<MatrixXd> es(M);
    ev = Eigen::MatrixXcd::Zero(dim,dim);
    int i;
    for(i=0; i<dim; i++)
    {
        ev(i,i) = exp(es.eigenvalues()[i]);
    }
    Eigen::MatrixXcd V = es.eigenvectors();
    res = V*ev*V.inverse();
    return res;
}


void initialize(param_t &pp)
{
    ParticleData pdata;
    
    pdata.loadFile("script/test.json");
    
    procL = pdata.getProcessesFromQNumberViolation("L");
    procQ = pdata.getProcessesFromQNumberConservation("Q");
    
    Xptl.clear();
    SMptl.clear();
    
    multimap< double, std::string> tmp_map;
    
    std::vector<std::string> ptlth = pdata.getParticleNames();
    for (const auto &name : ptlth) {
        if(pdata.getQuantumNumberValue("Th", name).getValue()>0)
        {
            //Xptl.push_back(name);
            tmp_map.insert({ Massp(name,pp),name});
        }
        if(pdata.getQuantumNumberValue("Th", name).getValue()<1)
        {
            SMptl.push_back(name);
        }
        ptlmap.insert({ name, Massp(name,pp) });
    }
    
    for(auto it = tmp_map.rbegin(); it != tmp_map.rend(); ++it)Xptl.push_back(it->second);
    
    for (const auto &name : ptlth) gptl.insert({ name,  pdata.getQuantumNumberValue("gg", name).getValue()});

    vector<double> mass_scale;
    
    for(auto &name : Xptl) mass_scale.push_back(ptlmap[name]);
    
    Mscale = *max_element(mass_scale.begin(), mass_scale.end());
    
    int nn = Xptl.size();
    GamQa = Eigen::MatrixXi::Zero(nn,nn);
    GamQb = Eigen::MatrixXi::Zero(nn,nn);
    GamLa = Eigen::MatrixXi::Zero(nn,nn);
    GamLb = Eigen::MatrixXi::Zero(nn,nn);
    
    std::cout << "This is matrix GamQ :" << std::endl << GamQa << std::endl;
    
    double epoch1 = Mscale*0.1, epoch2 = Mscale, epoch = 10.*Mscale;
    /*
    multimap< double, int> scatta,scattb;
    for(size_t ii = 0;ii!=Xptl.size();++ii){
        for(size_t jj = 0;jj!=Xptl.size();++jj){
            
            std::cout << "X["<<ii<<"]["<<jj<<"]" << std::endl;
            for(const auto &prname : procQ){
                int N_proc = prname.inParticles.size();
                if(N_proc>1){
                    if((prname.inParticles[0].name.find(Xptl[ii]) != std::string::npos)&&(prname.inParticles[1].name.find(Xptl[jj] ) != std::string::npos))
                    {
                        int id;
                        findAmp(prname.name,id);
                        double Xsec = 0.;
                        double iter = 0;
                        do{
                            double TT = Mscale*pow(10.,-1. + iter);
                            double zz = pow(10.,-1. + iter);
                            double Hb = sqrt(4.*pow(M_PI,3.)*geffT(TT)/45.)*pow(TT,2.)/mpl;
                            
                            double s = geffT(TT)*2*pow(M_PI,2.)/45.*pow(TT,3.);
                            
                            complex<double> pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                            Xsec += (!isnan(pr.real())&&(pr.real()>0.)) ? pr.real() : 0.;
                            iter = iter + 1.;
                        }while(iter<3.);
                        scatta.insert({ Xsec,id});
                        std::cout << "Xsec : " << Xsec << " id: " << id << " N_proc = " << N_proc << std::endl;
                        
                    }
                    
                    if((prname.inParticles[0].name.find(Xptl[ii]) != std::string::npos)&&(prname.outParticles[1].name.find(Xptl[jj] ) != std::string::npos))
                    {
                        int id;
                        findAmp(prname.name,id);
                        double Xsec = 0.;
                        double iter = 0;
                        do{
                            double TT = Mscale*pow(10.,-1. + iter);
                            double zz = pow(10.,-1. + iter);
                            double Hb = sqrt(4.*pow(M_PI,3.)*geffT(TT)/45.)*pow(TT,2.)/mpl;
                            
                            double s = geffT(TT)*2*pow(M_PI,2.)/45.*pow(TT,3.);
                            
                            complex<double> pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                            Xsec += (!isnan(pr.real())&&(pr.real()>0.)) ? pr.real() : 0.;
                            iter = iter + 1.;
                        }while(iter<3.);
                        scattb.insert({ Xsec,id});
                        std::cout << "Xsec : " << Xsec << " id: " << id << " N_proc = " << N_proc << std::endl;
                        std::cout << " Process : " << prname.name << std::endl;
                    }
                    
                }
                
                
            }
            
            if(!scatta.empty()){
                std::cout << " Answer: " << scatta.rbegin()->second << std::endl;
                GamQa(ii,jj) = scatta.rbegin()->second;
                std::cout << "GamQa["<<ii<<"]["<<jj<<"] = " << GamQa(ii,jj) << std::endl;
                scatta.clear();
            }
            if(!scattb.empty()){
                std::cout << " Answer: " << scattb.rbegin()->second << std::endl;
                GamQb(ii,jj) = scattb.rbegin()->second;
                std::cout << "GamQb["<<ii<<"]["<<jj<<"] = " << GamQb(ii,jj) << std::endl;
                scattb.clear();
            }
            
        }
    }
    std::cout << "Now, This is matrix GamQ :" << std::endl << GamQa << std::endl;
    */
}

Eigen::MatrixXcd Mass(Eigen::VectorXcd &th,Eigen::VectorXcd &mm)
{
    Eigen::MatrixXcd res;
    Eigen::MatrixXcd O12(3,3);
    Eigen::MatrixXcd O13(3,3);
    Eigen::MatrixXcd O23(3,3);
    Eigen::MatrixXcd MM(3,3);
    
    cout << endl << th << endl;
    cout << endl << mm << endl;
    Eigen::MatrixXcd temp(3,3);
    
    O12 << cos(th(0)),-sin(th(0)),0.,
    sin(th(0)),cos(th(0)),0.,
    0.,0.,1.;
    O13 << cos(th(1)),0.,-sin(th(1)),
    0.,1.,0.,
    sin(th(1)),0.,cos(th(1));
    O23 << 1.,0.,0.,
    0.,cos(th(2)),-sin(th(2)),
    0.,sin(th(2)),cos(th(2));
    
    cout << endl << O12 << endl;
    cout << endl << O13 << endl;
    cout << endl << O23 << endl;
    MM << mm(0), 0., 0.,0.,mm(1),0.,0.,0.,mm(2);
    
    temp = O12 * O13 * O23;
    
    
    res = (temp.conjugate()).transpose() * MM * temp;
    
    return res;
}


complex_t L(complex_t xx,complex_t yy,complex_t zz)
{
    return xx*xx + yy*yy + zz*zz - 2.*xx*yy - 2.*xx*zz - 2.*yy*zz;
}

double jac(complex_t sqss,Eigen::VectorXd &mm,vector<double> const &rr,int const &intp,int const &outp)
{
    double res = 0.;
    if((intp==1)&&(outp==3)){
        double mD = mm(0);
        double m1 = mm(1);
        double m2 = mm(2);
        double m3 = mm(3);
        double s23 = pow(m2 + m3 + (-m1 - m2 - m3 + mD)*rr[0],2.);
        
        res = ((mD - m1 - m2 - m3) > 0.) ? 0.5/mD*1./(4.*M_PI)*sqrt(L(mD*mD,m1*m1,s23).real())/(8.*M_PI*mD*mD)*sqrt(L(s23,m2*m2,m3*m3)).real()/(8.*M_PI*s23) : 0.;
    }
    if((intp==1)&&(outp==2)){
        double mD = mm(0);
        double m1 = mm(1);
        double m2 = mm(2);
        
        res = ((mD - m1 - m2) > 0.) ? 1./(16.*M_PI*pow(mD,3.))*sqrt(L(mD*mD,m1*m1,m2*m2)).real() : 0.;
    }
    if((intp==2)&&(outp==2)){
        complex_t ss = pow(sqss,2.);
        double ma = mm(0);
        double mb = mm(1);
        double m2 = mm(2);
        double m3 = mm(3);
        
        res = 1./(64.*M_PI*M_PI*ss.real())*sqrt(L(ss,m2*m2,m3*m3)).real()/sqrt(L(ss,ma*ma,mb*mb)).real();
    }
    return res;
}

template<class Type>
Type product(array<Type, 3> const &p, array<Type, 3> const &q)
{
    return p[0]*q[0] + p[1]*q[1] + p[2]*q[2];
}

vector<vector<complex_t>> pmom_std(
        complex_t sqss,
        vector<double> &mm,
        vector<double> const &cthz,
        vector<double> const &rr,
        int const &inp, 
        int const &outp)
{
    vector<vector<complex_t>> res(inp+outp);
    for (auto &row : res) {
        row = vector<complex_t>(inp+outp, 0);
    }
    
    array<array<complex_t, 3>, 4> pmom = {0};
    
    if((inp==1)&&(outp==3))
    {
        double mD = mm[0];
        double m1 = mm[1];
        double m2 = mm[2];
        double m3 = mm[3];
        double cth1 = cthz[0];
        double cth2 = cthz[1];
        
        double s23 = pow(m2 + m3 + (-m1 - m2 - m3 + mD)*rr[0],2.);
        
        pmom[0] = {mD, 0.,0.};
        pmom[1] = {(sqrt(pow(mD,2))*(1 + pow(m1,2)/pow(mD,2) - s23/pow(mD,2)))/2.,
        (Complex(0,-0.5)*sqrt(1 - pow(cth1,2))*sqrt(L(pow(mD,2),pow(m1,2),s23)))/sqrt(pow(mD,2)),
        (Complex(0,-0.5)*cth1*sqrt(L(pow(mD,2),pow(m1,2),s23)))/sqrt(pow(mD,2))};
        pmom[2] = {(-((pow(m1,2) - pow(mD,2) - s23)*(pow(m2,2) - pow(m3,2) + s23)) +
                    (cth1*cth2 + sqrt(1 - pow(cth1,2))*sqrt(1 - pow(cth2,2)))*sqrt(L(pow(mD,2),pow(m1,2),s23))*sqrt(L(s23,pow(m2,2),pow(m3,2))))/
        (4.*sqrt(pow(mD,2))*s23),(Complex(0,0.25)*(sqrt(1 - pow(cth1,2))*(pow(m2,2) - pow(m3,2) + s23)*sqrt(L(pow(mD,2),pow(m1,2),s23)) +
                                                   (pow(cth1,2)*sqrt(1 - pow(cth2,2))*(pow(m1,2) - pow(mD,2) + 2*sqrt(pow(mD,2))*sqrt(s23) - s23) +
                                                    sqrt(1 - pow(cth2,2))*(-pow(m1,2) + pow(mD,2) + s23) +
                                                    cth1*sqrt(1 - pow(cth1,2))*cth2*(-pow(m1,2) + pow(mD,2) - 2*sqrt(pow(mD,2))*sqrt(s23) + s23))*sqrt(L(s23,pow(m2,2),pow(m3,2)))
                                                   ))/(sqrt(pow(mD,2))*s23),(Complex(0,0.25)*(cth1*(pow(m2,2) - pow(m3,2) + s23)*sqrt(L(pow(mD,2),pow(m1,2),s23)) +
                                                                                              (2*cth2*sqrt(pow(mD,2))*sqrt(s23) + pow(cth1,2)*cth2*(-pow(m1,2) + pow(mD,2) - 2*sqrt(pow(mD,2))*sqrt(s23) + s23) +
                                                                                               cth1*sqrt(1 - pow(cth1,2))*sqrt(1 - pow(cth2,2))*(-pow(m1,2) + pow(mD,2) - 2*sqrt(pow(mD,2))*sqrt(s23) + s23))*
                                                                                              sqrt(L(s23,pow(m2,2),pow(m3,2)))))/(sqrt(pow(mD,2))*s23)};
        pmom[3] = {(-((pow(m2,2) - pow(m3,2) - s23)*(-pow(m1,2) + pow(mD,2) + s23)) -
                    (cth1*cth2 + sqrt(1 - pow(cth1,2))*sqrt(1 - pow(cth2,2)))*sqrt(L(pow(mD,2),pow(m1,2),s23))*sqrt(L(s23,pow(m2,2),pow(m3,2))))/
        (4.*sqrt(pow(mD,2))*s23),(Complex(0,0.25)*(sqrt(1 - pow(cth1,2))*(-pow(m2,2) + pow(m3,2) + s23)*sqrt(L(pow(mD,2),pow(m1,2),s23)) -
                                                   (pow(cth1,2)*sqrt(1 - pow(cth2,2))*(pow(m1,2) - pow(mD,2) + 2*sqrt(pow(mD,2))*sqrt(s23) - s23) +
                                                    sqrt(1 - pow(cth2,2))*(-pow(m1,2) + pow(mD,2) + s23) +
                                                    cth1*sqrt(1 - pow(cth1,2))*cth2*(-pow(m1,2) + pow(mD,2) - 2*sqrt(pow(mD,2))*sqrt(s23) + s23))*sqrt(L(s23,pow(m2,2),pow(m3,2)))
                                                   ))/(sqrt(pow(mD,2))*s23),(Complex(0,0.25)*(cth1*(-pow(m2,2) + pow(m3,2) + s23)*sqrt(L(pow(mD,2),pow(m1,2),s23)) +
                                                                                              (pow(cth1,2)*cth2*(pow(m1,2) - pow(mD,2) + 2*sqrt(pow(mD,2))*sqrt(s23) - s23) +
                                                                                               cth1*sqrt(1 - pow(cth1,2))*sqrt(1 - pow(cth2,2))*(pow(m1,2) - pow(mD,2) + 2*sqrt(pow(mD,2))*sqrt(s23) - s23) -
                                                                                               2*cth2*sqrt(pow(mD,2))*sqrt(s23))*sqrt(L(s23,pow(m2,2),pow(m3,2)))))/(sqrt(pow(mD,2))*s23)};
        
        res[0] = {product(pmom[0], pmom[0]), product(pmom[0], pmom[1]), product(pmom[0], pmom[2]), product(pmom[0], pmom[3])};
        res[1] = {product(pmom[1], pmom[0]), product(pmom[1], pmom[1]), product(pmom[1], pmom[2]), product(pmom[1], pmom[3])};
        res[2] = {product(pmom[2], pmom[0]), product(pmom[2], pmom[1]), product(pmom[2], pmom[2]), product(pmom[2], pmom[3])};
        res[3] = {product(pmom[3], pmom[0]), product(pmom[3], pmom[1]), product(pmom[3], pmom[2]), product(pmom[3], pmom[3])};
    }
    if((inp==1)&&(outp==2)) {
        double mD = mm[0];
        double m1 = mm[1];
        double m2 = mm[2];
        
        pmom[0] = {mD, 0.,0.};
        pmom[1] = {(pow(m1,2) - pow(m2,2) + pow(mD,2))/(2.*sqrt(pow(mD,2))),0,
        -0.5*(sqrt(L(pow(mD,2),pow(m1,2),pow(m2,2))))/sqrt(pow(mD,2))*1i};
        pmom[2] = {(pow(m1,2) - pow(m2,2) + pow(mD,2))/(2.*sqrt(pow(mD,2))),0,
        0.5*(sqrt(L(pow(mD,2),pow(m1,2),pow(m2,2))))/sqrt(pow(mD,2))*1i};
        
        res[0] = {product(pmom[0], pmom[0]), product(pmom[0], pmom[1]), product(pmom[0], pmom[2])};
        res[1] = {product(pmom[1], pmom[0]), product(pmom[1], pmom[1]), product(pmom[1], pmom[2])};
        res[2] = {product(pmom[2], pmom[0]), product(pmom[2], pmom[1]), product(pmom[2], pmom[2])};
    }
    if((inp==2)&&(outp==2)) {
        double ss = pow(sqss.real(),2.);
        double m1 = mm[0];
        double m2 = mm[1];
        double m3 = mm[2];
        double m4 = mm[3];
        double cth = cthz[0];
        
        pmom[0] = {((1 + pow(m1,2.)/ss - pow(m2,2.)/ss)*sqrt(ss))/2.,0,-sqrt(L(ss,pow(m1,2.),pow(m2,2.)))/(2.*sqrt(ss))*1i};
        pmom[1] = {((1 - pow(m1,2.)/ss + pow(m2,2.)/ss)*sqrt(ss))/2.,0,sqrt(L(ss,pow(m1,2),pow(m2,2.)))/(2.*sqrt(ss))*1i};
        pmom[2] = {((1 + pow(m3,2)/ss - pow(m4,2)/ss)*sqrt(ss))/2.,-0.5*(sqrt(1 - pow(cth,2))*sqrt(L(ss,pow(m3,2.),pow(m4,2.))))/sqrt(ss)*1i,
        -0.5*(cth*sqrt(L(ss,pow(m3,2.),pow(m4,2.))))/sqrt(ss)*1i};
        pmom[3] = {((1 - pow(m3,2.)/ss + pow(m4,2.)/ss)*sqrt(ss))/2.,0.5*(sqrt(1 - pow(cth,2.))*sqrt(L(ss,pow(m3,2.),pow(m4,2.))))/sqrt(ss)*1i,
        0.5*(cth*sqrt(L(ss,pow(m3,2.),pow(m4,2.))))/sqrt(ss)*1i};
        
        res[0] = {product(pmom[0], pmom[0]), product(pmom[0], pmom[1]), product(pmom[0], pmom[2]), product(pmom[0], pmom[3])};
        res[1] = {product(pmom[1], pmom[0]), product(pmom[1], pmom[1]), product(pmom[1], pmom[2]), product(pmom[1], pmom[3])};
        res[2] = {product(pmom[2], pmom[0]), product(pmom[2], pmom[1]), product(pmom[2], pmom[2]), product(pmom[2], pmom[3])};
        res[3] = {product(pmom[3], pmom[0]), product(pmom[3], pmom[1]), product(pmom[3], pmom[2]), product(pmom[3], pmom[3])};
    }
    return res;
}

Eigen::MatrixXcd pmom(
        complex_t sqss,
        Eigen::VectorXd &mm,
        vector<double> const &cthz,
        vector<double> const &rr,
        int const &inp, 
        int const &outp)
{
    Eigen::MatrixXcd res(inp+outp,inp+outp);
    
    std::array<Eigen::Vector3cd,4> pmom;
    
    for (int ii=0; ii!=4; ++ii){pmom[ii] = Eigen::VectorXcd::Zero(3);}
    
    if((inp==1)&&(outp==3))
    {
        double mD = mm(0);
        double m1 = mm(1);
        double m2 = mm(2);
        double m3 = mm(3);
        double cth1 = cthz[0];
        double cth2 = cthz[1];
        
        double s23 = pow(m2 + m3 + (-m1 - m2 - m3 + mD)*rr[0],2.);
        
        res = Eigen::MatrixXcd::Zero(4,4);
        
        pmom[0] << mD, 0.,0.;
        pmom[1] << (sqrt(pow(mD,2))*(1 + pow(m1,2)/pow(mD,2) - s23/pow(mD,2)))/2.,
        (Complex(0,-0.5)*sqrt(1 - pow(cth1,2))*sqrt(L(pow(mD,2),pow(m1,2),s23)))/sqrt(pow(mD,2)),
        (Complex(0,-0.5)*cth1*sqrt(L(pow(mD,2),pow(m1,2),s23)))/sqrt(pow(mD,2));
        pmom[2] << (-((pow(m1,2) - pow(mD,2) - s23)*(pow(m2,2) - pow(m3,2) + s23)) +
                    (cth1*cth2 + sqrt(1 - pow(cth1,2))*sqrt(1 - pow(cth2,2)))*sqrt(L(pow(mD,2),pow(m1,2),s23))*sqrt(L(s23,pow(m2,2),pow(m3,2))))/
        (4.*sqrt(pow(mD,2))*s23),(Complex(0,0.25)*(sqrt(1 - pow(cth1,2))*(pow(m2,2) - pow(m3,2) + s23)*sqrt(L(pow(mD,2),pow(m1,2),s23)) +
                                                   (pow(cth1,2)*sqrt(1 - pow(cth2,2))*(pow(m1,2) - pow(mD,2) + 2*sqrt(pow(mD,2))*sqrt(s23) - s23) +
                                                    sqrt(1 - pow(cth2,2))*(-pow(m1,2) + pow(mD,2) + s23) +
                                                    cth1*sqrt(1 - pow(cth1,2))*cth2*(-pow(m1,2) + pow(mD,2) - 2*sqrt(pow(mD,2))*sqrt(s23) + s23))*sqrt(L(s23,pow(m2,2),pow(m3,2)))
                                                   ))/(sqrt(pow(mD,2))*s23),(Complex(0,0.25)*(cth1*(pow(m2,2) - pow(m3,2) + s23)*sqrt(L(pow(mD,2),pow(m1,2),s23)) +
                                                                                              (2*cth2*sqrt(pow(mD,2))*sqrt(s23) + pow(cth1,2)*cth2*(-pow(m1,2) + pow(mD,2) - 2*sqrt(pow(mD,2))*sqrt(s23) + s23) +
                                                                                               cth1*sqrt(1 - pow(cth1,2))*sqrt(1 - pow(cth2,2))*(-pow(m1,2) + pow(mD,2) - 2*sqrt(pow(mD,2))*sqrt(s23) + s23))*
                                                                                              sqrt(L(s23,pow(m2,2),pow(m3,2)))))/(sqrt(pow(mD,2))*s23);
        pmom[3] << (-((pow(m2,2) - pow(m3,2) - s23)*(-pow(m1,2) + pow(mD,2) + s23)) -
                    (cth1*cth2 + sqrt(1 - pow(cth1,2))*sqrt(1 - pow(cth2,2)))*sqrt(L(pow(mD,2),pow(m1,2),s23))*sqrt(L(s23,pow(m2,2),pow(m3,2))))/
        (4.*sqrt(pow(mD,2))*s23),(Complex(0,0.25)*(sqrt(1 - pow(cth1,2))*(-pow(m2,2) + pow(m3,2) + s23)*sqrt(L(pow(mD,2),pow(m1,2),s23)) -
                                                   (pow(cth1,2)*sqrt(1 - pow(cth2,2))*(pow(m1,2) - pow(mD,2) + 2*sqrt(pow(mD,2))*sqrt(s23) - s23) +
                                                    sqrt(1 - pow(cth2,2))*(-pow(m1,2) + pow(mD,2) + s23) +
                                                    cth1*sqrt(1 - pow(cth1,2))*cth2*(-pow(m1,2) + pow(mD,2) - 2*sqrt(pow(mD,2))*sqrt(s23) + s23))*sqrt(L(s23,pow(m2,2),pow(m3,2)))
                                                   ))/(sqrt(pow(mD,2))*s23),(Complex(0,0.25)*(cth1*(-pow(m2,2) + pow(m3,2) + s23)*sqrt(L(pow(mD,2),pow(m1,2),s23)) +
                                                                                              (pow(cth1,2)*cth2*(pow(m1,2) - pow(mD,2) + 2*sqrt(pow(mD,2))*sqrt(s23) - s23) +
                                                                                               cth1*sqrt(1 - pow(cth1,2))*sqrt(1 - pow(cth2,2))*(pow(m1,2) - pow(mD,2) + 2*sqrt(pow(mD,2))*sqrt(s23) - s23) -
                                                                                               2*cth2*sqrt(pow(mD,2))*sqrt(s23))*sqrt(L(s23,pow(m2,2),pow(m3,2)))))/(sqrt(pow(mD,2))*s23);
        
        res << pmom[0].transpose() * pmom[0], pmom[0].transpose() * pmom[1],      pmom[0].transpose() * pmom[2], pmom[0].transpose() * pmom[3],
        pmom[1].transpose() * pmom[0], pmom[1].transpose() * pmom[1], pmom[1].transpose() * pmom[2], pmom[1].transpose() * pmom[3],
        pmom[2].transpose() * pmom[0], pmom[2].transpose() * pmom[1], pmom[2].transpose() * pmom[2], pmom[2].transpose() * pmom[3],
        pmom[3].transpose() * pmom[0], pmom[3].transpose() * pmom[1], pmom[3].transpose() * pmom[2], pmom[3].transpose() * pmom[3];
    }
    if((inp==1)&&(outp==2)) {
        double mD = mm(0);
        double m1 = mm(1);
        double m2 = mm(2);
        
        pmom[0] << mD, 0.,0.;
        pmom[1] << (pow(m1,2) - pow(m2,2) + pow(mD,2))/(2.*sqrt(pow(mD,2))),0,
        -0.5*(sqrt(L(pow(mD,2),pow(m1,2),pow(m2,2))))/sqrt(pow(mD,2))*1i;
        pmom[2] << (pow(m1,2) - pow(m2,2) + pow(mD,2))/(2.*sqrt(pow(mD,2))),0,
        0.5*(sqrt(L(pow(mD,2),pow(m1,2),pow(m2,2))))/sqrt(pow(mD,2))*1i;
        
        res = Eigen::MatrixXcd::Zero(3,3);
        
        res << pmom[0].transpose() * pmom[0], pmom[0].transpose() * pmom[1],      pmom[0].transpose() * pmom[2],
        pmom[1].transpose() * pmom[0], pmom[1].transpose() * pmom[1], pmom[1].transpose() * pmom[2],
        pmom[2].transpose() * pmom[0], pmom[2].transpose() * pmom[1], pmom[2].transpose() * pmom[2];
    }
    if((inp==2)&&(outp==2)) {
        double ss = pow(sqss.real(),2.);
        double m1 = mm(0);
        double m2 = mm(1);
        double m3 = mm(2);
        double m4 = mm(3);
        double cth = cthz[0];
        
        pmom[0] << ((1 + pow(m1,2.)/ss - pow(m2,2.)/ss)*sqrt(ss))/2.,0,-sqrt(L(ss,pow(m1,2.),pow(m2,2.)))/(2.*sqrt(ss))*1i;
        pmom[1] << ((1 - pow(m1,2.)/ss + pow(m2,2.)/ss)*sqrt(ss))/2.,0,sqrt(L(ss,pow(m1,2),pow(m2,2.)))/(2.*sqrt(ss))*1i;
        pmom[2] << ((1 + pow(m3,2)/ss - pow(m4,2)/ss)*sqrt(ss))/2.,-0.5*(sqrt(1 - pow(cth,2))*sqrt(L(ss,pow(m3,2.),pow(m4,2.))))/sqrt(ss)*1i,
        -0.5*(cth*sqrt(L(ss,pow(m3,2.),pow(m4,2.))))/sqrt(ss)*1i;
        pmom[3] << ((1 - pow(m3,2.)/ss + pow(m4,2.)/ss)*sqrt(ss))/2.,0.5*(sqrt(1 - pow(cth,2.))*sqrt(L(ss,pow(m3,2.),pow(m4,2.))))/sqrt(ss)*1i,
        0.5*(cth*sqrt(L(ss,pow(m3,2.),pow(m4,2.))))/sqrt(ss)*1i;
        
        res = Eigen::MatrixXcd::Zero(4,4);
        
        res << pmom[0].transpose() * pmom[0], pmom[0].transpose() * pmom[1],      pmom[0].transpose() * pmom[2], pmom[0].transpose() * pmom[3],
        pmom[1].transpose() * pmom[0], pmom[1].transpose() * pmom[1], pmom[1].transpose() * pmom[2], pmom[1].transpose() * pmom[3],
        pmom[2].transpose() * pmom[0], pmom[2].transpose() * pmom[1], pmom[2].transpose() * pmom[2], pmom[2].transpose() * pmom[3],
        pmom[3].transpose() * pmom[0], pmom[3].transpose() * pmom[1], pmom[3].transpose() * pmom[2], pmom[3].transpose() * pmom[3];
    }
    return res;
}

void findAmp(std::string pname,int &id)
{
    id = f_G.size() + 1;
    for(int ll=0; ll!=f_G.size();++ll)
    {
        if (f_G[ll].name == pname) {
            id = ll;
        }
    }
    
}

double epsilon1(std::string ll,std::string ptl,param_t &pp)
{
    std::vector<Process> procL;
    ParticleData pdata;
    pdata.loadFile("script/test.json");
    procL = pdata.getProcessesFromQNumberViolation(ll);
    std::string ptl_mass = "m_";
    complex<double> ampLp = 0.;
    complex<double> ampLap = 0.;
    
    complex<double> tree = 0.; // Only for testing
    for(const auto &prc : procL){
        int in_proc = prc.inParticles.size();
        int out_proc = prc.outParticles.size();
        double DeltaL=prc.qNumbers.begin()->second;
        if((in_proc<2)&&(out_proc>1)){
            if(prc.inParticles[0].name.find(ptl)!=std::string::npos){
                
                Eigen::MatrixXcd pmat(1+out_proc,1+out_proc);
                
                vector<double> mmax;
                for(const auto &prname : prc.inParticles)mmax.push_back(ptlmap[prname.name]);
                
                for(const auto &prname : prc.outParticles)mmax.push_back(ptlmap[prname.name]);
                
                Eigen::VectorXd mass = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(mmax.data(), mmax.size());
                
                
                vector<double> cth,rr;
                cth.push_back(1.);
                rr.push_back(1.);
                
                pmat = pmom(mass(0),mass,cth,rr,1,out_proc);
                pp.s_12 = pmat(0,1).real();
                pp.s_13 = pmat(0,2).real();
                pp.s_23 = pmat(1,2).real();
                int pid;
                if(prc.name.find("AsymL")!=std::string::npos){
                    findAmp(prc.name,pid);
                    if(DeltaL==1){
                        ampLp += f_G[pid](pp);
                    }
                }
                if(prc.name.find("AsymT")!=std::string::npos){
                    findAmp(prc.name,pid);
                    if(DeltaL==1){
                        tree += f_G[pid](pp);
                    }
                }
                if(prc.name.find("AsymL")!=std::string::npos){
                    findAmp(prc.name,pid);
                    if(DeltaL==-1){
                        ampLap += f_G[pid](pp);
                    }
                }
                if(prc.name.find("AsymT")!=std::string::npos){
                    findAmp(prc.name,pid);
                    if(DeltaL==-1){
                        tree += f_G[pid](pp);
                    }
                }
            }
        }
    }
    
    return 2.*((ampLp - ampLap)).real()/(tree.real());
}


complex_t Decay(int &pid,param_t &pp,Process prc,int loop)
{
    complex_t res = 0.;
    
    std::string pname = f_G[pid].name;
    std::string ptl_mass = "m_";
    int out = prc.outParticles.size();
    
    
    complex_t jacD = 0.,decay = 0.;
    Eigen::MatrixXcd pmat(1+out,1+out);
    
    vector<double> mmax;
    for(const auto &prname : prc.inParticles) mmax.push_back(ptlmap[prname.name]);
    
    for(const auto &prname : prc.outParticles) mmax.push_back(ptlmap[prname.name]);
    
    Eigen::VectorXd mass = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(mmax.data(), mmax.size());
    
    if(out==3){
        decay = 0.;
        for(int ii = 0;ii<13;ii++){
            for(int jj = 0;jj<13;jj++){
                for(int kk = 0;kk<13;kk++){
                    
                    vector<double> cth,rr;
                    cth.push_back(xs[jj]);
                    cth.push_back(xs[kk]);
                    rr.push_back(xs1[ii]);
                    pmat = pmom(mass(0),mass,cth,rr,1,out);
                    pp.s_12 = pmat(0,1).real();
                    pp.s_13 = pmat(0,2).real();
                    pp.s_14 = pmat(0,3).real();
                    pp.s_23 = pmat(1,2).real();
                    pp.s_24 = pmat(1,3).real();
                    pp.s_34 = pmat(2,3).real();
                    
                    jacD = jac(mass(0),mass,rr,1,out);
                    
                    if(loop==0){
                        if(f_G[pid].name.find("Tree")!=std::string::npos){
                            decay += f_G[pid](pp)*jacD*as[jj]*as[kk]*as1[ii];
                        }
                    }
                    if(loop==1){
                        if(f_G[pid].name.find("AsymL")!=std::string::npos){
                            decay += f_G[pid](pp)*jacD*as[jj]*as[kk]*as1[ii];
                        }
                    }
                }
            }
        }
    }
    if(out==2){
        decay = 0.;
        vector<double> cth,rr;
        cth.push_back(1.);
        rr.push_back(1.);
        pmat = pmom(mass(0),mass,cth,rr,1,out);
        pp.s_12 = pmat(0,1).real();
        pp.s_13 = pmat(0,2).real();
        pp.s_23 = pmat(1,2).real();
        
        jacD = jac(mass(0),mass,rr,1,out);
        
        if(loop==0){
            if(f_G[pid].name.find("Tree")!=std::string::npos){
                decay += f_G[pid](pp)*jacD;
            }
        }
        if(loop==1){
            if(f_G[pid].name.find("AsymL")!=std::string::npos){
                decay += f_G[pid](pp)*jacD;
            }
        }
        
        
        
    }
    res = decay;
    return res;
}

double gamma_std(int &pid,double const &T,param_t &data,Process prc,int loop)
{
    double result = 0;
    
    double gg = 1.;
    vector< pair <int,double> > gvect;
    
    std::string pname = f_G[pid].name;
    std::string ptl_mass = "m_";
    
    int in = prc.inParticles.size(),out = prc.outParticles.size();
    
    for(const auto &inpp : prc.inParticles) gg *= gptl[inpp.name];
    
    vector<double> mass;
    for(const auto &prname : prc.inParticles)mass.push_back(ptlmap[prname.name]);
    
    for(const auto &prname : prc.outParticles)mass.push_back(ptlmap[prname.name]);
    
    if((in==2)&&(out==2)){
        
        double s,pin,pout;
        
        double ma,mb,m1,m2,ms,Tsc;
        
        ms = *max_element(mass.begin(), mass.end());
        Tsc = T/ms;
        ma = abs(mass[0]/ms);
        mb = abs(mass[1]/ms);
        m1 = abs(mass[2]/ms);
        m2 = abs(mass[3]/ms);
        int ii,jj;
        double jac = 0.,mi=0.,mj=0.;
        complex_t amp = 0.;
        if((pow(ma+mb, 2.))>=(pow(m1+m2, 2.)))
        {
            mi=ma;
            mj=mb;
            
        }
        if((pow(ma+mb, 2.))<(pow(m1+m2, 2.)))
        {
            mi=m1;
            mj=m2;
        }
        for (ii=0; ii<12; ii++)
        {
            for (jj=0; jj<13; jj++)
            {
                s = pow(mi-mj,2.) + 2.*mi*mj*(2. + xki[ii]);//xs1[ii](1.- xs1[ii]);
                pin = 0.5*sqrt(L(s,ma*ma,mb*mb)/s).real();
                pout = 0.5*sqrt(L(s,m1*m1,m2*m2)/s).real();
                
                vector<double> cth,rr;
                cth.push_back(xs[jj]);
                rr.push_back(1.);
                auto pmat = pmom_std(sqrt(s)*ms,mass,cth,rr,in,out);
                
                data.s_12 = pmat[0][1].real();
                data.s_13 = pmat[0][2].real();
                data.s_14 = pmat[0][3].real();
                data.s_23 = pmat[1][2].real();
                data.s_24 = pmat[1][3].real();
                data.s_34 = pmat[2][3].real();
                
                jac = gg*pow(ms,4.)*2.*mi*mj*aki[ii]*as[jj];//*(1./(pow(xs1[ii] - 1.,2.)));
                amp = 0.;
                if(loop==0){
                    if(f_G[pid].name.find("Tree")!=std::string::npos){
                        amp += f_G[pid](data)*jac;
                    }
                }
                if(loop==1){
                    if(f_G[pid].name.find("AsymL")!=std::string::npos){
                        amp += f_G[pid](data)*jac;
                    }
                }
                result += exp(xki[ii])*Tsc/(256.*pow(M_PI,5.))*pin*pout*amp.real()*1./sqrt(s)*cyl_bessel_k(1, sqrt(s)/Tsc);
                cth.clear();
                rr.clear();
            }
        }
    }
    
    if(in==1)
    {
        complex_t decay = Decay(pid,data,prc,loop);
        
        double mD = mass[0];
        result = gg*mD*mD/(2.*pow(M_PI,2.))*T*cyl_bessel_k(1, abs(mD)/T)*decay.real();
    }
    return result;
}

double gamma(int &pid,double const &T,param_t &data,Process prc,int loop)
{
    static constexpr bool use_std_gamma = true;
    if constexpr (use_std_gamma)
    {
        return gamma_std(pid, T, data, prc, loop);
    }
    double result = 0;
    
    double gg = 1.;
    vector< pair <int,double> > gvect;
    
    std::string pname = f_G[pid].name;
    std::string ptl_mass = "m_";
    
    int in = prc.inParticles.size(),out = prc.outParticles.size();
    
    for(const auto &inpp : prc.inParticles) gg *= gptl[inpp.name];
    
    vector<double> mmax;
    for(const auto &prname : prc.inParticles)mmax.push_back(ptlmap[prname.name]);
    
    for(const auto &prname : prc.outParticles)mmax.push_back(ptlmap[prname.name]);
    
    Eigen::MatrixXcd pmat(in+out,in+out);
    
    Eigen::VectorXd mass = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(mmax.data(), mmax.size());
    
    if((in==2)&&(out==2)){
        
        double s,pin,pout;
        
        double ma,mb,m1,m2,ms,Tsc;
        
        ms = *max_element(mmax.begin(), mmax.end());
        Tsc = T/ms;
        ma = abs(mass(0)/ms);
        mb = abs(mass(1)/ms);
        m1 = abs(mass(2)/ms);
        m2 = abs(mass(3)/ms);
        int ii,jj;
        double jac = 0.,mi=0.,mj=0.;
        complex_t amp = 0.;
        if((pow(ma+mb, 2.))>=(pow(m1+m2, 2.)))
        {
            mi=ma;
            mj=mb;
            
        }
        if((pow(ma+mb, 2.))<(pow(m1+m2, 2.)))
        {
            mi=m1;
            mj=m2;
        }
        for (ii=0; ii<12; ii++)
        {
            for (jj=0; jj<13; jj++)
            {
                s = pow(mi-mj,2.) + 2.*mi*mj*(2. + xki[ii]);//xs1[ii](1.- xs1[ii]);
                pin = 0.5*sqrt(L(s,ma*ma,mb*mb)/s).real();
                pout = 0.5*sqrt(L(s,m1*m1,m2*m2)/s).real();
                
                vector<double> cth,rr;
                cth.push_back(xs[jj]);
                rr.push_back(1.);
                pmat = pmom(sqrt(s)*ms,mass,cth,rr,in,out);
                
                data.s_12 = pmat(0,1).real();
                data.s_13 = pmat(0,2).real();
                data.s_14 = pmat(0,3).real();
                data.s_23 = pmat(1,2).real();
                data.s_24 = pmat(1,3).real();
                data.s_34 = pmat(2,3).real();
                
                jac = gg*pow(ms,4.)*2.*mi*mj*aki[ii]*as[jj];//*(1./(pow(xs1[ii] - 1.,2.)));
                amp = 0.;
                if(loop==0){
                    if(f_G[pid].name.find("Tree")!=std::string::npos){
                        amp += f_G[pid](data)*jac;
                    }
                }
                if(loop==1){
                    if(f_G[pid].name.find("AsymL")!=std::string::npos){
                        amp += f_G[pid](data)*jac;
                    }
                }
                result += exp(xki[ii])*Tsc/(256.*pow(M_PI,5.))*pin*pout*amp.real()*1./sqrt(s)*cyl_bessel_k(1, sqrt(s)/Tsc);
                cth.clear();
                rr.clear();
            }
        }
    }
    
    if(in==1)
    {
        complex_t decay = Decay(pid,data,prc,loop);
        
        double mD = mass(0);
        result = gg*mD*mD/(2.*pow(M_PI,2.))*T*cyl_bessel_k(1, abs(mD)/T)*decay.real();
    }
    return result;
}


double washout(double zz,param_t &pp,std::string name)
{
    double TT = Mscale/zz;
    //std::cout << "TT = " << TT << std::endl;
    double Hb = sqrt(4.*pow(M_PI,3.)*geffT(TT)/45.)*pow(TT,2.)/mpl;
    double sT = geffT(TT)*2*pow(M_PI,2.)/45.;
    double yeql = 4./pow(M_PI,2.)*pow(sT,-1.);
    double s = geffT(TT)*2*pow(M_PI,2.)/45.*pow(TT,3.);
    int id = 0;
    
    double gammaQ = 0.;
    double result = 0.;
    
    if(name.find("L")!= std::string::npos)
    {
        for(const auto &prname : procL)
        {
            findAmp(prname.name,id);
            int DeltaL = prname.qNumbers.begin()->second;
            double gamma_temp = gamma(id,TT,pp,prname,0);
            result += (DeltaL>0) ? DeltaL*gamma_temp/(s*Hb*zz*yeql) : 0.;
        }
    }
    if(name.find("Q")!= std::string::npos)
    {
        for(const auto &prname : procQ){
            findAmp(prname.name,id);
            if((prname.name.find("Tree")!=std::string::npos)){
                double gamma_temp = gamma(id,TT,pp,prname,0);
                gammaQ += gamma_temp;
            }
        }
        result = gammaQ/(s*Hb*zz);
    }
    
    
    
    return result;
    
}

double bisection(double a, double b,double z_cut,param_t &pp,std::string name)
{
    double fa,fb,fc;
    
    fa = washout(pow(10.,a),pp,name) - z_cut;
    
    fb = washout(pow(10.,b),pp,name) - z_cut;
    
    if (fa * fb >= 0) {
        //cout << "You have not assumed right a and b\n";
        return 0;
    }
    double cc = a;
    while ((b-a) >= EP) {
        // Find middle point
        cc = (a+b)/2;
        // Check if middle point is root
        fc = washout(pow(10.,cc),pp,name) - z_cut;
        if (fc == 0.0)
            break;
        // Decide the side to repeat the steps
        else if (fc*fa < 0)
            b = cc;
        else
            a = cc;
    }
    //cout << "The value of root is : " << cc << endl;
    
    return cc;
}


double freezeout(param_t &pp)
{
    
    double steps = 40.;
    double sph = 28./79.;
    double hh = 4./steps;
    
    double TT;
    int id = 0;
    double zi = 0.;
    double zz_in = pow(10.,bisection(-1., 2.,1.,pp,"L"));//zi;
    double zz_out = pow(10.,bisection(-1., 2.,1.E-03,pp,"L"));//zi;Mscale/100.;//
     
    double tolerance = 0.;//1.E-05;
    std::cout << "Wash-out L :: " << zz_in << std::endl;
    std::cout << "Wash-out L 1.E-3 :: " << zz_out << std::endl;
    int N_proc;
    int out_proc;
    
    size_t nn = Xptl.size();
    
    Eigen::VectorXcd yn(nn),eps1(Xptl.size());
    
    Eigen::MatrixXcd eps2a(Xptl.size(),Xptl.size()),eps2b(Xptl.size(),Xptl.size());
    
    Eigen::VectorXcd tempf(nn),bth(nn),yth(nn),ypth(nn),yeq(nn);
    
    complex<double> a_asym,b_asym ,temp_asym, y_asym,yp_asym, y_temp;
    
    Eigen::MatrixXcd tempfp(nn,nn),athinv(nn,nn),exath(nn,nn);
    Eigen::MatrixXd ath(nn,nn);
    
    yn = Eigen::VectorXcd::Zero(nn);
    yn << 0.,0.,0.,0.;
    
    y_asym = 0.;
    
    double eps = epsilon1("L","N_3",pp);
    double factor = sph/(8.6E-11);
    double YeqN3 = Yeq("N_3",Mscale/zz_out);
    std::cout << "asym = " << factor*eps*YeqN3 << std::endl;
    //for(const auto &Xn : Xptl) std::cout << "X>> " << Xn << std::endl;
    
    for(const auto &Xn : SMptl) std::cout << "SM>> " << Xn << std::endl;
    
    for(size_t it = 0;it!=Xptl.size();++it) std::cout << "X["<<it<<"] ="<< Xptl[it] << std::endl;
    
    hh = (log10(zz_out) - log10(zz_in))/steps;
    
    zi = 0.;
    do{
        double zz = pow(10., log10(zz_in) + hh*zi);
        TT = Mscale/zz;
        
        complex<double> Xsec = 0., pr = 0.;
        
        double Hb = sqrt(4.*pow(M_PI,3.)*geffT(TT)/45.)*pow(TT,2.)/mpl;
        double sT = geffT(TT)*2*pow(M_PI,2.)/45.;
        double yeql = 4./pow(M_PI,2.)*pow(sT,-1.);
        double s = geffT(TT)*2*pow(M_PI,2.)/45.*pow(TT,3.);
        
        tempf = Eigen::VectorXcd::Zero(nn);
        tempfp = Eigen::MatrixXcd::Zero(nn,nn);
        a_asym = 0.;
        temp_asym = 0.;
        
        yth = yn;
        y_temp = y_asym;
        
        size_t si = 0;
        
        auto k1byk2 = [&](const double& zx) { return zx/(1.5 + zx); };
        
        for(const auto &mname : Xptl){
            yeq(si) = Yeq(mname,TT);
            tempf(si) += (yth(si)+1.)*k1byk2(ptlmap[mname]/TT);
            tempfp(si,si) += k1byk2(ptlmap[mname]/TT);
            ++si;
        }
        
        for(size_t ii = 0;ii!=Xptl.size();++ii){
            for(size_t jj = 0;jj!=Xptl.size();++jj){
                
                for(const auto &prname : procQ){
                    N_proc = prname.inParticles.size();
                    if(N_proc>1){
                        if((prname.inParticles[0].name.find(Xptl[ii]) != std::string::npos)&&(prname.inParticles[1].name.find(Xptl[jj] ) != std::string::npos))
                        {
                            Xsec = 0.;
                            findAmp(prname.name,id);
                            pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                            Xsec += (!isnan(pr.real())&&(pr.real()>0.)) ? pr : 0.;
                        }
                        tempf(ii) += -log(10)*zz/yeq(ii)*((yth(ii)*yth(jj) + yth(ii) + yth(jj))*Xsec);
                        
                        for(size_t kk = 0;kk!=Xptl.size();++kk){
                            double del_ki = (kk==ii) ? 1. : 0;
                            double del_kj = (kk==jj) ? 1. : 0;
                            tempfp(ii,kk) += -log(10)*zz/yeq(ii)*(del_ki*yth(jj) + del_ki + del_kj + del_kj*yth(ii))*Xsec;
                        }
                    }
                    
                    for(const auto &pin : prname.inParticles){
                        for(const auto &pout : prname.outParticles){
                            if((pin.name.find(Xptl[ii]) != std::string::npos)&&(pout.name.find(Xptl[jj] ) != std::string::npos)){
                                Xsec = 0.;
                                
                                findAmp(prname.name,id);
                                pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                
                                Xsec += (!isnan(pr.real())&&(pr.real()>0.)) ? pr : 0.;
                            }
                        }
                        tempf(ii) += -log(10)*zz/yeq(ii)*((yth(ii) - yth(jj))*Xsec);
                        
                        for(size_t kk = 0;kk!=Xptl.size();++kk){
                            double del_ki = (kk==ii) ? 1. : 0;
                            double del_kj = (kk==jj) ? 1. : 0;
                            tempfp(ii,kk) += -log(10)*zz/yeq(ii)*((del_ki - del_kj)*Xsec);
                        }
                    }
                    
                }
                for(const auto &prname : procL){
                    complex<double> XsecL = 0.;
                    
                    //====================================
                    double DeltaL = prname.qNumbers.begin()->second;
                    N_proc = prname.inParticles.size();
                    out_proc = prname.outParticles.size();
                    if((N_proc<2)&&(out_proc<3)){
                        for(const auto &outptl : prname.outParticles){
                            if((prname.inParticles[0].name.find(Xptl[ii])!=std::string::npos)&&(outptl.name.find(Xptl[jj])!=std::string::npos)){
                                Xsec = 0.;
                                XsecL = 0.;
                                findAmp(prname.name,id);
                                pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                Xsec += (!isnan(pr.real())&&(pr.real()>0.)) ? 0.5*abs(DeltaL)*pr : 0.;
                                if(DeltaL>0){
                                    pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,1)/(s*Hb*zz) : 0.;
                                    XsecL += (!isnan(pr.real())) ? 2.*pr : 0.;
                                }
                                if(DeltaL<0){
                                    pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,1)/(s*Hb*zz) : 0.;
                                    XsecL -= (!isnan(pr.real())) ? 2.*pr : 0.;
                                }
                                temp_asym += log(10)*zz*((yth(ii) - yth(jj))*XsecL - y_temp/yeql*Xsec*(yth(jj) + 1.) );
                                a_asym += -log(10)*zz/yeql*Xsec*(yth(jj) + 1.);
                                
                            }
                        }
                        
                    }
                    
                    if((N_proc>1)&&(out_proc<3)){
                        for(const auto &inptl : prname.inParticles){
                            for(const auto &outptl : prname.outParticles){
                                if((inptl.name.find(Xptl[ii]) != std::string::npos)&&(outptl.name.find(Xptl[jj]) != std::string::npos)){
                                    Xsec = 0.;
                                    XsecL = 0.;
                                    
                                    findAmp(prname.name,id);
                                    pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                    Xsec += (!isnan(pr.real())&&(pr.real()>0.)) ? 0.5*abs(DeltaL)*pr : 0.;
                                    /*
                                     if(DeltaL>0){
                                     pr = (id < ((int) f_G.size())) ? abs(gamma(id,TT,pp,prname,1)/(s*Hb*zz)) : 0.;
                                     XsecL += (pr.real()>tolerance) ? 2.*pr : 0.;
                                     }
                                     if(DeltaL<0){
                                     pr = (id < ((int) f_G.size())) ? abs(gamma(id,TT,pp,prname,1)/(s*Hb*zz)) : 0.;
                                     XsecL -= (pr.real()>tolerance) ? 2.*pr : 0.;
                                     }
                                     
                                     */
                                    temp_asym += log(10)*zz*( (yth(ii) - yth(jj))*XsecL - y_temp/yeql*Xsec*(yth(ii) + yth(jj) + 2.));
                                    a_asym += -log(10)*zz/yeql*Xsec*(yth(ii) + yth(jj) + 2.);
                                    
                                }
                            }
                        }
                        
                        for(const auto &inptl1 : prname.inParticles){
                            for(const auto &inptl2 : prname.inParticles){
                                if((inptl1.name.find(Xptl[ii]) != std::string::npos)&&(inptl2.name.find(Xptl[jj]) != std::string::npos)){
                                    Xsec = 0.;
                                    XsecL = 0.;
                                    
                                    findAmp(prname.name,id);
                                    pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                    Xsec += (!isnan(pr.real())&&(pr.real()>0.)) ? 0.5*abs(DeltaL)*pr : 0.;
                                    /*
                                     if(DeltaL>0){
                                     pr = (id < ((int) f_G.size())) ? abs(gamma(id,TT,pp,prname,1)/(s*Hb*zz)) : 0.;
                                     XsecL += (pr.real()>tolerance) ? 2.*pr : 0.;
                                     }
                                     if(DeltaL<0){
                                     pr = (id < ((int) f_G.size())) ? abs(gamma(id,TT,pp,prname,1)/(s*Hb*zz)) : 0.;
                                     XsecL -= (pr.real()>tolerance) ? 2.*pr : 0.;
                                     }
                                     */
                                    temp_asym += log(10)*zz*( (yth(ii)*yth(jj) + yth(ii) + yth(jj) + 1.)*XsecL - y_temp/yeql*Xsec);
                                    a_asym += -log(10)*zz/yeql*Xsec;
                                    
                                }
                            }
                        }
                        
                        
                    }
                }
                
                    //======================================
                
            }
        }

        ath = tempfp.real();
        bth = tempf - tempfp*yth;
        
        b_asym = temp_asym - a_asym*y_temp;
        
        int dim = Xptl.size();
        
        if((ath.determinant()!=0)&&(!isnan(ath.determinant()))&&(!isinf(ath.determinant())))
        {
            ypth = expM(ath*hh, dim)*yth + tempfp.inverse()*(expM(ath*hh, dim) -  MatrixXd::Identity(dim, dim))*bth;
        }
        else
        {
            ypth = yth + MatrixXd::Identity(dim, dim)*bth*hh;
        }
        
        if(a_asym!=0.){
            yp_asym = exp(a_asym*hh)*y_temp + 1./a_asym*(exp(a_asym*hh) - 1.)*b_asym;
        }
        if(a_asym==0.){
            yp_asym = y_temp + b_asym*hh;
        }
        
        yn = ypth;
        y_asym = yp_asym;
        //double washout_Q = washout(zz,pp,"Q");
        double washout_L = washout(zz,pp,"L");
        
        std::cout << TT << "\t"<< zz << "\t" << washout_L << "\t" << sph*y_asym.real()/(8.6E-11) << std::endl;
        ++zi;
    }while(zi<=steps);//while((TT>=100.));//
    
    return sph*y_asym.real()/(8.6E-11);
    
}

void BESolver(std::ofstream &ofile,param_t &pp)
{
    int N_proc;
    int out_proc;
    int id = 0;
    int DeltaL1 = 1,DeltaL2 = 2;
    
    size_t nn = Xptl.size();
    
    Eigen::VectorXcd yn(nn),eps1(Xptl.size());
    
    Eigen::MatrixXcd eps2a(Xptl.size(),Xptl.size()),eps2b(Xptl.size(),Xptl.size());
    
    Eigen::VectorXcd tempf(nn),bth(nn),yth(nn),ypth(nn),yeq(nn);
    
    complex<double> a_asym,b_asym ,temp_asym, y_asym, yp_asym, y_temp;
    
    Eigen::MatrixXcd tempfp(nn,nn),athinv(nn,nn),exath(nn,nn);
    Eigen::MatrixXd ath(nn,nn);
    
    yn = Eigen::VectorXcd::Zero(nn);
    
    yn << 0.,0.,0.,0.;
    
    y_asym = 0.;
    
    eps1 << epsilon1("L","N_3",pp),epsilon1("L","N_2",pp),epsilon1("L","N_1",pp),0.;
    
    std::cout << "eps =  " << eps1(0) << " " << eps1(1) << " " << eps1(2) << std::endl;
    eps2a << 0., 0.,0.,0.,
             0., 0.,0.,0.,
             0., 0.,0.,0.,
             0., 0.,0.,0.;
    
    eps2b << 0., 0.,0.,0.,
             0., 0.,0.,0.,
             0., 0.,0.,0.,
             0., 0.,0.,0.;

    double steps = 1000.;
    double sph = 28./79.;
    
    double zi = 0.1;
    double TT;
    double tolerance = 0.;//1.E-02;
    double zz_in = zi;//pow(10.,bisection(-1., 2.,1.,pp,"Q"));//
    double hh = (2. - log10(zz_in))/steps;
    
    std::cout << "Wash-out " << pow(10.,bisection(-1., 2.,1.,pp,"Q")) << std::endl;
    //double hh = 2./steps;
    do{
        double zz = pow(10., log10(zz_in) + hh*zi);
        //double zz = pow(10., hh*zi);
        TT = Mscale/zz;
        
        complex<double> Xsec = 0., pr = 0.;
        
        double Hb = sqrt(4.*pow(M_PI,3.)*geffT(TT)/45.)*pow(TT,2.)/mpl;
        
        double s = geffT(TT)*2*pow(M_PI,2.)/45.*pow(TT,3.);
        
        double sT = geffT(TT)*2*pow(M_PI,2.)/45.;
        
        double yeql = 4./pow(M_PI,2.)*pow(sT,-1.);
        /**/
        tempf = Eigen::VectorXcd::Zero(nn);
        tempfp = Eigen::MatrixXcd::Zero(nn,nn);
        a_asym = 0.;
        temp_asym = 0.;
        
        yth = yn;
        
        y_temp = y_asym;
        
        
        size_t si = 0;
        
        auto k1byk2 = [&](const double& zx) { return zx/(1.5 + zx); };
        
        for(const auto &mname : Xptl){
            yeq(si) = Yeq(mname,TT);
            tempf(si) += (yth(si)+1.)*k1byk2(ptlmap[mname]/TT);
            tempfp(si,si) += k1byk2(ptlmap[mname]/TT);
            ++si;
        }
        
        std::cout << std::endl;
        
        for(size_t ii = 0;ii!=Xptl.size();++ii){
            for(size_t jj = 0;jj!=Xptl.size();++jj){
                
                for(const auto &prname : procQ){
                    N_proc = prname.inParticles.size();
                    if(N_proc>1){
                        if((prname.inParticles[0].name.find(Xptl[ii]) != std::string::npos)&&(prname.inParticles[1].name.find(Xptl[jj] ) != std::string::npos))
                        {
                            Xsec = 0.;
                            findAmp(prname.name,id);
                            pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                            //Xsec += (pr.real()>0.1*k1byk2(ptlmap[Xptl[ii]]/TT)*yeq(ii).real()) ? pr : 0.;
                            Xsec += (!isnan(pr.real())&&(pr.real()>0.)) ? pr : 0.;
                        }
                        tempf(ii) += -log(10)*zz/yeq(ii)*((yth(ii)*yth(jj) + yth(ii) + yth(jj))*Xsec);
                        
                        for(size_t kk = 0;kk!=Xptl.size();++kk){
                            double del_ki = (kk==ii) ? 1. : 0;
                            double del_kj = (kk==jj) ? 1. : 0;
                            tempfp(ii,kk) += -log(10)*zz/yeq(ii)*(del_ki*yth(jj) + del_ki + del_kj + del_kj*yth(ii))*Xsec;
                        }
                    }
                    
                    for(const auto &pin : prname.inParticles){
                        for(const auto &pout : prname.outParticles){
                            if((pin.name.find(Xptl[ii]) != std::string::npos)&&(pout.name.find(Xptl[jj] ) != std::string::npos)){
                                Xsec = 0.;
                                
                                findAmp(prname.name,id);
                                pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                //Xsec += (pr.real()>0.1*k1byk2(ptlmap[Xptl[ii]]/TT)*yeq(ii).real()) ? pr : 0.;
                                Xsec += (!isnan(pr.real())&&(pr.real()>0.)) ? pr : 0.;
                            }
                        }
                        tempf(ii) += -log(10)*zz/yeq(ii)*((yth(ii) - yth(jj))*Xsec);
                        
                        for(size_t kk = 0;kk!=Xptl.size();++kk){
                            double del_ki = (kk==ii) ? 1. : 0;
                            double del_kj = (kk==jj) ? 1. : 0;
                            tempfp(ii,kk) += -log(10)*zz/yeq(ii)*((del_ki - del_kj)*Xsec);
                        }
                    }
                    
                }
                //============= Asymmetry ===================================
                for(const auto &prname : procL){
                    complex<double> XsecL = 0.;
                    
                    //double DeltaL = 1.;
                    //====================================
                    double DeltaL = prname.qNumbers.begin()->second;
                    N_proc = prname.inParticles.size();
                    out_proc = prname.outParticles.size();
                    
                    if((N_proc<2)&&(out_proc<3)){
                        for(const auto &outptl : prname.outParticles){
                            if((prname.inParticles[0].name.find(Xptl[ii])!=std::string::npos)&&(outptl.name.find(Xptl[jj])!=std::string::npos)){
                                Xsec = 0.;
                                XsecL = 0.;
                                findAmp(prname.name,id);
                                pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                Xsec += (!isnan(pr.real())&&(pr.real()>0.)) ? 0.5*abs(DeltaL)*pr : 0.;
                                if(DeltaL>0){
                                    pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,1)/(s*Hb*zz) : 0.;
                                    XsecL += (!isnan(pr.real())) ? 2.*pr : 0.;
                                }
                                if(DeltaL<0){
                                    pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,1)/(s*Hb*zz) : 0.;
                                    XsecL -= (!isnan(pr.real())) ? 2.*pr : 0.;
                                }
                                temp_asym += log(10)*zz*((yth(ii) - yth(jj))*XsecL - y_temp/yeql*Xsec*(yth(jj) + 1.) );
                                a_asym += -log(10)*zz/yeql*Xsec*(yth(jj) + 1.);
                                
                            }
                        }
                        
                    }
                    
                    /**/
                    if((N_proc>1)&&(out_proc<3)){
                        
                        complex<double> ytemp = 0.,yX = 0.;
                        for(const auto &inptl : prname.inParticles){
                            for(const auto &outptl : prname.outParticles){
                                if((inptl.name.find(Xptl[ii]) != std::string::npos)&&(outptl.name.find(Xptl[jj]) != std::string::npos)){
                                    Xsec = 0.;
                                    XsecL = 0.;
                                    
                                    findAmp(prname.name,id);
                                    pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                    Xsec += (!isnan(pr.real())&&(pr.real()>0.)) ? 0.5*abs(DeltaL)*pr : 0.;
                                    /*
                                     if(DeltaL>0){
                                     pr = (id < ((int) f_G.size())) ? abs(gamma(id,TT,pp,prname,1)/(s*Hb*zz)) : 0.;
                                     XsecL += (pr.real()>tolerance) ? 2.*pr : 0.;
                                     }
                                     if(DeltaL<0){
                                     pr = (id < ((int) f_G.size())) ? abs(gamma(id,TT,pp,prname,1)/(s*Hb*zz)) : 0.;
                                     XsecL -= (pr.real()>tolerance) ? 2.*pr : 0.;
                                     }*/
                                    
                                    temp_asym += log(10)*zz*( (yth(ii) - yth(jj))*XsecL - y_temp/yeql*Xsec*(yth(ii) + yth(jj) + 2.));
                                    a_asym += -log(10)*zz/yeql*Xsec*(yth(ii) + yth(jj) + 2.);
                                    
                                }
                            }
                        }
                        
                        for(const auto &inptl1 : prname.inParticles){
                            for(const auto &inptl2 : prname.inParticles){
                                if((inptl1.name.find(Xptl[ii]) != std::string::npos)&&(inptl2.name.find(Xptl[jj]) != std::string::npos)){
                                    Xsec = 0.;
                                    XsecL = 0.;
                                    
                                    findAmp(prname.name,id);
                                    pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                    Xsec += (!isnan(pr.real())&&(pr.real()>0.)) ? 0.5*abs(DeltaL)*pr : 0.;
                                    /*
                                     if(DeltaL>0){
                                     pr = (id < ((int) f_G.size())) ? abs(gamma(id,TT,pp,prname,1)/(s*Hb*zz)) : 0.;
                                     XsecL += (pr.real()>tolerance) ? 2.*pr : 0.;
                                     }
                                     if(DeltaL<0){
                                     pr = (id < ((int) f_G.size())) ? abs(gamma(id,TT,pp,prname,1)/(s*Hb*zz)) : 0.;
                                     XsecL -= (pr.real()>tolerance) ? 2.*pr : 0.;
                                     }*/
                                    
                                    temp_asym += log(10)*zz*( (yth(ii)*yth(jj) + yth(ii) + yth(jj) + 1.)*XsecL - y_temp/yeql*Xsec);
                                    a_asym += -log(10)*zz/yeql*Xsec;
                                    
                                }
                            }
                        }
                    }
                }
                
                //======================================
                
            }
        }
        
        ath = tempfp.real();
        bth = tempf - tempfp*yth;
        
        b_asym = temp_asym - a_asym*y_temp;
        
        int dim = Xptl.size();
        
        if((ath.determinant()!=0)&&(!isnan(ath.determinant()))&&(!isinf(ath.determinant())))
        {
            ypth = expM(ath*hh, dim)*yth + tempfp.inverse()*(expM(ath*hh, dim) -  MatrixXd::Identity(dim, dim))*bth;
        }
        else
        {
            ypth = yth + MatrixXd::Identity(dim, dim)*bth*hh;
        }
        
        if(a_asym!=0.){
            yp_asym = exp(a_asym*hh)*y_temp + 1./a_asym*(exp(a_asym*hh) - 1.)*b_asym;
        }
        if(a_asym==0.){
            yp_asym = y_temp + b_asym*hh;
        }
        
        yn = ypth;
        y_asym = yp_asym;
        
        double washout_L = washout(zz,pp,"L");
        
        ofile << TT << "\t"<< zz << "\t" << washout_L << "\t" << y_asym.real() << std::endl;
        std::cout << TT << "\t"<< zz << "\t" << washout_L << "\t" << sph*y_asym.real()/(8.6E-11) << std::endl;
        
        std::cout << yn.transpose() << std::endl;
      
        ++zi;
    }while((zi<=steps));
    
}


}
