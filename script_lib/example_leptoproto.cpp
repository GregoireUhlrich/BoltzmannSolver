
#include "leptoproto.h"
#include "BE.h"

using namespace std;
using namespace csl;
using namespace mty;


using namespace mty::lib;


using namespace leptoproto;

/**/
void analysis(int i_start,int i_end)
{
    double ii = ((double) i_start);
    
    std::string filename = "output-" + std::to_string(i_start) + std::to_string(i_end) + ".txt";
    std::ofstream output_file(filename);
    
    do
    {
        //double MN3 = ((double) ii)*1000.;
        //for(auto jj=1;jj!=6;++jj)
        {
            //int xx = rand()%20 + 1;
            //int yy = rand()%10 + 1;
            
            //double MN3 = ((double) xx)*1000.;
            double MN3 = ((double) ii)*100. + 1000.;
            
            double diff = 100.;//((double) jj)*100.;
            
            for(auto kk=1;kk!=20;++kk)
            {
                Eigen::VectorXd MN(3);
                
                MN << 100.,MN3 - 2.*diff,MN3;
                double meta = MN3 - diff;
                double muSet = ((double) kk)*0.5*MN3 + MN3;
                std::cout << MN3 << "\t" << muSet/MN3 << "\t" << diff << "\t" << std::endl;
                //double muSet = ((double) yy)*MN3;
                
                //std::cout << "K1 = " << cyl_bessel_k(1, 100.) << std::endl;
                double alph_e = 1./137.;
                double vev = 246.;
                double lam5 = 1.e+00;
                double lamet = 0.5;
                
                auto F = [&](double const& x){ return x/(x-1.)*log(x);};
                
                double th12 = M_PI*33.8/180.;
                double th13 = M_PI*8.61/180.;
                double th23 = M_PI*46.2/180.;
                double del = 281.2/180.*M_PI;
                
                Eigen::MatrixXcd mnu(3,3);
                mnu = Eigen::MatrixXd::Zero(3,3);
                
                double mnu1 = 1.0E-18*1.0E-09;
                double Dm21 = 7.5E-05*1.0E-18;
                double meff = 2.4*1.0E-3*1.0E-18;
                double mnu2 = std::sqrt(Dm21 + mnu1*mnu1);
                double mnu3N = std::sqrt(meff + Dm21 + mnu1*mnu1);
                //ouble mnu3I = std::sqrt(meff + mnu1*mnu1);
                
                // NH
                mnu << sqrt(mnu1), 0., 0.,
                0., sqrt(mnu2), 0.,
                0., 0., sqrt(mnu3N);
                //  IH
                // mnu << std::sqrt(mnu3I), 0., 0.,
                // 0., std::sqrt(mnu1), 0.,
                // 0., 0., std::sqrt(mnu2);
                
                //std::cout << "mnu >> " << std::endl << mnu<< std::endl;
                Eigen::MatrixXcd O12(3,3);
                Eigen::MatrixXcd O13(3,3);
                Eigen::MatrixXcd O23(3,3);
                Eigen::MatrixXcd MM(3,3);
                
                Eigen::MatrixXcd UPMNS(3,3);
                
                O12 << cos(th12),-sin(th12),0.,
                sin(th12),cos(th12),0.,
                0.,0.,1.;
                O13 << cos(th13),0.,-sin(th13)*std::exp(-del*1i),
                0.,1.,0.,
                sin(th13)*exp(del*1i),0.,cos(th13);
                O23 << 1.,0.,0.,
                0.,cos(th23),-sin(th23),
                0.,sin(th23),cos(th23);
                
                Eigen::MatrixXcd sqLam(3,3);
                
                double MetaR = pow(meta,2.) + 2.*lam5*vev*vev;
                double MetaI = pow(meta,2.) - 2.*lam5*vev*vev;
                
                complex<double> lam1,lam2,lam3;
                
                lam1 = MN(0)/(32.*pow(M_PI,2.))*(F(MetaR/pow(MN(0),2.))
                                                 - F(MetaI/pow(MN(0),2.)));
                
                lam2 = MN(1)/(32.*pow(M_PI,2.))*(F(MetaR/pow(MN(2),2.))
                                                 - F(MetaI/pow(MN(2),2.)));
                
                lam3 = MN(2)/(32.*pow(M_PI,2.))*(F(MetaR/pow(MN(0),2.))
                                                 - F(MetaI/pow(MN(0),2.)));
                
                sqLam << sqrt(lam1), 0., 0.,
                0., sqrt(lam2), 0.,
                0., 0., sqrt(lam3);
                
                Eigen::MatrixXcd RO23(3,3);
                complex<double> Rth23 = M_PI/8. - M_PI/8.i;
                //complex<double> Rth23 = -M_PI/4.i;
                RO23 << 1.,0.,0.,
                0.,cos(Rth23),-sin(Rth23),
                0.,sin(Rth23),cos(Rth23);
                
                UPMNS = O23 * O13 * O12;
                Eigen::MatrixXcd Yuk;
                Yuk = (sqLam.inverse()*RO23* mnu * (UPMNS.conjugate()).transpose());
                
                param_t pp;
                
                pp.e = std::sqrt(4.*M_PI*alph_e);
                pp.m_N_1 = MN(0);
                pp.m_N_2 = MN(1);
                pp.m_N_3 = MN(2);
                pp.muet2 = -std::pow(meta,2.);
                
                pp.thetaW = std::asin(std::sqrt(0.23116));
                pp.reg_prop = 1.;
                pp.muS2 = -1.;
                pp.lSet = 0.1;
                pp.lamet = lamet;
                pp.Finite = 1.;
                pp.muSet = muSet;//12000.;
                pp.Y_et_00 = Yuk(0,0);
                pp.Y_et_01 = Yuk(0,1);
                pp.Y_et_10 = Yuk(1,0);//0.001*factor;
                pp.Y_et_02 = Yuk(0,2);
                pp.Y_et_20 = Yuk(2,0);//factor;
                pp.Y_et_11 = Yuk(1,1);//0.01i*factor;
                pp.Y_et_12 = Yuk(1,2);//0.01i*factor;
                pp.Y_et_21 = Yuk(2,1);//factor;
                pp.Y_et_22 = Yuk(2,2);//factor;
                
                pp.Y_S_00 = 0.;
                pp.Y_S_01 = 0.0;
                pp.Y_S_10 = 0.0;
                pp.Y_S_02 = 0.0;
                pp.Y_S_20 = 0.0;
                pp.Y_S_11 = 0.;
                pp.Y_S_12 = 0.5;//*std::exp(pi/4.i);
                pp.Y_S_21 = 0.5;//*std::exp(pi/4.i);
                pp.Y_S_22 = 0.0;
                
                updateSpectrum(pp);
                
                //pp.print();
                
                //double yB_ratio = freezeout(pp);
                double yB_ratio = freezeout(pp);
                
                output_file << MN3 << "\t" << muSet/MN3 << "\t" << diff << "\t" << yB_ratio << std::endl;
                
            }
        }
        
        ++ii;
        
    }while(ii<=i_end);
    output_file.close();
}

int main(int argc, char* argv[]) {
    
    Eigen::VectorXd MN(3);
    
    MN << 100.,atof(argv[3]),atof(argv[1]);
    double meta = atof(argv[2]);
    double muSet = atof(argv[4]);
    
    std::cout << "K1 = " << cyl_bessel_k(1, 100.) << std::endl;
    double alph_e = 1./137.;
    double vev = 246.;
    double lam5 = 1.e+00;
    double lamet = 0.5;
    
    auto F = [&](double const& x){ return x/(x-1.)*log(x);};
    
    double th12 = M_PI*33.8/180.;
    double th13 = M_PI*8.61/180.;
    double th23 = M_PI*46.2/180.;
    double del = 281.2/180.*M_PI;
    
    Eigen::MatrixXcd mnu(3,3);
    mnu = Eigen::MatrixXd::Zero(3,3);
    
    double mnu1 = 1.0E-18*1.0E-09;
    double Dm21 = 7.5E-05*1.0E-18;
    double meff = 2.4*1.0E-3*1.0E-18;
    double mnu2 = std::sqrt(Dm21 + mnu1*mnu1);
    double mnu3N = std::sqrt(meff + Dm21 + mnu1*mnu1);
    //double mnu3I = std::sqrt(meff + mnu1*mnu1);
    
    /* NH*/
    mnu << sqrt(mnu1), 0., 0.,
              0., sqrt(mnu2), 0.,
          0., 0., sqrt(mnu3N);
    /*  IH
    mnu << std::sqrt(mnu3I), 0., 0.,
              0., std::sqrt(mnu1), 0.,
              0., 0., std::sqrt(mnu2);*/
    
    std::cout << "mnu >> " << std::endl << mnu<< std::endl;
    Eigen::MatrixXcd O12(3,3);
    Eigen::MatrixXcd O13(3,3);
    Eigen::MatrixXcd O23(3,3);
    Eigen::MatrixXcd MM(3,3);
    
    Eigen::MatrixXcd UPMNS(3,3);
    
    O12 << cos(th12),-sin(th12),0.,
    sin(th12),cos(th12),0.,
    0.,0.,1.;
    O13 << cos(th13),0.,-sin(th13)*std::exp(-del*1i),
    0.,1.,0.,
    sin(th13)*exp(del*1i),0.,cos(th13);
    O23 << 1.,0.,0.,
    0.,cos(th23),-sin(th23),
    0.,sin(th23),cos(th23);
    
    Eigen::MatrixXcd sqLam(3,3);
    
    double MetaR = pow(meta,2.) + 2.*lam5*vev*vev;
    double MetaI = pow(meta,2.) - 2.*lam5*vev*vev;
    
    complex<double> lam1,lam2,lam3;
    
    lam1 = MN(0)/(32.*pow(M_PI,2.))*(F(MetaR/pow(MN(0),2.))
     - F(MetaI/pow(MN(0),2.)));
    
    lam2 = MN(1)/(32.*pow(M_PI,2.))*(F(MetaR/pow(MN(2),2.))
     - F(MetaI/pow(MN(2),2.)));
    
    lam3 = MN(2)/(32.*pow(M_PI,2.))*(F(MetaR/pow(MN(0),2.))
     - F(MetaI/pow(MN(0),2.)));
    
    sqLam << sqrt(lam1), 0., 0.,
     0., sqrt(lam2), 0.,
     0., 0., sqrt(lam3);
    
    Eigen::MatrixXcd RO23(3,3);
    complex<double> Rth23 = M_PI/8. - M_PI/8.i;
    //complex<double> Rth23 = -M_PI/4.i;
    RO23 << 1.,0.,0.,
    0.,cos(Rth23),-sin(Rth23),
    0.,sin(Rth23),cos(Rth23);
    
    UPMNS = O23 * O13 * O12;
    Eigen::MatrixXcd Yuk;
    Yuk = (sqLam.inverse()*RO23* mnu * (UPMNS.conjugate()).transpose());
    
    std::cout << "sqLam >> " << std::endl << sqLam << std::endl;
    std::cout << "UPMNS >> " << std::endl << UPMNS << std::endl;
    
    std::cout << "sqLam >> " << std::endl << sqLam.inverse() << std::endl;
    std::cout << "Yukawa >> " << std::endl << Yuk << std::endl;
    
    param_t pp;
    
    pp.e = std::sqrt(4.*M_PI*alph_e);
    pp.m_N_1 = MN(0);
    pp.m_N_2 = MN(1);
    pp.m_N_3 = MN(2);
    pp.muet2 = -std::pow(meta,2.);
    
    pp.thetaW = std::asin(std::sqrt(0.23116));
    pp.reg_prop = 1.;
    pp.muS2 = -1.;
    pp.lSet = 0.1;
    pp.lamet = lamet;
    pp.Finite = 1.;
    pp.muSet = muSet;//12000.;
    pp.Y_et_00 = Yuk(0,0);
    pp.Y_et_01 = Yuk(0,1);
    pp.Y_et_10 = Yuk(1,0);//0.001*factor;
    pp.Y_et_02 = Yuk(0,2);
    pp.Y_et_20 = Yuk(2,0);//factor;
    pp.Y_et_11 = Yuk(1,1);//0.01i*factor;
    pp.Y_et_12 = Yuk(1,2);//0.01i*factor;
    pp.Y_et_21 = Yuk(2,1);//factor;
    pp.Y_et_22 = Yuk(2,2);//factor;
    
    pp.Y_S_00 = 0.;
    pp.Y_S_01 = 0.0;
    pp.Y_S_10 = 0.0;
    pp.Y_S_02 = 0.0;
    pp.Y_S_20 = 0.0;
    pp.Y_S_11 = 0.;
    pp.Y_S_12 = 0.5;//*std::exp(pi/4.i);
    pp.Y_S_21 = 0.5;//*std::exp(pi/4.i);
    pp.Y_S_22 = 0.0;
    
    updateSpectrum(pp);
    
    pp.print();
    
    initialize(pp);
    
    std::cout << "eps = " << epsilon1("L","N_3",pp) << std::endl;
    //std::cout << "eps = " << epsilon1("L","N_3c",pp) << std::endl;
    std::vector<Process> procL,procQ;
    
    ParticleData pdata;
    
    pdata.loadFile("script/test.json");
    
    //data.load("test.json");
    std::cout << std::endl << " ======= L =========" << std::endl;
    std::vector<std::string> particleNames = pdata.getParticleNames();
    for (const auto &name : particleNames) {
        std::cout << name << " : " << pdata.getQuantumNumberValue("L", name).getValue() << std::endl;
    }
    std::cout << std::endl << " ======= Q =========" << std::endl;
    std::vector<std::string> particleName = pdata.getParticleNames();
    for (const auto &name : particleName) {
        std::cout << name << " : " << pdata.getQuantumNumberValue("Q", name).getValue() << std::endl;
    }
    /*
    std::cout << std::endl << " ======= Th =========" << std::endl;
    std::vector<std::string> ptlth = pdata.getParticleNames();
    for (const auto &name : ptlth) {
        std::cout << name << " : " << pdata.getQuantumNumberValue("Th", name).getValue() << std::endl;
        std::cout << "Mass  >> " << Massp(name,pp) << std::endl;
    }
     */
    
    procL = pdata.getProcessesFromQNumberViolation("L");
    procQ = pdata.getProcessesFromQNumberConservation("Q");
    
    
    
    auto startTime = std::chrono::high_resolution_clock::now();
    double yB_ratio = freezeout(pp);
    //double yB_ratio = freezeout(pp);
    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    cout << "Sec: " << duration.count()/1000000 << endl;
    std::cout << "yB =  " << yB_ratio << std::endl;
    
    
    /**/
    std::string pname = argv[5];
    std::ofstream ofile(pname);
    
    BESolver(ofile,pp);
    ofile.close();
     
     
     
    //analysis(0,10);
    /*
    thread t1(analysis,1,50,ap);
    thread t2(analysis,51,100,ap);
    thread t3(analysis,101,150,ap);
    thread t4(analysis,151,200,ap);
    thread t5(analysis,21,25);
    thread t6(analysis,50,59);
    thread t7(analysis,60,69);
    thread t8(analysis,70,79);
    thread t9(analysis,80,89);
    thread t10(analysis,90,99);
     */
     
    /*
    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    t6.join();
    t7.join();
    t8.join();
    t9.join();
    t10.join();
               */
    
    //std::cout << "Process for Q = 0" << std::endl;
    //proc = pdata.getProcessesFromQNumberConservation("Q");
    
    //for(const auto &prname : proc)
    //std::cout << "name = " << prname.name << std::endl;
    
    return 0;
}
