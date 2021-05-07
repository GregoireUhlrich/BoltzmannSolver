#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_xx_df2_to_uf3c_df3c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XX_df2_to_uf3c_df3c(
        param_t const &param
        )
{
    clearcache();
    auto const &gw = param.gw;
    auto const &mX = param.mX;
    auto const &m_b = param.m_b;
    auto const &m_s = param.m_s;
    auto const &m_t = param.m_t;
    auto const &GbRt = param.GbRt;
    auto const &GcRt = param.GcRt;
    auto const &GdRt = param.GdRt;
    auto const &GsRt = param.GsRt;
    auto const &GtRt = param.GtRt;
    auto const &GuRt = param.GuRt;
    auto const &gwdR = param.gwdR;
    auto const &gwuR = param.gwuR;
    auto const &s_23 = param.s_23;
    auto const &s_24 = param.s_24;
    auto const &s_25 = param.s_25;
    auto const &s_34 = param.s_34;
    auto const &s_35 = param.s_35;
    auto const &s_45 = param.s_45;
    auto const &m_dtR1 = param.m_dtR1;
    auto const &m_dtR2 = param.m_dtR2;
    auto const &m_dtR3 = param.m_dtR3;
    auto const &m_utR1 = param.m_utR1;
    auto const &m_utR2 = param.m_utR2;
    auto const &m_utR3 = param.m_utR3;
    auto const &lpp_012 = param.lpp_012;
    auto const &lpp_021 = param.lpp_021;
    auto const &lpp_112 = param.lpp_112;
    auto const &lpp_121 = param.lpp_121;
    auto const &lpp_201 = param.lpp_201;
    auto const &lpp_202 = param.lpp_202;
    auto const &lpp_210 = param.lpp_210;
    auto const &lpp_212 = param.lpp_212;
    auto const &lpp_220 = param.lpp_220;
    auto const &lpp_221 = param.lpp_221;
    auto const &U_dtR_00 = param.U_dtR_00;
    auto const &U_dtR_01 = param.U_dtR_01;
    auto const &U_dtR_02 = param.U_dtR_02;
    auto const &U_dtR_10 = param.U_dtR_10;
    auto const &U_dtR_11 = param.U_dtR_11;
    auto const &U_dtR_12 = param.U_dtR_12;
    auto const &U_dtR_20 = param.U_dtR_20;
    auto const &U_dtR_21 = param.U_dtR_21;
    auto const &U_dtR_22 = param.U_dtR_22;
    auto const &U_utR_00 = param.U_utR_00;
    auto const &U_utR_01 = param.U_utR_01;
    auto const &U_utR_02 = param.U_utR_02;
    auto const &U_utR_10 = param.U_utR_10;
    auto const &U_utR_11 = param.U_utR_11;
    auto const &U_utR_12 = param.U_utR_12;
    auto const &U_utR_20 = param.U_utR_20;
    auto const &U_utR_21 = param.U_utR_21;
    auto const &U_utR_22 = param.U_utR_22;
    const complex_t IT_0041 = lpp_201*U_dtR_00;
    const complex_t IT_0042 = lpp_221*U_dtR_20;
    const complex_t IT_0043 = IT_0041 + IT_0042;
    const complex_t IT_0044 = lpp_210*U_dtR_00;
    const complex_t IT_0045 = lpp_212*U_dtR_20;
    const complex_t IT_0046 = -IT_0044 + -IT_0045;
    const complex_t IT_0047 = IT_0043 + IT_0046;
    const complex_t IT_0048 = gw*gwdR*std::conj(U_dtR_20);
    const complex_t IT_0049 = IT_0047*IT_0048;
    const complex_t IT_0010 = 2*mX;
    const complex_t IT_0011 = std::pow(IT_0010, 2);
    const complex_t IT_0012 = std::pow(m_dtR1, 2);
    const complex_t IT_0050 = std::pow(m_b, 2);
    const complex_t IT_0051 = std::pow((-2)*s_25 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0011 + -IT_0012 + IT_0050, -1);
    const complex_t IT_0052 = (complex_t{0, 1})*IT_0049*IT_0051;
    const complex_t IT_0053 = lpp_201*U_dtR_01;
    const complex_t IT_0054 = lpp_221*U_dtR_21;
    const complex_t IT_0055 = IT_0053 + IT_0054;
    const complex_t IT_0056 = lpp_210*U_dtR_01;
    const complex_t IT_0057 = lpp_212*U_dtR_21;
    const complex_t IT_0058 = -IT_0056 + -IT_0057;
    const complex_t IT_0059 = IT_0055 + IT_0058;
    const complex_t IT_0060 = gw*gwdR*std::conj(U_dtR_21);
    const complex_t IT_0061 = IT_0059*IT_0060;
    const complex_t IT_0024 = std::pow(m_dtR2, 2);
    const complex_t IT_0062 = std::pow((-2)*s_25 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0011 + -IT_0024 + IT_0050, -1);
    const complex_t IT_0063 = (complex_t{0, 1})*IT_0061*IT_0062;
    const complex_t IT_0064 = lpp_201*U_dtR_02;
    const complex_t IT_0065 = lpp_221*U_dtR_22;
    const complex_t IT_0066 = IT_0064 + IT_0065;
    const complex_t IT_0067 = lpp_210*U_dtR_02;
    const complex_t IT_0068 = lpp_212*U_dtR_22;
    const complex_t IT_0069 = -IT_0067 + -IT_0068;
    const complex_t IT_0070 = IT_0066 + IT_0069;
    const complex_t IT_0071 = gw*gwdR*std::conj(U_dtR_22);
    const complex_t IT_0072 = IT_0070*IT_0071;
    const complex_t IT_0036 = std::pow(m_dtR3, 2);
    const complex_t IT_0073 = std::pow((-2)*s_25 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0011 + -IT_0036 + IT_0050, -1);
    const complex_t IT_0074 = (complex_t{0, 1})*IT_0072*IT_0073;
    const complex_t IT_0075 = IT_0052 + IT_0063 + IT_0074;
    const complex_t IT_0076 = s_25*s_34;
    const complex_t IT_0000 = lpp_202*U_dtR_00;
    const complex_t IT_0001 = lpp_212*U_dtR_10;
    const complex_t IT_0002 = IT_0000 + IT_0001;
    const complex_t IT_0003 = lpp_220*U_dtR_00;
    const complex_t IT_0004 = lpp_221*U_dtR_10;
    const complex_t IT_0005 = -IT_0003 + -IT_0004;
    const complex_t IT_0006 = IT_0002 + IT_0005;
    const complex_t IT_0007 = gw*gwdR*std::conj(U_dtR_10);
    const complex_t IT_0008 = IT_0006*IT_0007;
    const complex_t IT_0009 = std::pow(m_s, 2);
    const complex_t IT_0013 = std::pow(2*s_23 + (complex_t{0, 1})*GdRt*m_dtR1 
      + IT_0009 + IT_0011 + -IT_0012, -1);
    const complex_t IT_0014 = (complex_t{0, 1})*IT_0008*IT_0013;
    const complex_t IT_0015 = lpp_202*U_dtR_01;
    const complex_t IT_0016 = lpp_212*U_dtR_11;
    const complex_t IT_0017 = IT_0015 + IT_0016;
    const complex_t IT_0018 = lpp_220*U_dtR_01;
    const complex_t IT_0019 = lpp_221*U_dtR_11;
    const complex_t IT_0020 = -IT_0018 + -IT_0019;
    const complex_t IT_0021 = IT_0017 + IT_0020;
    const complex_t IT_0022 = gw*gwdR*std::conj(U_dtR_11);
    const complex_t IT_0023 = IT_0021*IT_0022;
    const complex_t IT_0025 = std::pow(2*s_23 + (complex_t{0, 1})*GsRt*m_dtR2 
      + IT_0009 + IT_0011 + -IT_0024, -1);
    const complex_t IT_0026 = (complex_t{0, 1})*IT_0023*IT_0025;
    const complex_t IT_0027 = lpp_202*U_dtR_02;
    const complex_t IT_0028 = lpp_212*U_dtR_12;
    const complex_t IT_0029 = IT_0027 + IT_0028;
    const complex_t IT_0030 = lpp_220*U_dtR_02;
    const complex_t IT_0031 = lpp_221*U_dtR_12;
    const complex_t IT_0032 = -IT_0030 + -IT_0031;
    const complex_t IT_0033 = IT_0029 + IT_0032;
    const complex_t IT_0034 = gw*gwdR*std::conj(U_dtR_12);
    const complex_t IT_0035 = IT_0033*IT_0034;
    const complex_t IT_0037 = std::pow(2*s_23 + (complex_t{0, 1})*GbRt*m_dtR3 
      + IT_0009 + IT_0011 + -IT_0036, -1);
    const complex_t IT_0038 = (complex_t{0, 1})*IT_0035*IT_0037;
    const complex_t IT_0039 = -IT_0014 + -IT_0026 + -IT_0038;
    const complex_t IT_0040 = s_23*s_45;
    const complex_t IT_0077 = s_24*s_35;
    const complex_t IT_0078 = -IT_0077;
    const complex_t IT_0079 = IT_0076 + IT_0078;
    const complex_t IT_0080 = IT_0040 + IT_0079;
    const complex_t IT_0081 = gw*gwuR*std::conj(U_utR_20);
    const complex_t IT_0082 = lpp_021*U_utR_00;
    const complex_t IT_0083 = lpp_121*U_utR_10;
    const complex_t IT_0084 = lpp_221*U_utR_20;
    const complex_t IT_0085 = IT_0082 + IT_0083 + IT_0084;
    const complex_t IT_0086 = lpp_012*U_utR_00;
    const complex_t IT_0087 = lpp_112*U_utR_10;
    const complex_t IT_0088 = lpp_212*U_utR_20;
    const complex_t IT_0089 = IT_0086 + IT_0087 + IT_0088;
    const complex_t IT_0090 = -IT_0089;
    const complex_t IT_0091 = IT_0085 + IT_0090;
    const complex_t IT_0092 = IT_0081*IT_0091;
    const complex_t IT_0093 = std::pow(m_t, 2);
    const complex_t IT_0094 = std::pow(m_utR1, 2);
    const complex_t IT_0095 = std::pow((-2)*s_24 + (complex_t{0, 1})*GuRt
      *m_utR1 + IT_0011 + IT_0093 + -IT_0094, -1);
    const complex_t IT_0096 = (complex_t{0, 1})*IT_0092*IT_0095;
    const complex_t IT_0097 = gw*gwuR*std::conj(U_utR_21);
    const complex_t IT_0098 = lpp_021*U_utR_01;
    const complex_t IT_0099 = lpp_121*U_utR_11;
    const complex_t IT_0100 = lpp_221*U_utR_21;
    const complex_t IT_0101 = IT_0098 + IT_0099 + IT_0100;
    const complex_t IT_0102 = lpp_012*U_utR_01;
    const complex_t IT_0103 = lpp_112*U_utR_11;
    const complex_t IT_0104 = lpp_212*U_utR_21;
    const complex_t IT_0105 = IT_0102 + IT_0103 + IT_0104;
    const complex_t IT_0106 = -IT_0105;
    const complex_t IT_0107 = IT_0101 + IT_0106;
    const complex_t IT_0108 = IT_0097*IT_0107;
    const complex_t IT_0109 = std::pow(m_utR2, 2);
    const complex_t IT_0110 = std::pow((-2)*s_24 + (complex_t{0, 1})*GcRt
      *m_utR2 + IT_0011 + IT_0093 + -IT_0109, -1);
    const complex_t IT_0111 = (complex_t{0, 1})*IT_0108*IT_0110;
    const complex_t IT_0112 = gw*gwuR*std::conj(U_utR_22);
    const complex_t IT_0113 = lpp_021*U_utR_02;
    const complex_t IT_0114 = lpp_121*U_utR_12;
    const complex_t IT_0115 = lpp_221*U_utR_22;
    const complex_t IT_0116 = IT_0113 + IT_0114 + IT_0115;
    const complex_t IT_0117 = lpp_012*U_utR_02;
    const complex_t IT_0118 = lpp_112*U_utR_12;
    const complex_t IT_0119 = lpp_212*U_utR_22;
    const complex_t IT_0120 = IT_0117 + IT_0118 + IT_0119;
    const complex_t IT_0121 = -IT_0120;
    const complex_t IT_0122 = IT_0116 + IT_0121;
    const complex_t IT_0123 = IT_0112*IT_0122;
    const complex_t IT_0124 = std::pow(m_utR3, 2);
    const complex_t IT_0125 = std::pow((-2)*s_24 + (complex_t{0, 1})*GtRt
      *m_utR3 + IT_0011 + IT_0093 + -IT_0124, -1);
    const complex_t IT_0126 = (complex_t{0, 1})*IT_0123*IT_0125;
    const complex_t IT_0127 = -IT_0096 + -IT_0111 + -IT_0126;
    const complex_t IT_0128 = -IT_0040;
    const complex_t IT_0129 = IT_0079 + IT_0128;
    const complex_t IT_0130 = 2*IT_0075*(std::conj(IT_0075)*IT_0076 + (-0.5)
      *std::conj(IT_0039)*IT_0080 + 0.5*std::conj(IT_0127)*(IT_0076 + IT_0077 +
       IT_0128)) + 2*IT_0127*(IT_0077*std::conj(IT_0127) + 0.5*std::conj(IT_0075
      )*(IT_0076 + IT_0077 + IT_0128) + (-0.5)*std::conj(IT_0039)*IT_0129) + 2
      *IT_0039*(std::conj(IT_0039)*IT_0040 + (-0.5)*std::conj(IT_0075)*IT_0080 +
       (-0.5)*std::conj(IT_0127)*IT_0129);
    return IT_0130;
}
} // End of namespace brparity
