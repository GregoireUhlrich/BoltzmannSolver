#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_xx_to_uf2c_df2c_df1c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XX_to_uf2c_df2c_df1c(
        param_t const &param
        )
{
    clearcache();
    auto const &gw = param.gw;
    auto const &mX = param.mX;
    auto const &m_c = param.m_c;
    auto const &m_d = param.m_d;
    auto const &m_s = param.m_s;
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
    auto const &lpp_001 = param.lpp_001;
    auto const &lpp_010 = param.lpp_010;
    auto const &lpp_101 = param.lpp_101;
    auto const &lpp_102 = param.lpp_102;
    auto const &lpp_110 = param.lpp_110;
    auto const &lpp_112 = param.lpp_112;
    auto const &lpp_120 = param.lpp_120;
    auto const &lpp_121 = param.lpp_121;
    auto const &lpp_201 = param.lpp_201;
    auto const &lpp_210 = param.lpp_210;
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
    const complex_t IT_0000 = gw*gwuR*std::conj(U_utR_10);
    const complex_t IT_0001 = lpp_010*U_utR_00;
    const complex_t IT_0002 = lpp_110*U_utR_10;
    const complex_t IT_0003 = lpp_210*U_utR_20;
    const complex_t IT_0004 = IT_0001 + IT_0002 + IT_0003;
    const complex_t IT_0005 = lpp_001*U_utR_00;
    const complex_t IT_0006 = lpp_101*U_utR_10;
    const complex_t IT_0007 = lpp_201*U_utR_20;
    const complex_t IT_0008 = IT_0005 + IT_0006 + IT_0007;
    const complex_t IT_0009 = -IT_0008;
    const complex_t IT_0010 = IT_0004 + IT_0009;
    const complex_t IT_0011 = IT_0000*IT_0010;
    const complex_t IT_0012 = std::pow(m_c, 2);
    const complex_t IT_0013 = 2*mX;
    const complex_t IT_0014 = std::pow(IT_0013, 2);
    const complex_t IT_0015 = std::pow(m_utR1, 2);
    const complex_t IT_0016 = std::pow((-2)*s_23 + (complex_t{0, 1})*GuRt
      *m_utR1 + IT_0012 + IT_0014 + -IT_0015, -1);
    const complex_t IT_0017 = (complex_t{0, 1})*IT_0011*IT_0016;
    const complex_t IT_0018 = gw*gwuR*std::conj(U_utR_11);
    const complex_t IT_0019 = lpp_010*U_utR_01;
    const complex_t IT_0020 = lpp_110*U_utR_11;
    const complex_t IT_0021 = lpp_210*U_utR_21;
    const complex_t IT_0022 = IT_0019 + IT_0020 + IT_0021;
    const complex_t IT_0023 = lpp_001*U_utR_01;
    const complex_t IT_0024 = lpp_101*U_utR_11;
    const complex_t IT_0025 = lpp_201*U_utR_21;
    const complex_t IT_0026 = IT_0023 + IT_0024 + IT_0025;
    const complex_t IT_0027 = -IT_0026;
    const complex_t IT_0028 = IT_0022 + IT_0027;
    const complex_t IT_0029 = IT_0018*IT_0028;
    const complex_t IT_0030 = std::pow(m_utR2, 2);
    const complex_t IT_0031 = std::pow((-2)*s_23 + (complex_t{0, 1})*GcRt
      *m_utR2 + IT_0012 + IT_0014 + -IT_0030, -1);
    const complex_t IT_0032 = (complex_t{0, 1})*IT_0029*IT_0031;
    const complex_t IT_0033 = gw*gwuR*std::conj(U_utR_12);
    const complex_t IT_0034 = lpp_010*U_utR_02;
    const complex_t IT_0035 = lpp_110*U_utR_12;
    const complex_t IT_0036 = lpp_210*U_utR_22;
    const complex_t IT_0037 = IT_0034 + IT_0035 + IT_0036;
    const complex_t IT_0038 = lpp_001*U_utR_02;
    const complex_t IT_0039 = lpp_101*U_utR_12;
    const complex_t IT_0040 = lpp_201*U_utR_22;
    const complex_t IT_0041 = IT_0038 + IT_0039 + IT_0040;
    const complex_t IT_0042 = -IT_0041;
    const complex_t IT_0043 = IT_0037 + IT_0042;
    const complex_t IT_0044 = IT_0033*IT_0043;
    const complex_t IT_0045 = std::pow(m_utR3, 2);
    const complex_t IT_0046 = std::pow((-2)*s_23 + (complex_t{0, 1})*GtRt
      *m_utR3 + IT_0012 + IT_0014 + -IT_0045, -1);
    const complex_t IT_0047 = (complex_t{0, 1})*IT_0044*IT_0046;
    const complex_t IT_0048 = -IT_0017 + -IT_0032 + -IT_0047;
    const complex_t IT_0049 = s_23*s_45;
    const complex_t IT_0050 = lpp_101*U_dtR_10;
    const complex_t IT_0051 = lpp_102*U_dtR_20;
    const complex_t IT_0052 = IT_0050 + IT_0051;
    const complex_t IT_0053 = lpp_110*U_dtR_10;
    const complex_t IT_0054 = lpp_120*U_dtR_20;
    const complex_t IT_0055 = -IT_0053 + -IT_0054;
    const complex_t IT_0056 = IT_0052 + IT_0055;
    const complex_t IT_0057 = gw*gwdR*std::conj(U_dtR_10);
    const complex_t IT_0058 = IT_0056*IT_0057;
    const complex_t IT_0059 = std::pow(m_s, 2);
    const complex_t IT_0060 = std::pow(m_dtR1, 2);
    const complex_t IT_0061 = std::pow((-2)*s_24 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0014 + IT_0059 + -IT_0060, -1);
    const complex_t IT_0062 = (complex_t{0, 1})*IT_0058*IT_0061;
    const complex_t IT_0063 = lpp_101*U_dtR_11;
    const complex_t IT_0064 = lpp_102*U_dtR_21;
    const complex_t IT_0065 = IT_0063 + IT_0064;
    const complex_t IT_0066 = lpp_110*U_dtR_11;
    const complex_t IT_0067 = lpp_120*U_dtR_21;
    const complex_t IT_0068 = -IT_0066 + -IT_0067;
    const complex_t IT_0069 = IT_0065 + IT_0068;
    const complex_t IT_0070 = gw*gwdR*std::conj(U_dtR_11);
    const complex_t IT_0071 = IT_0069*IT_0070;
    const complex_t IT_0072 = std::pow(m_dtR2, 2);
    const complex_t IT_0073 = std::pow((-2)*s_24 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0014 + IT_0059 + -IT_0072, -1);
    const complex_t IT_0074 = (complex_t{0, 1})*IT_0071*IT_0073;
    const complex_t IT_0075 = lpp_101*U_dtR_12;
    const complex_t IT_0076 = lpp_102*U_dtR_22;
    const complex_t IT_0077 = IT_0075 + IT_0076;
    const complex_t IT_0078 = lpp_110*U_dtR_12;
    const complex_t IT_0079 = lpp_120*U_dtR_22;
    const complex_t IT_0080 = -IT_0078 + -IT_0079;
    const complex_t IT_0081 = IT_0077 + IT_0080;
    const complex_t IT_0082 = gw*gwdR*std::conj(U_dtR_12);
    const complex_t IT_0083 = IT_0081*IT_0082;
    const complex_t IT_0084 = std::pow(m_dtR3, 2);
    const complex_t IT_0085 = std::pow((-2)*s_24 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0014 + IT_0059 + -IT_0084, -1);
    const complex_t IT_0086 = (complex_t{0, 1})*IT_0083*IT_0085;
    const complex_t IT_0087 = IT_0062 + IT_0074 + IT_0086;
    const complex_t IT_0088 = s_24*s_35;
    const complex_t IT_0089 = s_25*s_34;
    const complex_t IT_0090 = -IT_0089;
    const complex_t IT_0091 = IT_0088 + IT_0090;
    const complex_t IT_0092 = lpp_101*U_dtR_00;
    const complex_t IT_0093 = lpp_121*U_dtR_20;
    const complex_t IT_0094 = IT_0092 + IT_0093;
    const complex_t IT_0095 = lpp_110*U_dtR_00;
    const complex_t IT_0096 = lpp_112*U_dtR_20;
    const complex_t IT_0097 = -IT_0095 + -IT_0096;
    const complex_t IT_0098 = IT_0094 + IT_0097;
    const complex_t IT_0099 = gw*gwdR*std::conj(U_dtR_00);
    const complex_t IT_0100 = IT_0098*IT_0099;
    const complex_t IT_0101 = std::pow(m_d, 2);
    const complex_t IT_0102 = std::pow((-2)*s_25 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0014 + -IT_0060 + IT_0101, -1);
    const complex_t IT_0103 = (complex_t{0, 1})*IT_0100*IT_0102;
    const complex_t IT_0104 = lpp_101*U_dtR_01;
    const complex_t IT_0105 = lpp_121*U_dtR_21;
    const complex_t IT_0106 = IT_0104 + IT_0105;
    const complex_t IT_0107 = lpp_110*U_dtR_01;
    const complex_t IT_0108 = lpp_112*U_dtR_21;
    const complex_t IT_0109 = -IT_0107 + -IT_0108;
    const complex_t IT_0110 = IT_0106 + IT_0109;
    const complex_t IT_0111 = gw*gwdR*std::conj(U_dtR_01);
    const complex_t IT_0112 = IT_0110*IT_0111;
    const complex_t IT_0113 = std::pow((-2)*s_25 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0014 + -IT_0072 + IT_0101, -1);
    const complex_t IT_0114 = (complex_t{0, 1})*IT_0112*IT_0113;
    const complex_t IT_0115 = lpp_101*U_dtR_02;
    const complex_t IT_0116 = lpp_121*U_dtR_22;
    const complex_t IT_0117 = IT_0115 + IT_0116;
    const complex_t IT_0118 = lpp_110*U_dtR_02;
    const complex_t IT_0119 = lpp_112*U_dtR_22;
    const complex_t IT_0120 = -IT_0118 + -IT_0119;
    const complex_t IT_0121 = IT_0117 + IT_0120;
    const complex_t IT_0122 = gw*gwdR*std::conj(U_dtR_02);
    const complex_t IT_0123 = IT_0121*IT_0122;
    const complex_t IT_0124 = std::pow((-2)*s_25 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0014 + -IT_0084 + IT_0101, -1);
    const complex_t IT_0125 = (complex_t{0, 1})*IT_0123*IT_0124;
    const complex_t IT_0126 = IT_0103 + IT_0114 + IT_0125;
    const complex_t IT_0127 = -IT_0049;
    const complex_t IT_0128 = IT_0091 + IT_0127;
    const complex_t IT_0129 = IT_0088 + IT_0089 + IT_0127;
    const complex_t IT_0130 = (-6)*IT_0129;
    const complex_t IT_0131 = IT_0048*(12*std::conj(IT_0048)*IT_0049 + (-6)
      *std::conj(IT_0087)*(IT_0049 + IT_0091) + 6*std::conj(IT_0126)*IT_0128) +
       IT_0126*(12*IT_0089*std::conj(IT_0126) + 6*std::conj(IT_0048)*IT_0128 +
       std::conj(IT_0087)*IT_0130) + IT_0087*(12*std::conj(IT_0087)*IT_0088 + (
      -6)*std::conj(IT_0048)*(IT_0049 + IT_0091) + std::conj(IT_0126)*IT_0130);
    return IT_0131;
}
} // End of namespace brparity
