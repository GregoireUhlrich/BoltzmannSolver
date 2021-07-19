#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_dtR1_to_uf3c_df1c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_dtR1_to_uf3c_df1c(
        param_t const &param
        )
{
    clearcache();
    const real_t s_23 = param.s_23;
    const complex_t lpp_201 = param.lpp_201;
    const complex_t lpp_202 = param.lpp_202;
    const complex_t lpp_210 = param.lpp_210;
    const complex_t lpp_220 = param.lpp_220;
    const complex_t U_dtR_10 = param.U_dtR_10;
    const complex_t U_dtR_20 = param.U_dtR_20;
    const complex_t IT_0000 = lpp_201*U_dtR_10;
    const complex_t IT_0001 = lpp_202*U_dtR_20;
    const complex_t IT_0002 = lpp_210*U_dtR_10;
    const complex_t IT_0003 = lpp_220*U_dtR_20;
    const complex_t IT_0004 = (complex_t{0, -1})*(IT_0000 + IT_0001 + -IT_0002
       + -IT_0003);
    return 4*s_23*IT_0004*std::conj(IT_0004);
}
} // End of namespace brparity