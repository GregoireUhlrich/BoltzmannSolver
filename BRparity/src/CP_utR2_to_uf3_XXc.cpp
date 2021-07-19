#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_utR2_to_uf3_XXc.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_utR2_to_uf3_XXc(
        param_t const &param
        )
{
    clearcache();
    const real_t gw = param.gw;
    const real_t gwuR = param.gwuR;
    const real_t s_23 = param.s_23;
    const complex_t U_utR_21 = param.U_utR_21;
    const complex_t IT_0000 = (complex_t{0, -1})*gw*gwuR*U_utR_21;
    return 2*s_23*IT_0000*std::conj(IT_0000);
}
} // End of namespace brparity