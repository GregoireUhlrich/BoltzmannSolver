#include <fstream>
#include <iostream>
#include <istream>
#include <marty.h>
#include <stdlib.h>
#include <string>

using namespace std;
using namespace csl;
using namespace mty;

// To change by the user !
std::string path_to_generated_library = ".";

Expr cc(Expr const &expr) { return GetComplexConjugate(expr); }

int main()
{
    // Model building

    // Dirac space
    Index a = DiracIndex();
    Index b = DiracIndex();

    Model toyModel;
    toyModel.addGaugedGroup(group::Type::SU, "L", 2);
    toyModel.addGaugedGroup(group::Type::U1, "Y");
    toyModel.addFlavorGroup("SM_flavor", 3, true);
    toyModel.init();

    toyModel.renameParticle("A_L", "W");
    toyModel.renameParticle("A_Y", "B");

    Particle et = scalarboson_s("et", toyModel);
    et->setGroupRep("L", {1});
    et->setGroupRep("Y", {-1, 2});
    toyModel.addParticle(et);

    Particle Li = diracfermion_s("L_L", toyModel);
    Li->setGroupRep("L", {1});
    Li->setGroupRep("Y", {-1, 2});
    Li->setFundamentalFlavorRep("SM_flavor");
    toyModel.addParticle(Li);

    Particle Ni = weylfermion_s("N", toyModel, Chirality::Right);
    Ni->setGroupRep("Y", {0});
    Ni->setFundamentalFlavorRep("SM_flavor");
    toyModel.addParticle(Ni);

    Particle SS = scalarboson_s("SS", toyModel);
    SS->setGroupRep("Y", {0});
    SS->setSelfConjugate(true);
    toyModel.addParticle(SS);

    // SM flavor space
    Index I = toyModel.generateIndex("SM_flavor", "L_L");
    Index J = toyModel.generateIndex("SM_flavor", "L_L");

    // SU(2)L space
    Index i = toyModel.generateIndex("L", "L_L");
    Index j = toyModel.generateIndex("L", "L_L");
    Tensor eps = i.getSpace()->getEpsilon();

    // Tensor pauli = toyModel.getGenerator("L", "L_L");
    Tensor pauli = GetGenerator(toyModel, "L", "L_L");
    auto ii = toyModel.generateIndices(3, "L", "L_L");

    auto flavorSpace = toyModel.getVectorSpace("SM_flavor", "L_L");
    Tensor Ye = csl::Tensor("Y_l", {flavorSpace, flavorSpace});
    Ye->setComplexProperty(csl::ComplexProperty::Complex);

    Tensor Yet = csl::Tensor("Y_et", {flavorSpace, flavorSpace});
    Yet->setComplexProperty(csl::ComplexProperty::Complex);

    Tensor Ynu = csl::Tensor("Y_nu", {flavorSpace, flavorSpace});
    Ynu->setComplexProperty(csl::ComplexProperty::Complex);

    Tensor YS = csl::Tensor("Y_S", {flavorSpace, flavorSpace});
    YS->setComplexProperty(csl::ComplexProperty::Complex);

    Tensor MN = csl::Tensor("M_N", {flavorSpace, flavorSpace});
    MN->setComplexProperty(csl::ComplexProperty::Real);

    Expr muet2 = constant_s("muet2");
    Expr muS2 = constant_s("muS2");
    Expr muSet = constant_s("muSet");
    muSet->setComplexProperty(csl::ComplexProperty::Complex);
    Expr let = constant_s("lamet");
    Expr lS = constant_s("lS");
    Expr lSet = constant_s("lSet");

    auto al = DiracIndices(3);
    Tensor C = dirac4.C_matrix;
    Tensor PL = dirac4.P_L;
    Tensor PR = dirac4.P_R;

    toyModel.addLagrangianTerm(-Ni({I, al[0]}) * C({al[0], al[1]}) * PL({al[1], al[2]}) * Yet({I, J}) * cc(et(i)) * Li({J, i, al[2]}),
                               true // Add also the complex conjugate of this term
    );

    toyModel.addLagrangianTerm(-Ni({I, al[0]}) * C({al[0], al[1]}) * MN({I, J}) * Ni({J, al[1]}), true // Add also the complex conjugate of this term
    );

    toyModel.addLagrangianTerm(-SS * Ni({I, al[0]}) * C({al[0], al[1]}) * YS({I, J}) * Ni({J, al[1]}));

    toyModel.addLagrangianTerm(muet2 * cc(et(i)) * et(i));

    toyModel.addLagrangianTerm(muS2 * SS * SS);

    toyModel.addLagrangianTerm(-lS * SS * SS * SS * SS);

    Expr prod1 = cc(et(i)) * et(i);
    toyModel.addLagrangianTerm(let * prod1 * RenamedIndices(prod1));

    toyModel.addLagrangianTerm(lSet * SS * SS * cc(et(i)) * et(i));

    toyModel.addLagrangianTerm(muSet * SS * cc(et(i)) * et(i));

    cout << toyModel << endl;

    Expr m_N_1 = constant_s("m_N_1");
    Expr m_N_2 = constant_s("m_N_2");
    Expr m_N_3 = constant_s("m_N_3");

    csl::Tensor M_N = csl::tensor_s("M_N", {flavorSpace, flavorSpace}, csl::matrix_s({{m_N_1, 0, 0}, {0, m_N_2, 0}, {0, 0, m_N_3}}));

    Index f_i = GetIndex(flavorSpace);
    Index f_j = GetIndex(flavorSpace);
    toyModel.replace(MN, M_N({f_i, f_j}));

    toyModel.breakFlavorSymmetry("SM_flavor");
    toyModel.renameParticle("L_L_1", "lL_1; e_L");
    toyModel.renameParticle("L_L_2", "lL_2 ; \\mu_L");
    toyModel.renameParticle("L_L_3", "lL_3 ; \\tau_L");

    toyModel.getParticle("N_1")->setSelfConjugate(true);
    toyModel.getParticle("N_2")->setSelfConjugate(true);
    toyModel.getParticle("N_3")->setSelfConjugate(true);
    toyModel.promoteToMajorana("N_1");
    toyModel.promoteToMajorana("N_2");

    Expr thetaW = constant_s("thetaW");
    Expr cW = cos_s(thetaW);
    Expr sW = sin_s(thetaW);

    Expr e = constant_s("e");
    Expr gY = toyModel.getScalarCoupling("g_Y");
    Expr gL = toyModel.getScalarCoupling("g_L");
    toyModel.replace(gY, e / cW);
    toyModel.replace(gL, e / sW);

    // Expr M_W = constant_s("M_W");
    // Expr M_Z = constant_s("M_Z");

    toyModel.refresh();
    std::cout << toyModel << std::endl;

    ///////////////////////////////////////////////////
    //
    ///////////////////////////////////////////////////

    auto rules = toyModel.getFeynmanRules();
    Display(rules); // Displays expressions in terminal
    // Show(rules); // Shows diagrams in the application
    //--------------------------

    mty::Library myLib("LeptoProto_test", path_to_generated_library);
    mty::option::decomposeInOperators = true;
    mty::option::decomposeInLocalOperator = false;

    FeynOptions options;
    options.addFilters(
            [&](FeynmanDiagram const &diag) {
                return !diag.isInLoop("lL_1")
                    && !diag.isInLoop("lL_2")
                    && !diag.isInLoop("N_1")
                    && !diag.isInLoop("N_2")
                    // && !diag.isInLoop("et")
                && !diag.isInLoop("lL_3")
                && !diag.isInLoop("W")
                && !diag.isInLoop("B");
            }
        );
    options.setTopology(mty::Topology::Mass);
    auto loop_ampl = toyModel.computeAmplitude(
        mty::OneLoop,
        {Incoming("et"), Incoming("lL_1"), Outgoing("et"), Outgoing("lL_1")},
        options
    );
    auto tree_ampl = toyModel.computeAmplitude(
        mty::TreeLevel,
        {Incoming("et"), Incoming("lL_1"), Outgoing("et"), Outgoing("lL_1")}
    );
    Display(loop_ampl);
    // Show(loop_ampl);
    Display(tree_ampl);
    return 0;
    auto squared = toyModel.computeSquaredAmplitude(tree_ampl, loop_ampl);
    myLib.addFunction("test", squared);

    ///////////////////////////////////////////////////
    // Saving JSON data to test.json
    ///////////////////////////////////////////////////

    // myLib.applyDiagonalizationData(toyModel);
    myLib.generateSpectrum(toyModel);

    myLib.cleanExistingSources();
    myLib.setGccCompiler();
    myLib.build(4);
    myLib.print();

    return 0;
}
