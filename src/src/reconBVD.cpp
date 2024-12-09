#include "reconBVD.hpp"
#include "interScheme.hpp"
#include <fstream>

void ReconerBVD::boundaryDimishing()
{
    std::copy(candidateValues2.begin(), candidateValues2.end(), data->begin());
    iter = data->begin() + 4 * nvar;
    candidateIter = candidateValues.begin() + 4 * nvar;
    candidate2Iter = candidateValues2.begin() + 4 * nvar;
    // candidate3Iter = candidateValues3.begin() + 4 * nvar;

    auto iterEdn = data->end() - 4 * nvar;
    auto candidateIterEdn = candidateValues.end() - 4 * nvar;
    auto candidate2IterEdn = candidateValues2.end() - 4 * nvar;
    // auto candidate3IterEdn = candidateValues3.end() - 4 * nvar;

    // std::ofstream file("output.dat", std::ios::out);
    // int count = 0;
    // for (auto ivar = candidateValues.begin(); ivar != candidateValues.end(); ivar += 2 * nvar) {
    //     file << std::format("{} {} {} {} {} {} {} \n", count++, ivar[0], ivar[1], ivar[2], ivar[3], ivar[4], ivar[5]);
    // }
    // file << std::endl;
    // file.close();

    // std::ofstream file2("output2.dat", std::ios::out);
    // count = 0;
    // for (auto ivar = candidateValues2.begin(); ivar != candidateValues2.end(); ivar += 2 * nvar) {
    //     file2 << std::format("{} {} {} {} {} {} {} \n", count++, ivar[0], ivar[1], ivar[2], ivar[3], ivar[4], ivar[5]);
    // }
    // file2 << std::endl;
    // file2.close();

    // std::ofstream file3("output3.dat", std::ios::out);
    // count = 0;
    // for (auto ivar = candidateValues3.begin(); ivar != candidateValues3.end(); ivar += 2 * nvar) {
    //     file3 << std::format("{} {} {} {} {} {} {} \n", count++, ivar[0], ivar[1], ivar[2], ivar[3], ivar[4], ivar[5]);
    // }
    // file3 << std::endl;
    // file3.close();

    for (; iter != iterEdn; iter += nvar, candidateIter += nvar, candidate2Iter += nvar) { //, candidate3Iter += nvar
        for (int i = 0; i < nvar; i++) {
            real bvd1 = std::abs(*(candidate2Iter - 2 * nvar) - *(candidate2Iter - nvar))
                + std::abs(*(candidate2Iter + 0 * nvar) - *(candidate2Iter + 1 * nvar));

            real bvd2 = std::abs(*(candidateIter - 2 * nvar) - *(candidateIter - nvar))
                + std::abs(*(candidateIter + 0 * nvar) - *(candidateIter + 1 * nvar));

            // real bvd3 = std::abs(*(candidate3Iter - 2 * nvar) - *(candidate3Iter - nvar))
            //     + std::abs(*(candidate3Iter + 0 * nvar) - *(candidate3Iter + 1 * nvar));
            if (bvd2 < bvd1) {
                // if (true) {
                // if (bvd2 < bvd3) {
                // *(iter - 3 * nvar) = *(candidateIter - 3 * nvar);
                // *(iter - 2 * nvar) = *(candidateIter - 2 * nvar);
                *(iter - 1 * nvar) = *(candidateIter - 1 * nvar);
                *(iter + 0 * nvar) = *(candidateIter + 0 * nvar);
                // *(iter + 1 * nvar) = *(candidateIter + 1 * nvar);
                // *(iter + 2 * nvar) = *(candidateIter + 2 * nvar);
                // } else {
                //     // *(iter - 3 * nvar) = *(candidate3Iter - 3 * nvar);
                //     // *(iter - 2 * nvar) = *(candidate3Iter - 2 * nvar);
                //     *(iter - 1 * nvar) = *(candidate3Iter - 1 * nvar);
                //     *(iter + 0 * nvar) = *(candidate3Iter + 0 * nvar);
                //     // *(iter + 1 * nvar) = *(candidate3Iter + 1 * nvar);
                //     // *(iter + 2 * nvar) = *(candidate3Iter + 2 * nvar);
                // }
            }

            // if (bvd2 > bvd1) {
            //     *(iter - 1 * nvar) = *(candidate2Iter - 1 * nvar);
            //     *(iter + 0 * nvar) = *(candidate2Iter + 0 * nvar);
            // }
            iter++;
            candidateIter++;
            candidate2Iter++;
            // candidate3Iter++;
        }
        // std::cout << candidate3IterEdn - candidate3Iter << '\n';
    }

    // std::ofstream file4("output4.dat", std::ios::out);
    // count = 0;
    // for (auto ivar = data->begin(); ivar != data->end(); ivar += 2 * nvar) {
    //     file4 << std::format("{} {} {} {} {} {} {} \n", count++, ivar[0], ivar[1], ivar[2], ivar[3], ivar[4], ivar[5]);
    // }
    // file4 << std::endl;
    // file4.close();

    assert(iter == iterEdn);
    assert(candidateIter == candidateIterEdn);
    // assert(candidate3Iter == candidate3IterEdn);
}
void ReconerBVD::initVarsR()
{
    int nPointR = bndL->getN() + n + bndR->getN() - 5;
    if (!data) {
        data = std::make_shared<Data>(nPointR, nvar * 2);
    }
    if (candidateValues.size() == 0) {
        candidateValues.init(nPointR, nvar * 2);
    }
    if (candidateValues2.size() == 0) {
        candidateValues2.init(nPointR, nvar * 2);
    }
    // if (candidateValues3.size() == 0) {
    //     candidateValues3.init(nPointR, nvar * 2);
    // }
}

void ReconerBVD::reconI(std::array<std::vector<real>::iterator, 6> input)
{
    for (int i = 0; i < nvar; i++) {
        // *candidate2Iter = Linear5({
        //     input[0][i],
        //     input[1][i],
        //     input[2][i],
        //     input[3][i],
        //     input[4][i],
        // });
        // *(candidate2Iter + nvar) = Linear5({
        //     input[5][i],
        //     input[4][i],
        //     input[3][i],
        //     input[2][i],
        //     input[1][i],
        // });

        // for candidate

        *candidate2Iter = THINC2(input[1][i],
            input[2][i],
            input[3][i]);
        *(candidate2Iter + nvar) = THINC2(input[4][i],
            input[3][i],
            input[2][i]);

        *candidateIter = THINC1(input[1][i],
            input[2][i],
            input[3][i]);
        *(candidateIter + nvar) = THINC1(input[4][i],
            input[3][i],
            input[2][i]);

        // *candidate3Iter = THINC2(input[1][i],
        //     input[2][i],
        //     input[3][i]);
        // *(candidate3Iter + nvar) = THINC2(input[4][i],
        //     input[3][i],
        //     input[2][i]);

        candidate2Iter++;
        // candidate3Iter++;

        candidateIter++;
    }
    // candidate3Iter += nvar;
    candidate2Iter += nvar;
    candidateIter += nvar;
    iter += 2 * nvar;
}