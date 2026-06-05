/***
The contents of this file were copied from the authors of PLINK with our appreciation.
The original file contains the following copyright

// This file is part of PLINK 1.9, copyright (C) 2005-2025 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

and can be found at https://github.com/chrchang/plink-ng/blob/master/1.9/plink_stats.c
***/


#include "stats_test.hpp"
#include "utils.hpp"

// using boost::multiprecision::cpp_dec_float_50;
// using boost::math::chi_squared_distribution;

// Fisher's Exact Test for 2x2 contingency table
#define DBL_MAX 1.7976931348623157e308

namespace stoat {
// ------------------------ Fisher exact test ------------------------

// Fisher's exact test for a 2x2 contingency table
// m11, m12, m21, m22 are the counts in the table
double FisherChi2::fastFishersExactTest(size_t m11, size_t m12,
                                 size_t m21, size_t m22) const {
    
    // Check for any full-zero row or column
    if ((m11 | m12) == 0 || (m21 | m22) == 0 || (m11 | m21) == 0 || (m12 | m22) == 0) {
        return std::nan("");
    }
    
    double tprob = (1 - kExactTestEpsilon2) * kExactTestBias;
    double cur_prob = tprob;
    double cprob = 0;
    size_t uii;
    double cur11, cur12, cur21, cur22;
    double preaddp;

    // Ensure we are left of the distribution center, m11 <= m22, and m12 <= m21.
    if (m12 > m21) {
        uii = m12;
        m12 = m21;
        m21 = uii;
    }

    if (m11 > m22) {
        uii = m11;
        m11 = m22;
        m22 = uii;
    }

    if ((static_cast<size_t>(m11) * m22) > (static_cast<size_t>(m12) * m21)) {
        uii = m11;
        m11 = m12;
        m12 = uii;
        uii = m21;
        m21 = m22;
        m22 = uii;
    }

    cur11 = m11;
    cur12 = m12;
    cur21 = m21;
    cur22 = m22;

    while (cur12 > 0.5) {
        cur11 += 1;
        cur22 += 1;
        cur_prob *= (cur12 * cur21) / (cur11 * cur22);
        cur12 -= 1;
        cur21 -= 1;
        if (cur_prob > DBL_MAX) {
        return 0.0;
        }
        if (cur_prob < kExactTestBias) {
        tprob += cur_prob;
        break;
        }
        cprob += cur_prob;
    }

    if (cprob == 0) {
        return 1.0;
    }

    while (cur12 > 0.5) {
        cur11 += 1;
        cur22 += 1;
        cur_prob *= (cur12 * cur21) / (cur11 * cur22);
        cur12 -= 1;
        cur21 -= 1;
        preaddp = tprob;
        tprob += cur_prob;
        if (tprob <= preaddp) {
        break;
        }
    }

    if (m11) {
        cur11 = m11;
        cur12 = m12;
        cur21 = m21;
        cur22 = m22;
        cur_prob = (1 - kExactTestEpsilon2) * kExactTestBias;
        do {
        cur12 += 1;
        cur21 += 1;
        cur_prob *= (cur11 * cur22) / (cur12 * cur21);
        cur11 -= 1;
        cur22 -= 1;
        preaddp = tprob;
        tprob += cur_prob;
        if (tprob <= preaddp) {
            return preaddp / (cprob + preaddp);
        }
        } while (cur11 > 0.5);
    }

    return tprob / (cprob + tprob);
}

}
