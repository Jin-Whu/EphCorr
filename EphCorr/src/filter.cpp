#include <algorithm>
#include <numeric>

#include "filter.h"

void filterOutlier(std::vector<double> &vec)
{
    double median = 0;
    int nums = vec.size();
    std::sort(vec.begin(), vec.end());

    if (nums % 2)
        median = vec[(nums - 2) / 2] + vec[nums / 2];
    else
        median = vec[nums / 2];

    std::vector<double> mads;
    for (const auto &num : vec)
    {
        mads.push_back(abs(num - median));
    }

    std::sort(mads.begin(), mads.end());

    double mad = 0;
    if (nums % 2)
        mad = mads[(nums - 2) / 2] + mads[nums / 2];
    else
        mad = mads[nums / 2];

    std::vector<int> outliersIndex;

    auto finder = [&](const double &x) {
        double value = 0.6745 * abs(x - median) / mad;
        return value > 3.5;
    };
    vec.erase(std::remove_if(vec.begin(), vec.end(), finder), vec.end());

}
