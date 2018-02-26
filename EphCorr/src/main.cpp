#include <numeric>
#include <fstream>
#include <sstream>

#include "cxxopts.hpp"
#include "rtklib.h"
#include <filesystem/path.h>

#include "filter.h"


struct Options
{
    std::string strPath;
    int iInterval;
    std::string strSys;
    double eps[6]; //start epoch
    double epe[6]; // end epoc
};

struct Corr
{
    double corr[4];
    std::string satID;
};

struct Corrs
{
    std::vector<double> corr[4]; // dr, dt, dn, dclk
};


bool operator ==(const Corr &corr1, const Corr&corr2)
{
    return corr1.satID == corr2.satID;
}

/*yyyymmdd to epoch[6]*/
void str2epoch(const char *str, double *ep)
{
    int year, month, day;
    std::string date{ str };

    year = std::stoi(date.substr(0, 4));
    month = std::stoi(date.substr(4, 2));
    day = std::stoi(date.substr(6, 2));
    ep[0] = year;
    ep[1] = month;
    ep[2] = day;
    ep[3] = 0;
    ep[4] = 0;
    ep[5] = 0;
}

void parseArg(int argc, char* argv[], Options &opt)
{
    cxxopts::Options options("EphCorr", "compare ephemeris corrections");
    options.add_options()
        ("p,path", "Path", cxxopts::value<std::string>())
        ("i,interval", "Interval", cxxopts::value<int>())
        ("s,sys", "GNSS System", cxxopts::value<std::string>())
        ("start", "Start Epoch", cxxopts::value<std::string>())
        ("end", "End epoch", cxxopts::value<std::string>());
    options.parse(argc, argv);

    if (options.count("path"))
    {
        opt.strPath = options["path"].as<std::string>();
    }

    if (options.count("interval"))
    {
        opt.iInterval = options["interval"].as<int>();
    }

    if (options.count("sys"))
    {
        opt.strSys = options["sys"].as<std::string>();
    }


    if (options.count("start"))
    {
        str2epoch(options["start"].as<std::string>().c_str(), opt.eps);
    }

    if (options.count("end"))
    {
        str2epoch(options["end"].as<std::string>().c_str(), opt.epe);
    }

}


void avediff(const std::string &path)
{
    std::ifstream iFile(path);
    if (!iFile)
        return;

    filesystem::path filePath{ path };
    filesystem::path dir = filePath.parent_path().parent_path() / filesystem::path("ave");

    if (!dir.exists())
        filesystem::create_directory(dir);

    std::string fileName = filePath.filename();

    std::string outFileName{ "ave" };
    outFileName.append(fileName);

    filesystem::path outPath = dir / filesystem::path{ outFileName };

    std::ofstream oFile(outPath.str());

    std::map<std::string, std::vector<Corr>> diffcorr;

    std::string line;
    std::string satID;
    int lag;
    double corr[4];
    Corr dcorr{};
    char satMID[7]{};
    while (std::getline(iFile, line))
    {
        std::istringstream ss(line);
        ss >> satID >> lag >> corr[0] >> corr[1] >> corr[2] >> corr[3];
        dcorr.satID = satID;
        std::copy(corr, corr + 4, dcorr.corr);
        snprintf(satMID, 7, "%s%03d", satID.c_str(), lag);
        diffcorr[satMID].push_back(dcorr);
    }



    for (const auto &item : diffcorr)
    {
        std::vector<double> dcorrs[4];

        oFile << item.first.substr(0, 3) << " " << std::stoi(item.first.substr(3, 3)) << " ";
        for (unsigned i = 0; i != 4; ++i)
        {
            std::transform(item.second.begin(), item.second.end(), std::back_inserter(dcorrs[i]), [&](const Corr &x) {return x.corr[i]; });
            filterOutlier(dcorrs[i]);
            std::fill(corr, corr + 4, 0);
            for (const auto &dc : dcorrs[i])
            {
                corr[i] += dc;
            }

            corr[i] /= dcorrs[i].size();
            oFile << corr[i] << " ";
        }
        oFile << "\n";

    }

    iFile.close();
    oFile.close();
    
}

void diff(const std::string &path, const int interval, const double *ep, const std::string &gnss)
{
    int doy = time2doy(epoch2time(ep));
    
    filesystem::path dir = filesystem::path{path}.parent_path() / filesystem::path("diff");

    if (!dir.exists())
        filesystem::create_directory(dir);

    char outFileName[21];
    snprintf(outFileName, 21, "diffcorr%d_%03d_%s", (int)ep[0], doy, gnss.c_str());
    filesystem::path outPath = dir / outFileName;

    std::ifstream iFile(path);
    std::ofstream oFile(outPath.str());

    if (!iFile)
        return;

    double ephCorr[4], epoch[6];
    std::string line, tmp, satID;
    gtime_t cur_gt{}, flagT{};

    std::vector<Corr> corrs;
    bool bFirst = true;
    Corr corr;
    while (std::getline(iFile, line))
    {
        std::istringstream ss(line);
        ss >> epoch[0] >> epoch[1] >> epoch[2] >> epoch[3] >> epoch[4] >> epoch[5] >> satID;
        ss >> tmp >> ephCorr[0] >> ephCorr[1] >> ephCorr[2] >> tmp >> tmp >> tmp >> ephCorr[3];
        std::copy(ephCorr, ephCorr + 4, corr.corr);
        corr.satID = satID;

        cur_gt = epoch2time(epoch);
        if (bFirst)
        {
            flagT = cur_gt;
            bFirst = false;
        }
        int period = static_cast<int>(timediff(cur_gt, flagT));

        if (period == interval)
        {
            corrs.clear();
            corrs.push_back(corr);
            flagT = cur_gt;
        }
        else if (period == 0)
        {
            corrs.push_back(corr);
        }
        else
        {
            auto it = std::find(corrs.begin(), corrs.end(), corr);
            if (it != corrs.end())
            {

                oFile << satID << " " << period << " ";
                for (unsigned i = 0; i != 4; ++i)
                {
                    oFile << abs(ephCorr[i] - it->corr[i]) << " ";
                }

            }
            oFile << "\n";
        }
    }

    iFile.close();
    oFile.close();

    avediff(outPath.str());
    printf("%d-%02d-%02d %s\n", (int)ep[0], (int)ep[1], (int)ep[2], gnss.c_str());
}


void finalGnss(const std::string &dir, const double *startEpoch, const double *endEpoch, const std::string &gnss)
{
    char file[28]{};
    int doy = 0, year = 0;
    double epoch[6]{};
    gtime_t start = epoch2time(startEpoch);
    gtime_t end = epoch2time(endEpoch);

    std::map <std::string, std::map<int, Corrs>> corrs; // <prn, <lag, values>>

    while (timediff(end, start) >= 0)
    {
        time2epoch(start, epoch);
        year = static_cast<int>(epoch[0]);
        doy = time2doy(start);
        start = timeadd(start, 86400);

        snprintf(file, 28, "avediffcorr%d_%03d_%s", year, doy, gnss.c_str());
        filesystem::path path = filesystem::path{ dir } / filesystem::path{ "ave" } / filesystem::path{ file };

        if (!path.exists()) continue;

        std::ifstream iFile(path.str());

        std::string line;
        std::string prn;
        int lag;
        double corr[4]{}; // dr, dt, dn, dclk
        while (std::getline(iFile, line))
        {
            std::istringstream ss(line);
            ss >> prn >> lag >> corr[0] >> corr[1] >> corr[2] >> corr[3];

            for (unsigned i = 0; i != 4; ++i)
            {
                corrs[prn][lag].corr[i].push_back(corr[i]);
            }
        }

        iFile.close();
    }

    std::string outFile{ "ave" };
    outFile.append(gnss);
    filesystem::path outPath = filesystem::path{dir} / filesystem::path{ outFile };
    std::ofstream out(outPath.str());
    out.setf(std::ios::fixed, std::ios::floatfield);
    out.precision(3);

    std::map<int, Corrs> tCorrs;
    for (auto &it : corrs)
    {
        std::string prn = it.first;
        for (auto &it : it.second)
        {
            int lag = it.first;
            out << prn << " " << lag << " ";

            for (unsigned i = 0; i != 4; ++i)
            {
                filterOutlier(it.second.corr[i]);

                tCorrs[lag].corr[i].insert(tCorrs[lag].corr[i].end(), it.second.corr[i].begin(), it.second.corr[i].end());

                double sum = std::accumulate(it.second.corr[i].begin(), it.second.corr[i].end(), 0.);
                double mean = sum / it.second.corr[i].size();
                out << mean << " ";
            }
            out << '\n';
        }

        out << '\n';
    }

    for (const auto &it : tCorrs)
    {
        int lag = it.first;

        out << lag << " ";
        for (unsigned i = 0; i != 4; ++i)
        {
            double sum = std::accumulate(it.second.corr[i].begin(), it.second.corr[i].end(), 0.);
            double mean = sum / it.second.corr[i].size();

            out << mean << " ";
        }

        for (unsigned i = 0; i != 4; ++i)
        {
            double max = *std::max_element(it.second.corr[i].begin(), it.second.corr[i].end());
            out << max << " ";
        }

        out << '\n';
    }

    out.close();

}

void aveFinal(const Options &opt)
{

    std::vector<std::string> gnss{ "GPS", "BDS", "GAL", "GLO" };
    if (opt.strSys == "ALL")
    {
        for (unsigned i = 0; i != 4; ++i)
        {
            finalGnss(opt.strPath, opt.eps, opt.epe, gnss[i]);
        }
    }
    else
    {
        finalGnss(opt.strPath, opt.eps, opt.epe, opt.strSys);
    }
}

void process(const Options &opt)
{
    char file[20]{};
    int doy = 0, year = 0;
    double epoch[6]{};

    gtime_t start = epoch2time(opt.eps);
    gtime_t end = epoch2time(opt.epe);

    std::vector<std::string> gnss{ "GPS", "BDS", "GLO", "GAL" };

    while (timediff(end, start) >= 0)
    {
        time2epoch(start, epoch);
        year = static_cast<int>(epoch[0]) % 100;
        doy = time2doy(start);
        if (opt.strSys == "ALL")
        {
            for (unsigned i = 0; i != 4; ++i)
            {
                snprintf(file, 20, "cmp%03d0_%02d_%s.txt", doy, year, gnss[i].c_str());
                filesystem::path path = filesystem::path{ opt.strPath } / filesystem::path{ file };
                diff(path.str(), opt.iInterval, epoch, gnss[i]);
            }
        }
        else
        {
                snprintf(file, 19, "cmp_%03d0_%02d_%s.txt", doy, year, opt.strSys.c_str());
                filesystem::path  path = filesystem::path{ opt.strPath } / filesystem::path{ file };
                diff(path.str(), opt.iInterval, epoch, opt.strSys);
        }
        start = timeadd(start, 86400);
    }

    aveFinal(opt);
}

int main(int argc, char* argv[])
{
    Options opt{};
    parseArg(argc, argv, opt);
    process(opt);
}