
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <stdlib.h>
#include <map>
#include <random>
#include <sstream>

#include <tools/SimplePlot/SimplePlot.hpp>
#include <IMCPricing/rate/RateCurve_impl.hpp>
#include <IMCPricing/date/Time.hpp>
#include <fitlib/NelsonSiegelSvenssonModel.hpp>

int K;

std::vector<double> readnss(char* fname) {
    std::ifstream fin(fname);
    std::vector<double> paras;
    do {
        std::string x;
        fin >> x;
        if (!x.empty()) {
            paras.push_back(atof(x.c_str()));
        }
    } while(fin);
    fin.close();
    if (paras.size() < 4) {
        return std::vector<double>();
    }
    if (paras.size() > 4) {
        std::vector<double> bla;
        bla.push_back(paras[0]);
        bla.push_back(paras[1]);
        bla.push_back(paras[2]);
        bla.push_back(paras[3]);
        return bla;
    }
    return paras;
}

std::vector<std::vector<double>> readdata(const char *fname) {
    std::vector<std::vector<double>> data;

    std::ifstream fin(fname);
    std::string line;

    while(!fin.eof() && std::getline(fin, line)) {
        std::stringstream strstr(line);

        std::vector<double> tmp;
        while (!strstr.eof()) {
            std::string a;
            strstr >> a;
            if (a=="") break;
            tmp.push_back(std::stod(a));
        }
        data.push_back(tmp);
    }
    fin.close();
    return data;
}


double distance(const std::vector<double>& a, const std::vector<double>& b) {
    double d = 0.0;
    for(size_t i=0; i<a.size(); i++) {
        d += (a[i]-b[i])*(a[i]-b[i]);
    }
    return std::sqrt(d);
}


std::vector<std::vector<double>> distances(const std::vector<std::vector<double>>& centers,
                                           const std::vector<std::vector<double>>& samples) {
    std::vector<std::vector<double>> dists(K);

    for(int i=0; i<K; i++) {
        auto& center = centers[i];
        for(auto e : samples) {
            dists[i].push_back(distance(center, e));
        }
    }
    return dists;
}

// random samples
std::vector<std::vector<double>> random2d() {
    std::random_device rd;
    std::mt19937 gen(rd());
    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    double stddev = 0.9;
    std::normal_distribution<double> d1(5, stddev);
    std::normal_distribution<double> d2(1, stddev);
    std::normal_distribution<double> d3(2, stddev);

    std::vector<std::vector<double>> samples;
    for(int i=0; i<1000; i++) {
        double a1 = d1(gen);
        double b1 = d1(gen);
        samples.push_back({a1, b1});
        double a2 = d1(gen);
        double b2 = d2(gen);
        samples.push_back({a2, b2});
        double a3 = d2(gen);
        double b3 = d3(gen);
        samples.push_back({a3, b3});
    }

    SimplePlot plt(".", "samples.png", 1, SimplePlot::PNG);
    plt.Lines(false);
    plt.Points(true);

    for(auto e: samples) {
        plt.Add(e[0], e[1]);
    }
    plt.Plot();

    return samples;
}

int main(int argc, char* argv[]) {

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <Number of Cluster> <Data File> <steps>" << std::endl;
        exit(EXIT_FAILURE);
    }

    K = std::stoi(argv[1]);
    int steps = std::stoi(argv[3]);

    const char *fname = argv[2];
    auto samples = readdata(fname);

    std::vector<std::vector<double>> centers(K);
    //auto samples = random2d();

    // select random centers, just the first K
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> d(0, samples.size());
    for(int i=0; i<K; i++) {
        //centers[i] = samples[i];
        // really randomize
        int idx = d(gen);
        std::cerr << "Center " << i << " : " << idx << std::endl;
        centers[i] = samples[idx];
    }

    std::cerr << "===================" << std::endl;
    for(int i=0; i<steps; i++) {

        for (int center=0; center<K; center++) {
            std::cerr << centers[center][0] << "\t" <<  centers[center][1] << "\t" << centers[center][2] << "\t" << centers[center][3]  << std::endl;
        }
        auto dists = distances(centers, samples);
        // now associate them to the closest one
        std::map<int, std::vector<int>> closest;
        for(int k=0; k<samples.size(); k++) {
            int min = 0;
            double d = dists[0][k];
            for (int j=1; j<K; j++) {
                if (dists[j][k] < d) {
                    min = j;
                    d = dists[j][k];
                }
            }
            closest[min].push_back(k);
        }
        // calculate new centroids as new centers;
        for (int center=0; center<K; center++) {
            for(size_t dim=0; dim<centers[center].size(); dim++) {
                double sum = 0.0;
                for(size_t k=0; k<closest[center].size(); k++) {
                    sum += samples[closest[center][k]][dim];
                }
                centers[center][dim] = sum/closest[center].size();
            }
        }

        std::cerr << "------------------" << std::endl;
    }

    // reassign
    auto dists = distances(centers, samples);
    // now associate them to the closest one
    std::map<int, std::vector<int>> closest;
    for(int k=0; k<samples.size(); k++) {
        int min = 0;
        double d = dists[0][k];
        for (int j=1; j<K; j++) {
            if (dists[j][k] < d) {
                min = j;
                d = dists[j][k];
            }
        }
        closest[min].push_back(k);
    }

    std::cerr << "+++++++++++++++++++++++++++++++++++++++" << std::endl;
    for (int center=0; center<K; center++) {
        std::cerr << centers[center][0] << "\t" <<  centers[center][1] << "\t" << centers[center][2] << "\t" << centers[center][3]  << std::endl;
    }

    Time tval(2016, 1, 1, 0.0);
    std::vector<Time> pillars;
    Time tx = tval;
    while (tx < Time(2018, 1, 1)) {
        pillars.push_back(tx);
        tx = Time(tx.julianDayNumber()+2, 0.0);
    }

    // plot cluster:
    SimplePlot pltx(".", "all.png", samples.size(), SimplePlot::PNG);
    int gx=0;
    for (int center=0; center<K; center++) {
        auto& cluster = closest[center];
        std::cerr << "\n++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::ostringstream os;
        os << "cluster_" << center;
        SimplePlot plt(".", os.str().c_str(), cluster.size(), SimplePlot::PNG);
        int days = 700;
        int g = 0;
        for(auto idx: cluster) {
            auto& e = samples[idx];
            NelsonSiegelSvenssonModel nss(e, false);
            auto rc = nss.buildCurve(tval, pillars);
            std::cerr << e[0] << "\t" << e[1] << "\t" << e[2] << "\t" << e[3] << std::endl;
            for(int d=0; d<days; d+=10) {
                Time tx(tval.julianDayNumber()+d, 0.0);
                double r = rc->forwardRate(tval, tx);
                plt.AddToGraph(g, d, r*100);
                pltx.AddToGraph(gx, d, r*100);
            }
            g++;
            gx++;
        }
        plt.Plot();
        pltx.Plot();
    }
}
