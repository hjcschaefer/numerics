
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <stdlib.h>
#include <map>
#include <random>
#include <sstream>
#include <iomanip>

#include <tools/SimplePlot/SimplePlot.hpp>
#include <IMCPricing/rate/RateCurve_impl.hpp>
#include <IMCPricing/date/Time.hpp>
#include <fitlib/NelsonSiegelSvenssonModel.hpp>


// at which days to evaluate the rate curve
std::vector<int> DAYS = {0, 1, 2, 3, 4, 5, 10, 30, 60, 90, 120, 365, 500, 700};

std::map<int, std::vector<double>> INDEX;

// read in the data:
// returns a vector where each element is the rate curve evaluated at DAYS
std::vector<std::vector<double>> readdata(const char *fname) {
    std::vector<std::vector<double>> data;

    Time tval(2016, 1, 1, 0.0);
    std::vector<Time> pillars;
    Time tx = tval;
    while (tx < Time(2018, 1, 1)) {
        pillars.push_back(tx);
        tx = Time(tx.julianDayNumber()+2, 0.0);
    }

    std::ifstream fin(fname);
    if (!fin.good()) {
        std::cerr << "COULD NOT OPEN " << fname << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;

    int j=0;
    while(!fin.eof() && std::getline(fin, line)) {
        std::stringstream strstr(line);

        std::vector<double> tmp;
        while (!strstr.eof()) {
            std::string a;
            strstr >> a;
            if (a=="") break;
            tmp.push_back(std::stod(a));
        }
        NelsonSiegelSvenssonModel nss(tmp, false);
        auto rc = nss.buildCurve(tval, pillars);
        std::vector<double> tmp2;
        for(int d : DAYS) {
            Time tx(tval.julianDayNumber()+d, 0.0);
            double r = rc->forwardRate(tval, tx);
            tmp2.push_back(r);
        }
        INDEX[j] = tmp;
        j++;
        data.push_back(tmp2);
    }
    fin.close();
    return data;
}


// calculate the distance between two rate curves, returns scalar
double distance(const std::vector<double>& a, const std::vector<double>& b) {
    double d = 0.0;
    for(size_t i=0; i<a.size(); i++) {
        d += (a[i]-b[i])*(a[i]-b[i]);
    }
    return std::sqrt(d);
}


// calculate all the distances between all samples and each cluster center
// returns: vector which for each cluster contains all the distances between that center and each sample.
// so, element k of the return vector has all the distances to all samples
std::vector<std::vector<double>> distances(int clusters,
                                           const std::vector<std::vector<double>>& centers,
                                           const std::vector<std::vector<double>>& samples) {
    std::vector<std::vector<double>> dists(clusters);

    for(int i=0; i<clusters; i++) {
        auto& center = centers[i];
        int j = 0;
        for(auto e : samples) {
            auto dst = distance(center, e);
            dists[i].push_back(dst);
            j++;
        }
    }
    return dists;
}

// returns random centers which are calculated by randomly assigning samples to a cluster
// and then calculating the centroids;
std::vector<std::vector<double>> randomCenters(int clusters, const std::vector<std::vector<double>>& samples) {

    size_t n = samples[0].size(); // vector length of each sample, basically on how many days evaluated
    std::vector<std::vector<double>> centers(clusters);
    for(int i=0; i<clusters; i++) {
        centers[i] = std::vector<double>(n, 0.0);
    }

    std::vector<int> counts(clusters, 0);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> d(0, clusters-1);

    for(auto sample : samples) {
        int whichCenter = d(gen); // randomly assign to a center
        counts[whichCenter]++;

        // already part of the centroid calculation
        for(size_t i=0; i<n; i++) {
            centers[whichCenter][i] += sample[i];
        }
    }

    // centroid calculation
    for(int i=0; i<clusters; i++) {
        for(size_t j=0; j<n; j++) {
            centers[i][j] /= counts[i];
        }
    }
    return centers;
}

int biggestCluster(int clusters, const std::map<int,std::vector<int>>& closest) {
    int n = closest.at(0).size();
    int j = 0;
    for(int k=1; k<clusters; k++) {
        if (closest.at(k).size() > n) {
            n = closest.at(k).size();
            j = k;
        }
    }
    return j;
}

// given all the centers we check for each sample how close it is to each center. We then
// assign it to the closest center
std::map<int, std::vector<int>> mapToCenter(int clusters,
                                            const std::vector<std::vector<double>>& centers,
                                            const std::vector<std::vector<double>>& samples) {
    // distances to all centers
    auto dists = distances(clusters, centers, samples);
    // now associate them to the closest one
    std::map<int, std::vector<int>> closest;
    for(int i=0; i<clusters; i++) {
        closest[i] = std::vector<int>();
    }

    for(int k=0; k<samples.size(); k++) { // loop over all samples
        int min = 0;
        double d = dists[0][k]; // distance to center 0 for this sample
        for (int j=1; j<centers.size(); j++) { // loop over the centers
            if (dists[j][k] < d) { // find smallest distance
                min = j;
                d = dists[j][k];
            }
        }
        closest[min].push_back(k);
    }
    // check if we have an empty cluster and deal with it
    for (int j=0; j<centers.size(); j++) { // loop over the centers
        if (closest[j].size() == 0) {
            std::cerr << "Empty Cluster!!!" << std::endl;
            // now we look at the biggest center and in there find the most outlying member and use this as a new center
            auto biggest = biggestCluster(clusters, closest);
            // go through all the distances in this cluster and find the biggest one
            double d = dists[biggest][0];
            int idx = 0;
            for (int k=1; k<dists[biggest].size(); k++) {
                if (dists[biggest][k] > d) {
                    d = dists[biggest][k];
                    idx = k;
                }
            }
            // so element idx gets its own cluster and needs to be moved:
            closest[j] = {idx};
            std::vector<int> tmp;
            for (int k : closest[biggest]) {
                if (k != idx) {
                    tmp.push_back(k);
                }
            }
            closest[biggest] = tmp;
            std::cerr << "Moved " << idx << " From cluster " << biggest << " to " << j << std::endl;
        }
    }

    return closest;
}

// given the samples and the map which associates each sample to a cluster we can calculate the centroids
// i.e. the new centers
std::vector<std::vector<double>> centroids(int clusters,
                                           const std::vector<std::vector<double>>& samples,
                                           const std::map<int, std::vector<int>>& closest) {

    size_t ndim = samples[0].size();
    std::vector<std::vector<double>> centers(clusters, std::vector<double>(ndim, 0.0));
    // calculate new centroids as new centers;
    for (int center=0; center<clusters; center++) {
        for(size_t dim=0; dim<ndim; dim++) {
            double sum = 0.0;
            for(size_t k=0; k<closest.at(center).size(); k++) {
                sum += samples[closest.at(center)[k]][dim];
            }
            centers[center][dim] = sum/closest.at(center).size();
        }
    }
    return centers;
}

// --------------- Main Algorithm ----------------------------
std::vector<std::vector<double>> cluster(const std::vector<std::vector<double>>& samples, int clusters, int steps) {
    auto centers = randomCenters(clusters, samples);

    double err = 10.0;
    for(int i=0; i<steps && err > 1e-5; i++) {

        auto closest = mapToCenter(clusters, centers, samples);
        std::cerr << "-----------------" << std::endl;
        for (int center=0; center<clusters; center++) {
            std::cerr << "\t" << center << "\t:\t" << closest[center].size() << std::endl;
        }
        auto tmp = centroids(clusters, samples, closest);

        // have we converged?
        err = 0.0;
        for(size_t j=0; j<tmp.size(); j++) {
            err += distance(tmp[j], centers[j]);
        }
        std::cerr << "Error : " << err << std::endl;
        centers = tmp;
    }

    // reassign
    auto closest = mapToCenter(clusters, centers, samples);
    return centers;
}

std::vector<double> variances(int clusters, const std::vector<std::vector<double>>& centers,
                              const std::vector<std::vector<double>>& samples,
                              const std::map<int, std::vector<int>>& closest) {

    std::vector<double> variance(clusters, 0.0);
    for (int center=0; center<clusters; center++) {
        auto& cluster = closest.at(center);
        for(auto idx: cluster) {
            auto& e = samples[idx];
            variance[center] += distance(e, centers.at(center));
        }
        variance[center]/= cluster.size();
    }
    return variance;
}


int main(int argc, char* argv[]) {

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <Number of Cluster> <Data File> <steps>" << std::endl;
        exit(EXIT_FAILURE);
    }

    int clusters = std::stoi(argv[1]);
    int steps = std::stoi(argv[3]);

    const char *fname = argv[2];
    auto samples = readdata(fname);
    std::cerr << "Read " << samples.size() << " Samples." << std::endl;

    auto centers = cluster(samples, clusters, steps);
    // reassign
    auto closest = mapToCenter(clusters, centers, samples);
    // calculate the variance of the clusters
    std::vector<double> variance = variances(clusters, centers, samples, closest);
    std::cerr << "Variances: " << std::endl;
    double sum = 0.0;
    for (int center=0; center<clusters; center++) {
        std::cerr << "\t" << center << "\t" << variance[center] << std::endl;;
        sum += variance[center];
    }
    std::cerr << "*** Overall " << sum << std::endl;

    Time tval(2016, 1, 1, 0.0);
    std::vector<Time> pillars;
    Time tx = tval;
    while (tx < Time(2018, 1, 1)) {
        pillars.push_back(tx);
        tx = Time(tx.julianDayNumber()+2, 0.0);
    }

    // plot cluster:
    SimplePlot pltx(".", "all.png", samples.size(), SimplePlot::PNG);
    SimplePlot plty(".", "nss.png", clusters, SimplePlot::PNG);
    plty.Lines(false);
    plty.Points(true);
    int gx=0;
    int days = 700;
    std::map<size_t,std::string> legend;
    for (int center=0; center<clusters; center++) {
        auto& cluster = closest[center];
        std::ostringstream ox;
        ox << "Var = " << std::setw(12) << variance[center] << " " << cluster.size();
        legend[gx] = ox.str();
        std::ostringstream os;
        os << "cluster_" << center;
        SimplePlot plt(".", os.str().c_str(), cluster.size(), SimplePlot::PNG);
        int g = 0;
        for(auto idx: cluster) {
            // plot nss parameters
            auto nss = INDEX[idx];
            plty.setColor(center, center);
            plty.AddToGraph(center, 0, 100*nss[0]);
            plty.AddToGraph(center, 1, 100*nss[1]);
            plty.AddToGraph(center, 2, 100*nss[2]);
            plty.AddToGraph(center, 3, nss[3]);
            pltx.setColor(gx, center);
            auto& e = samples[idx];
            int i=0;
            for(int d : DAYS) {
            //for(int d=0; d<days; d+=10) {
                plt.AddToGraph(g, d, e[i]*100);
                pltx.AddToGraph(gx, d, e[i]*100);
                i++;
            }
            g++;
            gx++;
        }
        plt.Plot();
        pltx.addText(legend);
        pltx.Plot();
        plty.Plot();
    }
}
