#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <iomanip>
#include <cmath>

using std::string;
using std::endl;
using std::cout;

void execPageRankAlgo(int TEST_CASE, int MAX_ITERATION, double SCALING_FACTOR, double TOLERANCE_VALUE);

void readURLs(const string &filePath, std::map<string, std::vector<string>> &pointFrom,
              std::map<string, std::vector<string>> &pointTo);

void removeR(string &line);

void writeOutput(const std::map<string, double> &pageRank, int TEST_CASE);

void writeReadme(int TEST_CASE, double scaling_factor, int max_iteration);

void generateGraph(int TEST_CASE, int nodesCount, float density);

unsigned long getRandomNumber(unsigned long max, unsigned long min);

template<typename T1, typename T2>
struct higher_second {
    typedef std::pair<T1, T2> type;

    bool operator()(type const &a, type const &b) const {
        return a.second > b.second;
    }
};

int main() {
    // Configuration
    const int TEST_CASE = 3;
    const double SCALING_FACTOR = 0.85;
    const int MAX_ITERATION = 30;
    const double TOLERANCE_VALUE = 0.000001;

    generateGraph(TEST_CASE, 100, 0.05);
//    execPageRankAlgo(TEST_CASE, MAX_ITERATION, SCALING_FACTOR, TOLERANCE_VALUE);

    return 0;
}

void execPageRankAlgo(int TEST_CASE, int MAX_ITERATION, double SCALING_FACTOR, double TOLERANCE_VALUE) {
    // Read input
    std::map<string, std::vector<string>> pointFrom;
    std::map<string, std::vector<string>> pointTo;
    readURLs("/Users/midnightblur/Documents/workspace/CLionProjects/COMP6651_Assignment03/test cases/test case " +
                     std::to_string(TEST_CASE) + "/links.txt", pointFrom, pointTo);

    // Init page rank to 1 for all URLs
    std::map<string, double> pageRank;
    for (auto &url : pointTo) {
        pageRank.insert(std::make_pair(url.first, 1));
    }

    // Calculate PageRank for all URLs
    clock_t clock_1 = clock();
    for (int i = 0; i < MAX_ITERATION; i++) {
        std::map<string, double> tmpPageRank;
        for (auto currentURL = pointTo.begin(); currentURL != pointTo.end(); currentURL++) {
            double currentPageRank = 1 - SCALING_FACTOR;
            for (auto &source : currentURL->second) {
                double pageRankDestination = pageRank.find(source)->second;
                size_t outboundLinksDestination = pointFrom.find(source)->second.size();
                currentPageRank += SCALING_FACTOR * pageRankDestination / outboundLinksDestination;
            }
            tmpPageRank.insert(std::make_pair(currentURL->first, currentPageRank));
        }

        // Check for convergence
        bool isConverged = true;
        for (auto &pair : pageRank) {
            if (fabs(pair.second - tmpPageRank.find(pair.first)->second) > TOLERANCE_VALUE) {
                isConverged = false;
                break;
            }
        }
        if (isConverged) {
            clock_t clock_2 = clock();
            double dict_diff((double) clock_2 - (double) clock_1);
            cout << "Converge time = " << dict_diff / (CLOCKS_PER_SEC / 1000) << " ms" << endl;
            cout << "Iteration no " << i + 1 << endl;
            break;
        } else {
            pageRank = tmpPageRank;
        }
    }

    // Sort the result based on PageRank value

    writeOutput(pageRank, TEST_CASE);
    writeReadme(TEST_CASE, SCALING_FACTOR, MAX_ITERATION);
}

void readURLs(const string &filePath, std::map<string, std::vector<string>> &pointFrom,
              std::map<string, std::vector<string>> &pointTo) {
    std::ifstream file(filePath);
    string line;

    std::hash<std::string> str_hash;
    while (std::getline(file, line)) {
        string delimiter = ", ";
        string url1 = line.substr(0, line.find(delimiter));
        line.erase(0, line.find(delimiter) + delimiter.length());
        string url2 = line;
        removeR(url2);

        // region Update pointFrom map
        auto fromURL1 = pointFrom.find(url1);
        auto fromURL2 = pointFrom.find(url2);
        if (fromURL1 == pointFrom.end()) {
            std::vector<string> destination;
            destination.insert(destination.end(), url2);
            pointFrom.insert(std::make_pair(url1, destination));
        } else {
            std::vector<string> destination = fromURL1->second;
            destination.insert(destination.end(), url2);
            pointFrom.erase(fromURL1);
            pointFrom.insert(std::make_pair(url1, destination));
        }

        if (fromURL2 == pointFrom.end()) {
            std::vector<string> destination;
            pointFrom.insert(std::make_pair(url2, destination));
        }
        // endregion

        // region Update pointTo map
        auto toURL1 = pointTo.find(url1);
        auto toURL2 = pointTo.find(url2);
        if (toURL2 == pointTo.end()) {
            std::vector<string> destination;
            destination.insert(destination.end(), url1);
            pointTo.insert(std::make_pair(url2, destination));
        } else {
            std::vector<string> destination = toURL2->second;
            destination.insert(destination.end(), url1);
            pointTo.erase(toURL2);
            pointTo.insert(std::make_pair(url2, destination));
        }

        if (toURL1 == pointTo.end()) {
            std::vector<string> destination;
            pointTo.insert(std::make_pair(url1, destination));
        }
        // endregion
    }

    file.close();

    // region Display input
//    cout << "Point From" << endl;
//    for (auto &urlIT : pointFrom) {
//        cout << urlIT.first << ": ";
//        for (auto &it : urlIT.second)
//            cout << it << " ";
//        cout << endl;
//    }
//
//    cout << "Point To" << endl;
//    for (auto &urlIT : pointTo) {
//        cout << urlIT.first << ": ";
//        for (auto &it : urlIT.second)
//            cout << it << " ";
//        cout << endl;
//    }
    // endregion
}

void removeR(string &line) {
    if (!line.empty() && line.at(line.size() - 1) == '\r')
        line.erase(line.size() - 1);
}

void writeOutput(const std::map<string, double> &pageRank, const int TEST_CASE) {
    std::vector<std::pair<string, double>> resultVector(pageRank.begin(), pageRank.end());
    sort(resultVector.begin(), resultVector.end(), higher_second<string, double>());

    std::ofstream outputFile;
    outputFile.open(
            "/Users/midnightblur/Documents/workspace/CLionProjects/COMP6651_Assignment03/test cases/test case " +
                    std::to_string(TEST_CASE) + "/Output.txt");
    for (auto &pair : resultVector) {
        cout << pair.first << ": " << pair.second << endl;
        outputFile << pair.first << ", " << pair.second << endl;
    }
    outputFile.close();
}

void writeReadme(const int TEST_CASE, const double scaling_factor, const int max_iteration) {
    std::ofstream outputFile;
    outputFile.open(
            "/Users/midnightblur/Documents/workspace/CLionProjects/COMP6651_Assignment03/test cases/test case " +
                    std::to_string(TEST_CASE) + "/readme.txt");
    outputFile << std::fixed;
    outputFile << "For this test case parametric values used are:" << endl;
    outputFile << "scaling factor=" << std::setprecision(2) << scaling_factor << endl;
    outputFile << "maximum iterations=" << max_iteration << endl;
    outputFile << endl << "Note:" << endl;
    outputFile << "For experiements to show the effect of the size of the input and the scaling factor on the "
            "convergence time, check for convergence with tolerance value=1.0e-6. "
            "If your graph is a large graph set maximum iterations to 100.";
    outputFile.close();
}

void generateGraph(int TEST_CASE, int nodesCount, float density) {
    int maxEdgesCount = nodesCount * (nodesCount - 1);
    auto edgesCount = static_cast<int>(maxEdgesCount * density);
    cout << "No of Nodes = " << nodesCount << endl;
    cout << "No of Edges = " << edgesCount << endl;
    cout << "Graph density = " << density << endl;

    // Generate nodes
    std::map<string, std::vector<string>> graph;
    for (int i = 1; i <= nodesCount; i++) {
        std::vector<string> destinationURLs;
        graph.insert(std::make_pair(std::to_string(i), destinationURLs));
    }

    // Make a fully connected graph
    for (auto &sourceNode : graph) {
        std::vector<string> desitnationURLs = sourceNode.second;
        for (auto &targetNode : graph) {
            if (targetNode != sourceNode) {
                desitnationURLs.insert(desitnationURLs.end(), targetNode.first);
            }
        }
        sourceNode.second = desitnationURLs;
    }

    // Delete edges to reach desired density
    while (maxEdgesCount > edgesCount) {
        // Pick a random sourceNode
        auto sourceNodeIter = graph.begin();
        std::advance(sourceNodeIter, getRandomNumber(graph.size() - 1, 0));

        if (sourceNodeIter->second.size() > 1) { // keep at least 1 edge so the graph is still a connected one
            // Pick a random targetNode to remove the edge
            auto targetNodeIter = sourceNodeIter->second.begin();
            std::advance(targetNodeIter, getRandomNumber(sourceNodeIter->second.size() - 1, 0));

            sourceNodeIter->second.erase(targetNodeIter);
            maxEdgesCount--;
        }
    }

    // Write graph to file
    std::ofstream outputFile;
    outputFile.open(
            "/Users/midnightblur/Documents/workspace/CLionProjects/COMP6651_Assignment03/test cases/test case " +
                    std::to_string(TEST_CASE) + "/links.txt");

    for (auto &sourceNode : graph) {
        for (auto &targetNode : sourceNode.second) {
            outputFile << sourceNode.first << ", " << targetNode << endl;
        }
    }

    outputFile.close();
}

unsigned long getRandomNumber(unsigned long max, unsigned long min) {
    struct timespec ts{};
    clock_gettime(CLOCK_MONOTONIC, &ts);
    srand(static_cast<unsigned int>((time_t) ts.tv_nsec));
    return rand() % (max - min + 1) + min;
}