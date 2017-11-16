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
void generateEdgarGilbertModelByDensity(int nodesCount, std::vector<double> densityList);
void generateEdgarGilbertModelByNodesCount(std::vector<int> nodesCountList, double density);
void generateStochasticBlockModel(std::vector<int> communitiesSize, double withinCommunitiesDensity);
void deleteNodes(int desiredNodesCount, int &currentNodesCount, int &currentEdgesCount,
                 std::map<string, std::vector<string>> &graph);
void deleteEdgesToDensity(double density, int nodesCount, int &currentEdgesCount,
                          std::map<string, std::vector<string>> &graph);
void writeGraph(int nodesCount, int edgesCount, const std::map<string, std::vector<string>> &graph);
unsigned long getRandomNumber(unsigned long max, unsigned long min);
void setupRandomSeed();
template<typename T1, typename T2>
struct higher_second {
    typedef std::pair<T1, T2> type;

    bool operator()(type const &a, type const &b) const {
        return a.second > b.second;
    }
};

int main() {
    setupRandomSeed();
    // Generate Edgar Gilbert Model graph
//    std::vector<int> nodesCountList {170};
//    double density = 0.093039;
//    generateEdgarGilbertModelByNodesCount(nodesCountList, density);

//    int nodesCount = 1000;
//    std::vector<double> densityList{0.9, 0.7, 0.5, 0.3, 0.1, 0.05};
//    generateEdgarGilbertModelByDensity(nodesCount, densityList);

    // Generate Stochastic Block Model graph
//    std::vector<int> communitiesSize {10, 12, 15, 15, 18, 20, 30, 50};
//    double communitiesDensity = 0.5;
//    generateStochasticBlockModel(communitiesSize, communitiesDensity);

    // Execute PageRank algorithm
    const int TEST_CASE = 28;
    const double SCALING_FACTOR = 0.85;
    const int MAX_ITERATION = 100;
    const double TOLERANCE_VALUE = 0.000001;
    execPageRankAlgo(TEST_CASE, MAX_ITERATION, SCALING_FACTOR, TOLERANCE_VALUE);

    return 0;
}

void execPageRankAlgo(int TEST_CASE, int MAX_ITERATION, double SCALING_FACTOR, double TOLERANCE_VALUE) {
    // Read input
    std::map<string, std::vector<string>> pointFrom;
    std::map<string, std::vector<string>> pointTo;
    clock_t clock_1 = clock();
    readURLs("/Users/midnightblur/Documents/workspace/CLionProjects/COMP6651_Assignment03/test cases/test case " +
                     std::to_string(TEST_CASE) + "/links.txt", pointFrom, pointTo);
    clock_t clock_2 = clock();
    double dict_diff((double) clock_2 - (double) clock_1);
//    cout << "Reading input = " << dict_diff / (CLOCKS_PER_SEC / 1000) << " ms" << endl;
    cout << dict_diff / (CLOCKS_PER_SEC / 1000) << endl;


    // Init page rank to 1 for all URLs
    std::map<string, double> pageRank;
    for (auto &url : pointTo) {
        pageRank.insert(std::make_pair(url.first, 1));
    }

    // Calculate PageRank for all URLs
    clock_1 = clock();
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
            clock_2 = clock();
            dict_diff = (double) clock_2 - (double) clock_1;
//            cout << "Converge time = " << dict_diff / (CLOCKS_PER_SEC / 1000) << " ms" << endl;
            cout << dict_diff / (CLOCKS_PER_SEC / 1000) << endl;
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
            destination.push_back(url2);
            pointFrom.insert(std::make_pair(url1, destination));
        } else {
            std::vector<string> destination = fromURL1->second;
            destination.push_back(url2);
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
            destination.push_back(url1);
            pointTo.insert(std::make_pair(url2, destination));
        } else {
            std::vector<string> destination = toURL2->second;
            destination.push_back(url1);
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
//        cout << pair.first << ": " << pair.second << endl;
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

void generateEdgarGilbertModelByDensity(int nodesCount, std::vector<double> densityList) {
    int currentEdgesCount = nodesCount * (nodesCount - 1);

    // Generate nodes
    std::map<string, std::vector<string>> graph;
    for (int i = 1; i <= nodesCount; i++) {
        std::vector<string> destinationURLs;
        graph.insert(std::make_pair(std::to_string(i), destinationURLs));
    }

    // Make a fully connected graph
    for (auto &sourceNode : graph) {
        for (auto &targetNode : graph) {
            if (targetNode != sourceNode) {
                sourceNode.second.push_back(targetNode.first);
            }
        }
    }

    // Delete edges to reach desired density
    for (auto &density : densityList) {
        cout << endl << "New Erdos Renyi Model" << endl;
        deleteEdgesToDensity(density, nodesCount, currentEdgesCount, graph);
        writeGraph(nodesCount, currentEdgesCount, graph);
    }
}

void generateEdgarGilbertModelByNodesCount(std::vector<int> nodesCountList, double density) {
    int currentNodesCount = nodesCountList.front();
    int currentEdgesCount = currentNodesCount * (currentNodesCount - 1);

    // Generate nodes
    std::map<string, std::vector<string>> graph;
    for (int i = 1; i <= currentNodesCount; i++) {
        std::vector<string> destinationURLs;
        graph.insert(std::make_pair(std::to_string(i), destinationURLs));
    }

    // Make a fully connected graph
    for (auto &sourceNode : graph) {
        for (auto &targetNode : graph) {
            if (targetNode != sourceNode) {
                sourceNode.second.push_back(targetNode.first);
            }
        }
    }

    // Delete edges to reach desired density
    deleteEdgesToDensity(density, currentNodesCount, currentEdgesCount, graph);

    // Delete nodes to reach desired number of nodes
    for (auto &nodeCount : nodesCountList) {
        cout << endl << "New Erdos Renyi Model" << endl;
        deleteNodes(nodeCount, currentNodesCount, currentEdgesCount, graph);
        writeGraph(currentNodesCount, currentEdgesCount, graph);
    }
}

void generateStochasticBlockModel(std::vector<int> communitiesSize, double withinCommunitiesDensity) {
    cout << endl << "New Stochastic Block Model" << endl;
    // Init communities and their members
    int nodesCount = 0;
    for (auto size : communitiesSize) {
        nodesCount += size;
    }

    int edgesCount = 0;
    std::map<string, std::vector<string>> communitiesMap;
    std::map<string, std::vector<string>> graph;
    for (int i = 1; i <= communitiesSize.size(); i++) {
        int size = communitiesSize.at(static_cast<unsigned long>(i - 1));
        string communityID = std::to_string(i);
        std::vector<string> communityMemberVector;
        for (int j = 0; j < size; j++) {
            string nodeID = communityID + "-" + std::to_string(j);
            communityMemberVector.push_back(nodeID);

            std::vector<string> destinationNodes;
            graph.insert(std::make_pair(nodeID, destinationNodes));
        }
        communitiesMap.insert(std::make_pair(communityID, communityMemberVector));
    }

    // Make all communities fully connected
    for (auto &community : communitiesMap) {
        for (auto &sourceNode : community.second) {
            for (auto &targetNode : community.second) {
                if (targetNode != sourceNode) {
                    graph.find(sourceNode)->second.push_back(targetNode);
                }
            }
        }
    }

    // Make all communities reach desired density by removing arbitrary edges within communities
    for (auto &community : communitiesMap) {
        int size = static_cast<int>(community.second.size());
        int currentCommunityEdgesCount = size * (size - 1);
        auto desiredCommunityEdgesCount = static_cast<int>(currentCommunityEdgesCount * withinCommunitiesDensity);
        while (currentCommunityEdgesCount > desiredCommunityEdgesCount) {
            // Pick a random sourceNode in the community
            auto sourceNodeIter = community.second.begin();
            std::advance(sourceNodeIter, getRandomNumber(community.second.size() - 1, 0));

            // Pick a random targetNode and remove the edge
            std::vector<string> &targetNodes = graph.find(*sourceNodeIter)->second;
            if (targetNodes.size() > 1) { // Only delete the edge if there're more than 1 edge coming from sourceNode
                auto targetNodeIter = targetNodes.begin();
                std::advance(targetNodeIter, getRandomNumber(targetNodes.size() - 1, 0));
                targetNodes.erase(targetNodeIter);

                currentCommunityEdgesCount--;
            }
        }
        edgesCount += desiredCommunityEdgesCount;
    }

    // Create some arbitrary edges between communities
    for (auto &sourceCommunity : communitiesMap) {
        for (auto &targetCommunity : communitiesMap) {
            if (sourceCommunity != targetCommunity) {
                auto noOfNewEdges = static_cast<int>(getRandomNumber(sourceCommunity.second.size() / 2, 1));
                for (int i = 0; i < noOfNewEdges; i++) {
                    // Pick a random sourceNode from sourceCommunity
                    auto sourceNodeIter = sourceCommunity.second.begin();
                    std::advance(sourceNodeIter, getRandomNumber(sourceCommunity.second.size() - 1, 0));

                    // Pick a random targetNode from targetCommunity
                    auto targetNodeIter = targetCommunity.second.begin();
                    std::advance(targetNodeIter, getRandomNumber(targetCommunity.second.size() - 1, 0));

                    // Create new edge between them
                    auto &destinationNodes = graph.find(*sourceNodeIter)->second;
                    if (std::find(destinationNodes.begin(), destinationNodes.end(), *targetNodeIter) == destinationNodes.end()) {
                        destinationNodes.push_back(*targetNodeIter);
                        edgesCount++;
                    }
                }
            }
        }
    }

    // Write the graph to file
    writeGraph(nodesCount, edgesCount, graph);
}

void deleteNodes(int desiredNodesCount, int &currentNodesCount, int &currentEdgesCount,
                 std::map<string, std::vector<string>> &graph) {
    while (currentNodesCount > desiredNodesCount) {
        // Pick a random node
        auto nodeIter = graph.begin();
        std::advance(nodeIter, getRandomNumber(graph.size() - 1, 0));

        // Delete the node and all edges coming from it
        string node = nodeIter->first;
        currentEdgesCount -= nodeIter->second.size();
        graph.erase(nodeIter);

        // Delete all edges coming to it
        for (auto &sourceNode : graph) {
            for (auto targetIter = sourceNode.second.begin(); targetIter != sourceNode.second.end(); targetIter++) {
                if (*targetIter == node) {
                    sourceNode.second.erase(targetIter);
                    currentEdgesCount--;
                    break;
                }
            }
        }
        currentNodesCount--;
    }
}

void deleteEdgesToDensity(double density, int nodesCount, int &currentEdgesCount,
                          std::map<string, std::vector<string>> &graph) {
    int maxEdgesCount = nodesCount * (nodesCount - 1);
    auto targetEdgesCount = static_cast<int>(maxEdgesCount * density);
    while (currentEdgesCount > targetEdgesCount) {
        // Pick a random sourceNode
        auto sourceNodeIter = graph.begin();
        std::advance(sourceNodeIter, getRandomNumber(graph.size() - 1, 0));

        if (sourceNodeIter->second.size() > 1) { // keep at least 1 edge so the graph is still a connected one
            // Pick a random targetNode to remove the edge
            auto targetNodeIter = sourceNodeIter->second.begin();
            std::advance(targetNodeIter, getRandomNumber(sourceNodeIter->second.size() - 1, 0));

            sourceNodeIter->second.erase(targetNodeIter);
            currentEdgesCount--;
        }
    }
}

void writeGraph(int nodesCount, int edgesCount, const std::map<string, std::vector<string>> &graph) {
    double density = (double) edgesCount / (nodesCount * (nodesCount - 1));
    std::ofstream outputFile;
    outputFile.open(
            "/Users/midnightblur/Documents/workspace/CLionProjects/COMP6651_Assignment03/test cases/generatedGraphs/"
                    "n" + std::to_string(nodesCount) + "-e" + std::to_string(edgesCount) + "-d" + std::to_string(
                    density) + ".txt");

    for (auto &sourceNode : graph) {
        for (auto &targetNode : sourceNode.second) {
            outputFile << sourceNode.first << ", " << targetNode << endl;
        }
    }

    outputFile.close();

    cout << "No of Nodes = " << nodesCount << endl;
    cout << "No of Edges = " << edgesCount << endl;
    cout << "Graph density = " << density << endl;
}

unsigned long getRandomNumber(unsigned long max, unsigned long min) {
    return rand() % (max - min + 1) + min;
}

void setupRandomSeed() {
    struct timespec ts{};
    clock_gettime(CLOCK_MONOTONIC, &ts);
    srand(static_cast<unsigned int>((time_t) ts.tv_nsec));
}