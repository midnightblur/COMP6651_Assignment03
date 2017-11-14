#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <iomanip>

using std::string;
using std::endl;
using std::cout;

void readURLs(const string &filePath, std::map<string, std::vector<string>> &pointFrom,
              std::map<string, std::vector<string>> &pointTo);
void removeR(string &line);
void writeOutput(const std::map<string, float> &pageRank, const int TEST_CASE);
void writeReadme(const int TEST_CASE, const float scaling_factor, const int max_iteration);

template<typename T1, typename T2>
struct higher_second {
    typedef std::pair<T1, T2> type;

    bool operator()(type const &a, type const &b) const {
        return a.second > b.second;
    }
};

int main() {
    // Configuration
    int test_case = 3;
    float scaling_factor = 0.5;
    int max_iteration = 12;

    // Read input
    std::map<string, std::vector<string>> pointFrom;
    std::map<string, std::vector<string>> pointTo;
    readURLs("/Users/midnightblur/Documents/workspace/CLionProjects/COMP6651_Assignment03/test cases/test case " +
                     std::to_string(test_case) + "/links.txt", pointFrom, pointTo);

    // Init page rank to 1 for all URLs
    std::map<string, float> pageRank;
    for (auto &url : pointTo) {
        pageRank.insert(std::make_pair(url.first, 1));
    }

    // Calculate PageRank for all URLs
    for (int i = 0; i < max_iteration; i++) {
        std::map<string, float> tmpPageRank;
        for (auto currentURL = pointTo.begin(); currentURL != pointTo.end(); currentURL++) {
            float currentPageRank = 1 - scaling_factor;
            for (auto &source : currentURL->second) {
                float pageRankDestination = pageRank.find(source)->second;
                size_t outboundLinksDestination = pointFrom.find(source)->second.size();
                currentPageRank += scaling_factor * pageRankDestination / outboundLinksDestination;
            }
            tmpPageRank.insert(std::make_pair(currentURL->first, currentPageRank));
        }
        pageRank = tmpPageRank;
    }

    // Sort the result based on PageRank value
    std::vector<std::pair<string, float>> resultVector(pageRank.begin(), pageRank.end());
    sort(resultVector.begin(), resultVector.end(), higher_second<string, float>());
    for (auto &pair : resultVector) {
        cout << pair.first << ": " << pair.second << endl;
    }
    writeOutput(pageRank, test_case);
    writeReadme(test_case, scaling_factor, max_iteration);

    return 0;
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

void writeOutput(const std::map<string, float> &pageRank, const int TEST_CASE) {
    std::ofstream outputFile;
    outputFile.open(
            "/Users/midnightblur/Documents/workspace/CLionProjects/COMP6651_Assignment03/test cases/test case " +
                    std::to_string(TEST_CASE) + "/Output.txt");
    for (auto &pair : pageRank) {
        outputFile << pair.first << ", " << pair.second << endl;
    }
    outputFile.close();
}

void writeReadme(const int TEST_CASE, const float scaling_factor, const int max_iteration) {
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