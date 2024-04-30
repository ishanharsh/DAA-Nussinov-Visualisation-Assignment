/**
 * @file nussinov.cpp
 * @brief Implementation of the Nussinov algorithm for RNA secondary structure prediction.
 * 
 * The Nussinov algorithm is a fundamental tool for RNA secondary structure prediction, providing insights into the structural properties and functional roles of RNA molecules. By leveraging dynamic programming principles, the algorithm efficiently computes the optimal secondary structure, enabling a wide range of applications in biological research and biomedical engineering.
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>

using namespace std;

/**
 * @brief Check if two nucleotides can pair with each other.
 * 
 * @param a First nucleotide
 * @param b Second nucleotide
 * @return True if the nucleotides can form a base pair, False otherwise
 */
bool pair_check(char a, char b) {
    return (a == 'A' && b == 'U') || (a == 'U' && b == 'A') || (a == 'C' && b == 'G') || (a == 'G' && b == 'C');
}

/**
 * @brief Calculate the optimal score of folding a subsequence from index i to j.
 * 
 * @param i Start index of the subsequence
 * @param j End index of the subsequence
 * @param sequence RNA sequence
 * @param DP Dynamic programming table
 * @return The optimal score of folding the subsequence
 */
int OPT(int i, int j, const string& sequence, vector<vector<int>>& DP) {
    if (i >= j - 4) {
        return 0;
    } else {
        int unpaired = OPT(i, j - 1, sequence, DP);
        vector<int> pairing;
        for (int t = i; t < j - 4; ++t) {
            if (pair_check(sequence[t], sequence[j])) {
                pairing.push_back(1 + OPT(i, t - 1, sequence, DP) + OPT(t + 1, j - 1, sequence, DP));
            }
        }
        if (pairing.empty()) {
            pairing.push_back(0);
        }
        int paired = *max_element(pairing.begin(), pairing.end());
        return max(unpaired, paired);
    }
}

/**
 * @brief Traceback to find the optimal secondary structure.
 * 
 * @param i Start index of the subsequence
 * @param j End index of the subsequence
 * @param structure Vector to store base pairs forming the secondary structure
 * @param DP Dynamic programming table
 * @param sequence RNA sequence
 */
void traceback(int i, int j, vector<pair<int, int>>& structure, const vector<vector<int>>& DP, const string& sequence) {
    if (j <= i) {
        return;
    } else if (DP[i][j] == DP[i][j - 1]) {
        traceback(i, j - 1, structure, DP, sequence);
    } else {
        for (int k = i; k < j - 4; ++k) {
            if (pair_check(sequence[k], sequence[j])) {
                if (k - 1 < 0) {
                    if (DP[i][j] == DP[k + 1][j - 1] + 1) {
                        structure.emplace_back(k, j);
                        traceback(k + 1, j - 1, structure, DP, sequence);
                        break;
                    }
                } else if (DP[i][j] == DP[i][k - 1] + DP[k + 1][j - 1] + 1) {
                    structure.emplace_back(k, j);
                    traceback(i, k - 1, structure, DP, sequence);
                    traceback(k + 1, j - 1, structure, DP, sequence);
                    break;
                }
            }
        }
    }
}

/**
 * @brief Generate dot-bracket notation for the predicted secondary structure.
 * 
 * @param sequence RNA sequence
 * @param structure Vector containing base pairs forming the secondary structure
 * @return Dot-bracket notation of the secondary structure
 */
string write_structure(const string& sequence, const vector<pair<int, int>>& structure) {
    string dot_bracket(sequence.length(), '.');
    for (auto s : structure) {
        dot_bracket[min(s.first, s.second)] = '(';
        dot_bracket[max(s.first, s.second)] = ')';
    }
    return dot_bracket;
}

/**
 * @brief Initialize the dynamic programming table.
 * 
 * Helps us generate the 2D dynamic programming table useful for memoization of similar recursion function calls;
 * @param N Length of the RNA sequence
 * @return Initialized dynamic programming table
 */
vector<vector<int>> initialize(int N) {
    vector<vector<int>> DP(N, vector<int>(N, 0));
    for (int k = 0; k < 4; ++k) {
        for (int i = 0; i < N - k; ++i) {
            int j = i + k;
            DP[i][j] = 0;
        }
    }
    return DP;
}

/**
 * @brief Predict the secondary structure of an RNA sequence using the Nussinov algorithm.
 * 
 * Its is the function where the whole algorithm is set up and recursion calls are made.
 * @param sequence RNA sequence
 * @return Vector containing base pairs forming the secondary structure
 */
vector<pair<int, int>> nussinov(const string& sequence) {
    int N = sequence.length();
    vector<vector<int>> DP = initialize(N);
    vector<pair<int, int>> structure;

    for (int k = 5; k < N; ++k) {
        for (int i = 0; i < N - k; ++i) {
            int j = i + k;
            DP[i][j] = OPT(i, j, sequence, DP);
        }
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < i; ++j) {
            DP[i][j] = DP[j][i];
        }
    }

    traceback(0, N - 1, structure, DP, sequence);
    for(int i=0;i<N;i++){
    for(int j=0;j<N;j++)
    cout<<DP[i][j]<<" ";
    cout<<endl;
    }
    return structure;
}

/**
 * @brief Main Driver function
 * 
 * This is the main function of the program and consists code to take input of the primary sequence and call the necessary algorithm functions and prints out the indices of the base pairs and the dot bracket notation for the secondary structure.
 * @return Exit status
 */
int main() 
{
    string sequence;
    cin >> sequence;
    auto structure = nussinov(sequence);
    string dot_brac = write_structure(sequence,structure);
    cout<<structure.size()<<endl;
    for (auto s : structure) {
        cout << "(" << s.first << ", " << s.second << ") ";
    }
    cout << endl << dot_brac << endl;
    return 0;
}
