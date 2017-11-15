#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string.h>
#include <vector>

using namespace std;

const int MAX_N = 60 + 1;
const int MAX_V = 500 + 1;
const int MAX_V_EMB = 3600 + 1;

struct Edge {
    int from;
    int to;
    int weight;

    Edge(int from = -1, int to = -1, int weight = -1) {
        this->from = from;
        this->to = to;
        this->weight = weight;
    }
};

struct Mapping {
    int from;
    int to;

    Mapping(int from = -1, int to = -1) {
        this->from = from;
        this->to = to;
    }
};

int V, V_emb;
int E, E_emb;

bool edgeMapG[MAX_V][MAX_V];
bool edgeMapGemb[MAX_V_EMB][MAX_V_EMB];
int vertexMapping[MAX_V_EMB];
int edgeWeight[MAX_V][MAX_V];

class AtCoder {
public:
    void init(vector <Edge> G, vector <Edge> G_emb) {
        memset(edgeMapG, false, sizeof(edgeMapG));
        memset(edgeMapGemb, false, sizeof(edgeMapGemb));
        memset(vertexMapping, -1, sizeof(vertexMapping));
        memset(edgeWeight, 0, sizeof(edgeWeight));

        for (int i = 1; i <= V; i++) {
            vertexMapping[i] = i;
        }

        for (int i = 0; i < E; i++) {
            Edge edge = G[i];
            edgeWeight[edge.from][edge.to] = edge.weight;
            edgeWeight[edge.to][edge.from] = edge.weight;
        }

        for (int i = 0; i < E_emb; i++) {
            Edge edge = G_emb[i];
            edgeMapGemb[edge.from][edge.to] = true;
            edgeMapGemb[edge.to][edge.from] = true;
        }
    }

    vector <Mapping> mapping(vector <Edge> G, vector <Edge> G_emb) {
        init(G, G_emb);

        int score = calcScore();
        fprintf(stderr, "Score = %d\n", score);

        return createAnswer();
    }

    int calcScore() {
        int score = 0;

        for (int i = 1; i < V; i++) {
            for (int j = i + 1; j <= V; j++) {
                if (edgeMapGemb[vertexMapping[i]][vertexMapping[j]]) {
                    score += edgeWeight[i][j];
                }
            }
        }

        return score;
    }

    vector <Mapping> createAnswer() {
        vector <Mapping> ret;

        for (int i = 1; i <= V; i++) {
            ret.push_back(Mapping(i, vertexMapping[i]));
        }

        return ret;
    }
};

int main() {
    int u, v, w;
    vector <Edge> G;
    cin >> V >> E;
    for (int i = 0; i < E; i++) {
        cin >> u >> v >> w;
        G.push_back(Edge(u, v, w));
    }

    int a, b;
    vector <Edge> G_emb;
    cin >> V_emb >> E_emb;
    for (int i = 0; i < E_emb; i++) {
        cin >> a >> b;
        G_emb.push_back(Edge(a, b));
    }

    AtCoder ac;
    vector <Mapping> mapping = ac.mapping(G, G_emb);

    int n = mapping.size();
    for (int i = 0; i < n; i++) {
        cout << mapping[i].from << " " << mapping[i].to << endl;
    }

    return 0;
}
