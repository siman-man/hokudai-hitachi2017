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

int V;
int E;

bool edgeMap[MAX_N][MAX_N];
int vertexMapping[MAX_V_EMB];

class AtCoder {
public:
    void init(vector <Edge> G, vector <Edge> G_emb) {
        memset(edgeMap, false, sizeof(edgeMap));
        memset(vertexMapping, -1, sizeof(vertexMapping));

        for (int i = 1; i <= V; i++) {
            vertexMapping[i] = i;
        }
    }

    vector <Mapping> mapping(vector <Edge> G, vector <Edge> G_emb) {
        init(G, G_emb);

        return createAnswer();
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

    int V_emb, E_emb;
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
