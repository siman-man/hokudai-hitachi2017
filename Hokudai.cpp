#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cassert>
#include <string.h>
#include <vector>

using namespace std;
typedef long long ll;

const int MAX_N = 60;
const int MAX_V = 500;
const int MAX_V_EMB = 3600;

double TIME_LIMIT = 10.0;
const ll CYCLE_PER_SEC = 2700000000;

const int DY[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
const int DX[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

unsigned long long xor128() {
    static unsigned long long rx = 123456789, ry = 362436069, rz = 521288629, rw = 88675123;
    unsigned long long rt = (rx ^ (rx << 11));
    rx = ry;
    ry = rz;
    rz = rw;
    return (rw = (rw ^ (rw >> 19)) ^ (rt ^ (rt >> 8)));
}

unsigned long long int getCycle() {
    unsigned int low, high;
    __asm__ volatile ("rdtsc" : "=a" (low), "=d" (high));
    return ((unsigned long long int) low) | ((unsigned long long int) high << 32);
}

double getTime(unsigned long long int begin_cycle) {
    return (double) (getCycle() - begin_cycle) / CYCLE_PER_SEC;
}

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

struct Node {
    vector<int> neighbors;
};

struct Mapping {
    int from;
    int to;

    Mapping(int from = -1, int to = -1) {
        this->from = from;
        this->to = to;
    }
};

struct Coord {
    int y;
    int x;

    Coord(int y = -1, int x = -1) {
        this->y = y;
        this->x = x;
    }
};

int V, V_emb;
int E, E_emb;
int N;
ll startCycle;

bool edgeMapG[MAX_V][MAX_V];
bool edgeMapGemb[MAX_V_EMB][MAX_V_EMB];
int vertexMapGemb[MAX_N][MAX_N];
int vertexMapping[MAX_V];
int edgeWeight[MAX_V][MAX_V];
Coord coordList[MAX_V_EMB];
vector<Node> nodeList;

class AtCoder {
public:
    void init(vector <Edge> G, vector <Edge> G_emb) {
        memset(edgeMapG, false, sizeof(edgeMapG));
        memset(edgeMapGemb, false, sizeof(edgeMapGemb));
        memset(vertexMapGemb, -1, sizeof(vertexMapGemb));
        memset(vertexMapping, -1, sizeof(vertexMapping));
        memset(edgeWeight, 0, sizeof(edgeWeight));
        N = (int) sqrt(V_emb);

        int n = (int) ceil(sqrt(V));

        fprintf(stderr, "n = %d, N = %d\n", n, N);

        for (int i = 0; i < V; i++) {
            nodeList.push_back(Node());
        }

        for (int i = 0; i < E; i++) {
            Edge edge = G[i];
            edgeWeight[edge.from][edge.to] = edge.weight;
            edgeWeight[edge.to][edge.from] = edge.weight;

            nodeList[edge.from].neighbors.push_back(edge.to);
            nodeList[edge.to].neighbors.push_back(edge.from);
        }

        for (int i = 0; i < E_emb; i++) {
            Edge edge = G_emb[i];
            edgeMapGemb[edge.from][edge.to] = true;
            edgeMapGemb[edge.to][edge.from] = true;
        }

        for (int i = 0; i < MAX_V_EMB; i++) {
            int y = i / N;
            int x = i % N;
            coordList[i] = Coord(y, x);
        }

        int y = 0;
        int x = 1;
        vertexMapping[0] = 0;
        vertexMapGemb[0][0] = 0;

        for (int i = 1; i < V; i++) {
            int z = y * N + x;
            int maxScore = -1;
            int maxId = -1;

            for (int j = 0; j < V; j++) {
                if (vertexMapping[j] != -1) continue;

                vertexMapping[j] = z;
                vertexMapGemb[y][x] = j;

                int score = calcScoreSub(j);
                assert(score >= 0);
                if (maxScore < score) {
                    maxScore = score;
                    maxId = j;
                }

                vertexMapping[j] = -1;
                vertexMapGemb[y][x] = -1;
            }

            vertexMapping[maxId] = z;
            vertexMapGemb[y][x] = maxId;

            x++;
            if (x == n) {
                x = 0;
                y++;
            }
        }
    }

    vector <Mapping> mapping(vector <Edge> G, vector <Edge> G_emb) {
        startCycle = getCycle();
        init(G, G_emb);

        mappingVertex();

        return createAnswer();
    }

    void mappingVertex() {
        int bestScore = calcScore();
        int bestVertexMapping[MAX_V];
        memcpy(bestVertexMapping, vertexMapping, sizeof(vertexMapping));
        double currentTime = getTime(startCycle);
        double remainTime = 0.0;
        ll tryCount = 0;
        int R = 500000;
        double k = 4.0;
        int DD[8] = {-N-1, -N, -N+1, -1, 1, N-1, N, N+1};

        int currentScore = bestScore;

        while (currentTime < TIME_LIMIT) {
            currentTime = getTime(startCycle);
            remainTime = (TIME_LIMIT - currentTime) / TIME_LIMIT;

            int v = xor128() % V;
            int t = vertexMapping[v];
            int i = xor128() % nodeList[v].neighbors.size();
            int j = xor128() % 8;
            int z = vertexMapping[nodeList[v].neighbors[i]] + DD[j];
            int y = z / N;
            int x = z % N;
            if (y < 0 || y >= N || x < 0 || x >= N) continue;
            int m = vertexMapGemb[y][x];

            int diffScore = calcScoreSub(v) + calcScoreSub(m);
            if (m == -1) {
                moveVertex(v, z);
            } else {
                swapVertexMapping(v, m);
            }
            diffScore -= calcScoreSub(v) + calcScoreSub(m);

            int score = currentScore - diffScore;

            if (bestScore < score) {
                bestScore = score;
                memcpy(bestVertexMapping, vertexMapping, sizeof(vertexMapping));
            }

            if (currentScore < score || (diffScore < 30 && xor128() % R < R * exp(-diffScore / (k * remainTime)))) {
                currentScore = score;
            } else {
                if (m == -1) {
                    moveVertex(v, t);
                } else {
                    swapVertexMapping(v, m);
                }
            }

            tryCount++;
        }

        memcpy(vertexMapping, bestVertexMapping, sizeof(bestVertexMapping));
        fprintf(stderr, "BestScore = %d, tryCount = %lld\n", bestScore, tryCount);
    }

    void moveVertex(int v, int z) {
        int t = vertexMapping[v];
        vertexMapping[v] = z;
        Coord c1 = coordList[z];
        Coord c2 = coordList[t];
        vertexMapGemb[c1.y][c1.x] = v;
        vertexMapGemb[c2.y][c2.x] = -1;
    }

    void swapVertexMapping(int i, int j) {
        int t = vertexMapping[i];
        int s = vertexMapping[j];
        vertexMapping[i] = s;
        vertexMapping[j] = t;

        Coord c1 = coordList[t];
        Coord c2 = coordList[s];

        vertexMapGemb[c1.y][c1.x] = j;
        vertexMapGemb[c2.y][c2.x] = i;
    }

    int calcScore() {
        int score = 0;

        for (int i = 0; i < V - 1; i++) {
            for (int j = i + 1; j < V; j++) {
                if (edgeMapGemb[vertexMapping[i]][vertexMapping[j]]) {
                    score += edgeWeight[i][j];
                }
            }
        }

        return score;
    }

    int calcScoreSub(int v) {
        if (v == -1) return 0;
        int score = 0;
        Coord c = coordList[vertexMapping[v]];

        for (int i = 0; i < 8; i++) {
            int ny = c.y + DY[i];
            int nx = c.x + DX[i];
            if (ny < 0 || ny >= N || nx < 0 || nx >= N) continue;
            if (vertexMapGemb[ny][nx] == -1) continue;
            score += edgeWeight[v][vertexMapGemb[ny][nx]];
        }

        return score;
    }

    vector <Mapping> createAnswer() {
        vector <Mapping> ret;

        for (int i = 0; i < V; i++) {
            ret.push_back(Mapping(i + 1, vertexMapping[i] + 1));
        }

        return ret;
    }
};

int main() {
    TIME_LIMIT = 2.0;
    int u, v, w;
    vector <Edge> G;
    cin >> V >> E;
    for (int i = 0; i < E; i++) {
        cin >> u >> v >> w;
        G.push_back(Edge(u - 1, v - 1, w));
    }

    int a, b;
    vector <Edge> G_emb;
    cin >> V_emb >> E_emb;
    for (int i = 0; i < E_emb; i++) {
        cin >> a >> b;
        G_emb.push_back(Edge(a - 1, b - 1));
    }

    AtCoder ac;
    vector <Mapping> mapping = ac.mapping(G, G_emb);

    int n = mapping.size();
    for (int i = 0; i < n; i++) {
        cout << mapping[i].from << " " << mapping[i].to << endl;
    }

    return 0;
}
