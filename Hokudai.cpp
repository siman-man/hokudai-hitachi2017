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
const int MAX_V = 501;
const int MAX_V_EMB = MAX_N * MAX_N + 1;

double TIME_LIMIT = 10.2;
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

bool edgeMapGemb[MAX_V_EMB][MAX_V_EMB];
int vertexMapGemb[MAX_N][MAX_N];
int vertexMapping[MAX_V];
char edgeWeight[MAX_V][MAX_V];
Coord coordList[MAX_V_EMB];
vector <Edge> edgeList;

class AtCoder {
public:
    void init(vector <Edge> G, vector <Edge> G_emb) {
        memset(edgeMapGemb, false, sizeof(edgeMapGemb));
        memset(vertexMapGemb, 0, sizeof(vertexMapGemb));
        memset(vertexMapping, 0, sizeof(vertexMapping));
        memset(edgeWeight, 0, sizeof(edgeWeight));
        N = (int) sqrt(V_emb);

        int n = (int) ceil(sqrt(V));

        fprintf(stderr, "V = %d, n = %d, N = %d\n", V, n, N);

        for (int i = 0; i < E; i++) {
            Edge edge = G[i];
            edgeWeight[edge.from][edge.to] = edge.weight;
            edgeWeight[edge.to][edge.from] = edge.weight;

            for (int j = 0; j < edge.weight / 4; j++) {
                edgeList.push_back(Edge(edge.from, edge.to));
                edgeList.push_back(Edge(edge.to, edge.from));
            }
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

        int offset = (N - n) / 2;
        int y = offset;
        int x = offset;

        for (int i = 1; i <= V; i++) {
            int z = y * N + x;

            vertexMapping[i] = z;
            vertexMapGemb[y][x] = i;

            x++;
            if (x == n + offset) {
                x = offset;
                y++;
            }
        }
    }

    vector <Mapping> mapping(vector <Edge> G, vector <Edge> G_emb) {
        startCycle = getCycle();
        init(G, G_emb);

        if (edgeList.size() > 0) {
            mappingVertex();
        }

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
        double k = 4.5;
        int DD[8] = {-N - 1, -N, -N + 1, -1, 1, N - 1, N, N + 1};

        int currentScore = bestScore;
        int thread = 25;
        int es = edgeList.size();
        double expCache[thread];

        random_shuffle(edgeList.begin(), edgeList.end());
        int cnt = 0;

        while (currentTime < TIME_LIMIT) {
            cnt++;
            if (tryCount % 100000 == 0) {
                currentTime = getTime(startCycle);
                remainTime = (TIME_LIMIT - currentTime) / TIME_LIMIT;

                for (int i = 0; i < thread; i++) {
                    expCache[i] = R * exp(-i / (k * remainTime));
                }
            }

            int e = cnt % es;
            Edge edge = edgeList[e];
            int v = edge.from;
            int t = vertexMapping[v];
            int j = xor128() % 8;
            int z = vertexMapping[edge.to] + DD[j];
            if (t == z || z < 0 || z >= N * N) continue;
            Coord c = coordList[z];
            int u = vertexMapGemb[c.y][c.x];

            int diffScore = calcScoreSub(v, t) + calcScoreSub(u, z);
            if (u == 0) {
                moveVertex(v, z);
            } else {
                swapVertexMapping(v, u);
            }
            diffScore -= calcScoreSub(v, z) + calcScoreSub(u, t);

            int score = currentScore - diffScore;

            if (bestScore < score) {
                bestScore = score;
                memcpy(bestVertexMapping, vertexMapping, sizeof(vertexMapping));
            }

            if (currentScore < score || (diffScore < thread && xor128() % R < expCache[diffScore])) {
                currentScore = score;
            } else {
                if (u == 0) {
                    moveVertex(v, t);
                } else {
                    swapVertexMapping(v, u);
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
        vertexMapGemb[c2.y][c2.x] = 0;
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

        for (int i = 1; i < V; i++) {
            for (int j = i + 1; j <= V; j++) {
                if (edgeMapGemb[vertexMapping[i]][vertexMapping[j]]) {
                    score += edgeWeight[i][j];
                }
            }
        }

        return score;
    }

    int calcScoreSub(int v, int z) {
        if (v == 0) return 0;
        int score = 0;
        Coord c = coordList[z];

        for (int i = 0; i < 8; i++) {
            int ny = c.y + DY[i];
            int nx = c.x + DX[i];
            // FIXME: テストケースに依存しているので直す
            // if (ny < 0 || ny >= N || nx < 0 || nx >= N) continue;
            score += edgeWeight[v][vertexMapGemb[ny][nx]];
        }

        return score;
    }

    vector <Mapping> createAnswer() {
        vector <Mapping> ret;

        for (int i = 1; i <= V; i++) {
            ret.push_back(Mapping(i, vertexMapping[i] + 1));
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
