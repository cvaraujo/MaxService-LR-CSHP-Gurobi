//
// Created by carlos on 27/05/19.
//

#ifndef MS_GRAPH_H
#define MS_GRAPH_H

#include "Arc.h"
#include "Include.h"

using namespace std;
using namespace boost;

struct SPPRC_Graph_Vert {
    SPPRC_Graph_Vert(int n = 0, int c = 0) : num(n), con(c) {}

    int num;
    // Resource consumed
    int con;
};

struct SPPRC_Graph_Arc {
    SPPRC_Graph_Arc(int n = 0, int c = 0, int r = 0) : num(n), cost(c), res(r) {}

    int num;
    // traversal cost
    int cost;
    // traversal resource
    int res;
};

typedef adjacency_list<vecS, vecS, directedS, SPPRC_Graph_Vert, SPPRC_Graph_Arc> SPPRCGraph;

class Graph {
    typedef graph_traits<SPPRCGraph>::vertex_descriptor vertex_descriptor;
    typedef graph_traits<SPPRCGraph>::edge_descriptor edge_descriptor;

    int n, nRemoved, m, paramDelay, paramJitter, paramVariation, paramBandwidth, root, bigMDelay = 0, bigMJitter = 0;

    BoostGraph graphDelaySP, graphJitterSP;
    SPPRCGraph graphDelay, graphJitter;

    vector<VertexDescriptor> predecessors;
    vector<int> distanceDelay, distanceJitter;

public:
    vector<vector<Arc *>> arcs;
    vector<int> terminals, nonTerminals, DuS, delayVector, jitterVector;
    vector<bool> removed;
    vector<vector<vector<bool>>> removedF;
    vector<bool> notAttend;

    Graph(string instance, string param, string outputName);

    void SAE();

    void MVE();

    void showGraph();

    int getN() const;

    int getNAfterRemove() const;

    void setN(int n);

    int getM() const;

    void setM(int m);

    int getParamDelay() const;

    void setParamDelay(int paramDelay);

    int getParamJitter() const;

    void setParamJitter(int paramJitter);

    int getParamVariation() const;

    void setParamVariation(int paramVariation);

    int getParamBandwidth() const;

    void setParamBandwidth(int paramBandwidth);

    int getRoot() const;

    void setRoot(int root);

    int getBigMDelay();

    int getBigMJitter();

    int getShpTerminal(int k);
};


#endif