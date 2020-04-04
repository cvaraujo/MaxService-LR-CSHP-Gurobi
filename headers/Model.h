//
// Created by carlos on 06/03/19.
//

#ifndef MRP_MODEL_H
#define MRP_MODEL_H

#include "Include.h"
#include "Graph.h"

class Model {
    Graph *graph;
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    vector<vector<vector<GRBVar>>> f;
    vector<vector<GRBVar>> y;
    vector<GRBVar> z;
    vector<double> multipliersDelay, multipliersJitter;
    vector<vector<double>> multipliersVar;
    vector<vector<vector<double>>> multipliersRel;
    vector<string> cause;
    vector<double> limitants;

    bool feasible = true;
    int preprocessingTime, maxIter = 450, upperBound;
    double lowerBound, lambda = 1.5;

    void preprocessing();

    void changeCoefObjective();

    void objectiveFunction();

    void rootFlow();

    void flowConservation();

    void terminalsFlow();

    void relXandY();

    void onlyWay();

    void maxArcs();

    void limDelayAndJitter();

    void limVariation();

    void primeToTerminals();

    void nonTerminalsLeafs();

    void getGradient(vector<double> &gradientDelay, vector<double> &gradientJitter);

    void getGradientVar(vector<vector<double>> &gradientVar);

    void getGradientRelation(vector<vector<vector<double>>> &gradientRel);

    double getNorm(vector<double> &gradient);

    double getNormRelation(vector<vector<vector<double>>> &gradient);

    double getNormVar(vector<vector<double>> &gradient);

    int getOriginalObjValue();

    bool isFeasible();

    void objectiveFunctionLrArb();

public:
    Model(Graph *graph);

    void initialize();

    void initModel();

    void initModelCshp();

    bool solve();

    void showSolution(string instance);

    int lagrangean();

};


#endif //MRP_MODEL_H
