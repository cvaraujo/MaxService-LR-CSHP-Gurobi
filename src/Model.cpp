//
// Created by carlos on 06/03/19.
//

#include <chrono>
#include "../headers/Model.h"

Model::Model(Graph *graph) {
    if (graph != nullptr) {
        this->graph = graph;
        // cout << "Begin the preprocessing" << endl;
        // auto start = chrono::steady_clock::now();
        // graph->MVE();
        // graph->SAE();
        // auto end = chrono::steady_clock::now();
        // this->preprocessingTime = chrono::duration_cast<chrono::seconds>(end - start).count();
        // cout << "Finish preprocessing" << endl;
        // initialize();
    } else exit(EXIT_FAILURE);
}

void Model::initialize() {
    int o, d, n = graph->getN(), m = graph->getM();
    try {

        env.set("LogFile", "MS_mip.log");
        env.start();

        f = vector<vector<vector<GRBVar>>>(n, vector<vector<GRBVar>>(n, vector<GRBVar>(n)));
        y = vector<vector<GRBVar>>(n, vector<GRBVar>(n));
        z = vector<GRBVar>(n);

        char name[20];
        for (o = 0; o < n; o++) {
            for (auto *arc : graph->arcs[o]) {
                d = arc->getD();
                sprintf(name, "y_%d_%d", o, d);
                y[o][d] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
                for (int k: graph->DuS) {
                    sprintf(name, "f_%d_%d_%d", o, d, k);
                    this->f[o][d][k] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
                }
            }
        }

        // for (auto q : graph->DuS) {
        //     sprintf(name, "y_0_%d", q);
        //     y[0][q] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
        //     for (auto e : graph->DuS) {
        //         sprintf(name, "f_0_%d_%d", e, q);
        //         this->f[0][e][q] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);    
        //     }
        // }

        for (auto i : graph->terminals) {
            sprintf(name, "z_%d", i);
            z[i] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
        }

        model.update();
    } catch (GRBException &ex) {
        cout << ex.getMessage() << endl;
        cout << ex.getErrorCode() << endl;
        exit(EXIT_FAILURE);
    }
}

void Model::initModelCshp() {
    cout << "Model Creation" << endl;
    objectiveFunctionLrArb();
    rootFlow(), flowConservation(), terminalsFlow();
    maxArcs(), limDelayAndJitter();
}

void Model::objectiveFunctionLrArb() {
    GRBLinExpr objective;
    int j, o, bigMK, bigML, n = graph->getN();

    for (auto k : graph->terminals)
        objective += z[k];

    for (auto k : graph->DuS) 
        for (int i = 0; i < n; i++)
            for (auto arc : graph->arcs[i])
                objective += multipliersRel[i][arc->getD()][k] * (f[i][arc->getD()][k] - y[i][arc->getD()]);

    for (auto k : graph->terminals) 
        for (auto l : graph->terminals)
            if (k != l) {
                bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
                objective += z[k] * -multipliersVar[k][l] * bigMK;

                for (int i = 0; i < n; i++)
                    for (auto arc : graph->arcs[i])
                        objective += multipliersVar[k][l] * (arc->getDelay() * (f[i][arc->getD()][k] - f[i][arc->getD()][l]));    
            }

    for (auto k : graph->terminals) 
        for (auto l : graph->terminals)
            if (k != l) {
                bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
                objective += z[l] * -multipliersVar[k][l] * bigML;
                objective += multipliersVar[k][l] * -graph->getParamVariation();  
            }

/*

    for (auto k : graph->terminals) {
        objective += z[k];
        for (auto l : graph->terminals) {
            if (k != l) {
                bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
                bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
                objective -= (z[k] * (multipliersVar[k][l] * bigMK) + z[l] * (multipliersVar[l][k] * bigML));
                objective -= (multipliersVar[k][l] * graph->getParamVariation());
            }
        }
        for (int i = 0; i < graph->getN(); i++) {
            for (auto arc : graph->arcs[i]) {
                j = arc->getD();
                for (auto l : graph->terminals) {
                    if (k != l) 
                        objective += f[i][j][k] * (arc->getDelay() * (multipliersVar[k][l] - multipliersVar[l][k]));
                }
            }
        }
    }   

*/
    model.setObjective(objective, GRB_MINIMIZE);
    model.update();
    cout << "Objective Function" << endl;
}

void Model::changeCoefObjective(){
	GRBLinExpr objective;
    int j, o, bigMK, bigML;
    for (auto k : graph->terminals) {
        objective += z[k];
        for (auto l : graph->terminals) {
            if (k != l) {
                bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
                bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
                objective -= (z[k] * (multipliersVar[k][l] * bigMK) + z[l] * (multipliersVar[l][k] * bigML));
                objective -= (multipliersVar[k][l] * graph->getParamVariation());
            }
        }
        for (int i = 0; i < graph->getN(); i++) {
            for (auto arc : graph->arcs[i]) {
                j = arc->getD();
                for (auto l : graph->terminals) {
                    if (k != l) objective += f[i][j][k] * (arc->getDelay() * (multipliersVar[k][l] - multipliersVar[l][k]));
                }
            }
        }
    }   

    model.setObjective(objective, GRB_MINIMIZE);
    model.update();
    cout << "Objective Function" << endl;
}

void Model::rootFlow() {
    int o, d, root = graph->getRoot();
    for (auto k : graph->DuS) {
        if (graph->removed[k]) continue;
        GRBLinExpr flowExpr, rootExpr;
        for (o = 0; o < graph->getN(); o++) {
            if (!graph->removed[o]) {
                for (auto *arc : graph->arcs[o]) {
                    d = arc->getD();
                    if (o == root) flowExpr += f[root][d][k];
                    else if (d == root) rootExpr += f[o][root][k];
                }
            }
        }
        model.addConstr((flowExpr - rootExpr) == 1, "root_flow_all_" + to_string(k));
    }
    model.update();
    cout << "Flow on root node" << endl;
}

void Model::flowConservation() {
    int o, d, root = graph->getRoot();
    for (auto k : graph->DuS) {
        if (graph->removed[k]) continue;
        for (int j = 0; j < graph->getN(); j++) {
            if (!graph->removed[j]) {
                if (j != root && j != k) {
                    GRBLinExpr flowIn, flowOut;
                    for (o = 0; o < graph->getN(); o++) {
                        for (auto *arc : graph->arcs[o]) {
                            d = arc->getD();
                            if (o == j) flowOut += f[j][d][k];
                            if (d == j) flowIn += f[o][j][k];
                        }
                    }
                    model.addConstr((flowIn - flowOut) == 0, "flow_conservation_" + to_string(j) + "_" + to_string(k));
                }
            }
        }
    }
    model.update();
    cout << "Flow conservation" << endl;
}

void Model::terminalsFlow() {
    int o, d;
    for (auto k : graph->DuS) {
        if (graph->removed[k]) continue;
        GRBLinExpr flowIn, flowOut;
        for (o = 0; o < graph->getN(); o++) {
            if (!graph->removed[o]) {
                for (auto *arc : graph->arcs[o]) {
                    d = arc->getD();
                    if (o == k) flowOut += f[k][d][k];
                    if (d == k) flowIn += f[o][k][k];
                }
            }
        }
        model.addConstr((flowOut - flowIn) == -1, "flow_on_terminals_" + to_string(k));
    }
    model.update();
    cout << "Flow on terminals" << endl;
}

void Model::maxArcs() {
    GRBLinExpr totalArcs;
    for (int o = 0; o < graph->getN(); o++)
        for (auto *arc : graph->arcs[o]) 
            totalArcs += y[arc->getO()][arc->getD()];

    model.addConstr(totalArcs == (graph->getNAfterRemove() - 1), "maximum_of_arcs");

    model.update();
    cout << "maximum of arcs in the tree" << endl;
}

void Model::limDelayAndJitter() {
    int o, d, paramDelay, paramJitter;
    for (auto k : graph->terminals) {
        GRBLinExpr limDelay, limJitter;
        for (o = 0; o < graph->getN(); o++) {
            for (auto *arc : graph->arcs[o]) {
                d = arc->getD();
                limDelay += arc->getDelay() * f[o][d][k];
                limJitter += arc->getJitter() * f[o][d][k];
            }
        }
        paramDelay = graph->getParamDelay(), paramJitter = graph->getParamJitter();
        model.addConstr(limDelay <= (paramDelay + (graph->getBigMDelay() - paramDelay) * z[k]),
                        "delay_limit_" + to_string(k));
        model.addConstr(limJitter <= (paramJitter + (graph->getBigMJitter() - paramJitter) * z[k]),
                        "jitter_limit_" + to_string(k));
    }
    model.update();
    cout << "Delay and Jitter limits" << endl;
}

bool Model::solve() {
    try {
        model.set("TimeLimit", "3600.0");
        // model.set("OutputFlag", "0");
        model.update();
        model.write("model.lp");
        model.optimize();
        return true;
    } catch (GRBException &ex) {
        cout << ex.getMessage() << endl;
        return false;   
    }
}

void Model::getGradientVar(vector<vector<double>> &gradientVar) {
    int i, j, n = graph->getN(), bigMK, bigML;
    for (int k : graph->terminals) {
        for (int l : graph->terminals) {
            if (k != l) {
                bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
                bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
                gradientVar[k][l] = -graph->getParamVariation() - (z[k].get(GRB_DoubleAttr_X) * bigMK + z[l].get(GRB_DoubleAttr_X) * bigML);
                for (i = 0; i < n; i++) {
                    for (auto *arc : graph->arcs[i]) {
                        j = arc->getD();
                        gradientVar[k][l] += arc->getDelay() * (f[i][j][k].get(GRB_DoubleAttr_X) - f[i][j][l].get(GRB_DoubleAttr_X));
                    }
                }
                if (gradientVar[k][l] > 0) feasible = false;
            }
        }
    }
    // if (!feasible) cause.push_back("v");
}

void Model::getGradientRelation(vector<vector<vector<double>>> &gradientRel) {
    int i, j, n = graph->getN();
    for (auto k : graph->DuS) {
        for (i = 0; i < n; i++) {
            for (auto *arc : graph->arcs[i]) {
                j = arc->getD();
                if (y[i][j].get(GRB_DoubleAttr_X) > 0.5) gradientRel[i][j][k] = -1;
                else gradientRel[i][j][k] = 0;

                if (f[i][j][k].get(GRB_DoubleAttr_X) > 0.5) gradientRel[i][j][k] += 1;
                
                if (gradientRel[i][j][k] > 0) feasible = false;
            }
        }
    }
    // if (!feasible) cause.push_back("r");
}

double Model::getNormRelation(vector<vector<vector<double>>> &gradient) {
    double sum = 0;
    int i, j, n = graph->getN();
    for (auto k : graph->DuS) {
        for (i = 0; i < n; i++) {
            for (auto *arc : graph->arcs[i]) {
                j = arc->getD();
                sum += pow(gradient[i][j][k], 2);
            }
        }
    }
    return sqrt(sum);
}

double Model::getNormVar(vector<vector<double>> &gradient) {
    double sum = 0;
    for (int k : graph->terminals)
        for (int l : graph->terminals)
            if (k != l) sum += pow(gradient[k][l], 2);
    return sqrt(sum);
}

int Model::getOriginalObjValue() {
    int foValue = 0;
    for (auto t : graph->terminals) 
        if (z[t].get(GRB_DoubleAttr_X) > 0.9) {
            // cout << t << ", ";
            foValue++; 
        }
    // cout << endl;
    // getchar();
    return foValue; 
}

bool Model::isFeasible() {
    if (feasible) return true; 
    feasible = true;
    return false;
 /*   else {
        // cout << "Feasible - Upper Bound = " << upperBound << ", (Relaxed) Lower Bound = " << lowerBound << endl;
        // if (upperBound != 6) getchar();
        feasible = true;
        return false;
    }
*/
}

int Model::lagrangean() {
    int progress = 0, iter = 0, n = graph->getN();
    double thetaDelay, normDelay, thetaJitter, normJitter, thetaRel, normRel, thetaVar, normVar, objPpl, originalObj;
    upperBound = int(graph->terminals.size());
    lowerBound = 0;
    // limitants = vector<double>();
    // cause = vector<string>();
    
    vector<vector<double>> gradientVar = vector<vector<double >>(n, vector<double>(n));;
    vector<vector<vector<double>>> gradientRel = vector<vector<vector<double >>>(n, vector<vector<double>>(n, vector<double>(n)));
    
    multipliersVar = vector<vector<double >>(n, vector<double>(n));
    multipliersRel = vector<vector<vector<double >>>(n, vector<vector<double>>(n, vector<double>(n)));
    initialize();
    initModelCshp(); 

    while(iter < maxIter) {

        // objectiveFunctionLrArb();
        if (solve()) {

            // for(auto k : graph->DuS) {
            //     for (int i = 0; i < n; i++) {
            //         for (auto arc : graph->arcs[i]) {
            //             if (f[i][arc->getD()][k].get(GRB_DoubleAttr_X) > 0.9)
            //                 cout << i << ", " << arc->getD() << ", " << k << endl;
            //         }
            //     }
            // }
            // cout << "--------------------" << endl;
            // for (int i = 0; i < n; i++) {
            //     for (auto arc : graph->arcs[i]) {
            //         if (y[i][arc->getD()].get(GRB_DoubleAttr_X) > 0.9)
            //             cout << i << ", " << arc->getD() << endl;
            //     }
            // }

            getGradientVar(gradientVar);
            getGradientRelation(gradientRel);

            objPpl = model.get(GRB_DoubleAttr_ObjVal);
            cout << "PPL: " << objPpl << endl;
            // limitants.push_back(objPpl);

            if (objPpl > lowerBound)
                lowerBound = ceil(objPpl), progress = 0;
             	cout << "(Feasible) Upper Bound = " << upperBound << ", (Relaxed) Lower Bound = " << lowerBound << endl;
                // cout << "gap: " << ((upperBound - lowerBound) / upperBound) * 100 << endl;
            	if ((upperBound - lowerBound) / upperBound <= 0.0001) return upperBound;
            else {
                progress++;
                if (progress == 20) {
                    lambda /= 2;
                	progress = 0;
                }
            }

            originalObj = getOriginalObjValue();
            cout << "Original Obj: " << originalObj << endl;

            if (isFeasible() && originalObj < upperBound) {
                upperBound = originalObj;
                cout << "(Feasible) Upper Bound = " << upperBound << ", (Relaxed) Lower Bound = " << lowerBound << endl;
                // cout << "gap: " << ((upperBound - lowerBound) / upperBound) * 100 << endl;
                // getchar();
                if ((upperBound - lowerBound) / upperBound <= 0.0001) return upperBound;
            }

            normVar = getNormVar(gradientVar);
            normRel = getNormRelation(gradientRel);

            if (normVar == 0) thetaVar = 0;
            else thetaVar = lambda * ((upperBound - objPpl) / pow(normVar, 2));

            if (normRel == 0) thetaRel = 0;
            else thetaRel = lambda * ((upperBound - objPpl) / pow(normRel, 2));

            for (int k : graph->terminals)
                for (int l : graph->terminals)
                    if (k != l) 
                        multipliersVar[k][l] = max(0.0, multipliersVar[k][l] + (gradientVar[k][l] * thetaVar));

            for (auto k : graph->DuS) 
                for (int i = 0; i < n; i++) 
                    for (auto *arc : graph->arcs[i]) 
                            multipliersRel[i][arc->getD()][k] = max(0.0, multipliersRel[i][arc->getD()][k] + gradientRel[i][arc->getD()][k] * thetaRel);

            cout << "(Feasible) Upper Bound = " << upperBound << ", (Relaxed) Lower Bound = " << lowerBound << endl;
            // model.reset(0);
            objectiveFunctionLrArb();
            iter++;
            getchar();
        }
    }

    return upperBound;
}

void Model::showSolution(string instance) {
    try {
        ofstream output;
        output.open(instance, ofstream::app);
        output << "Prep. Time: " << preprocessingTime << endl;
        // double ub = model.get(GRB_DoubleAttr_ObjVal), lb = model.get(GRB_DoubleAttr_ObjBound);
        output << "UB: " << upperBound << endl;
        output << "LB: " << lowerBound << endl;
        if (upperBound != 0) output << "gap: " << (upperBound - lowerBound) / upperBound << endl;
        
        for (auto v : limitants) {
            output << v << endl;
        }        
        output << "-------------------------------" << endl;
        for (auto c : cause) {
            output << c << endl;
        }
        // output << "N. Nodes: " << model.get(GRB_DoubleAttr_NodeCount) << endl;
        // output << "Runtime: " << model.get(GRB_DoubleAttr_Runtime) << endl;

        // output << "----- Solution -----" << endl;
        // for (auto i : graph->terminals)
        //     if (z[i].get(GRB_DoubleAttr_X) > 0.5)
        //         output << i << endl;
        // output << "----- Tree -----" << endl;
        // for (int o = 0; o < graph->getN(); o++) {
        //     for (auto *arc : graph->arcs[o]) {
        //         if (y[o][arc->getD()].get(GRB_DoubleAttr_X) > 0.5) {
        //             output << arc->getO() << ", " << arc->getD() << " = " << arc->getDelay() << endl;
        //         }
        //     }
        // }
        output.close();
    } catch (GRBException &ex) {
        cout << ex.getMessage() << endl;
    }

}
