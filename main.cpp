#include <sys/stat.h>
#include "headers/Graph.h"
#include "headers/Model.h"

int main(int argc, const char *argv[]) {
    if (argc < 4) {
        return 0;
    } else {
        mkdir("results", 0777);
        auto *graph = new Graph(argv[1], argv[2], argv[3]);
        auto *model = new Model(graph);
        model->lagrangean();
        // model->showSolution(argv[3]);
    }

    return 0;
}
