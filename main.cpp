#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <float.h>
#include <math.h>
#include <map>
#include <set>
#include <vector>
//#include <mpi.h>

#define DISTANCE_THRESHOLD 10
#define MODULARITY_CHANGE_THRESHOLD 0.001

int findInC(std::vector<std::set<int> > C, int v, int needle) {
    for (int i = 0; i < v; i++) {
        if (C[i].find(needle) != C[i].end()) {
            return i;
        }
    }

    assert(false);

    return -1;
}

void louvain(int v, double **L, std::vector<std::set<int> > C) {
    // For all V, determine the weighted degree of vertex i, `k[i]`
    double *k = (double *) calloc(static_cast<size_t>(v), sizeof(double));
    for (int i = 0; i < v; i++) { // for all vertices
        for (int j = 0; j < v; j++) {// for all connected vertices
            if (i != j && L[i][j] != -1) { // TODO: remove this for the general case
                k[i] += L[i][j];
            }
        }
    }

    // Calculate the sum of all edge weights, `m`
    double m = 0;
    for (int i = 0; i < v; i++) {
        m += k[i];
    }
    m = 0.5 * m;

    int *l = (int *) malloc(v * sizeof(int)); // Community Labels
    double *C_int = (double *) malloc(v * sizeof(double)); // intra-community edges
    double *C_tot = (double *) malloc(v * sizeof(double)); // total-community edges

    for (int i = 0; i < v; i++) {
        C[i].insert(i); // initialize each community to contain the node with the same namesake
        l[i] = i; // label each community
        C_int[i] = 0; // initialize intra-community edges
        C_tot[i] = 0; // initialize total-community edges

        for (int j = 0; j < v; j++) {
            if (i != j && L[i][j] != -1) { // TODO: remove this for the general case
                C_tot[i] += L[i][j];
            }
        }
    }

    double Q_c = 0;
    double Q_p = DBL_MAX;
    int epoch = 0;

    while (true) {
        epoch++;
        for (int i = 0; i < v; i++) {
            int C_old_index = findInC(C, v, i);
            std::set<int> N_i = std::set<int>(); // Communities to try
            N_i.insert(C_old_index);
            for (int j = 0; j < v; j++) { // add the community of each adjacent point `j` to `i`
                if (i != j && L[i][j] != -1) { // TODO: remove this for the general case
                    N_i.insert(findInC(C, v, j));
                }
            }
            double maxGain = 0;
            int C_new_index = C_old_index;
            for (auto it = N_i.begin(); it != N_i.end(); it++) {
                int c = *it;
                /*
                 * Caculate change in Q_{i->c}
                 */
                // TODO: probably remove these loops, or make them more compact (as in remove the need for them duplicated later on in the codebase
                for (int cc = 0; cc < C.size(); cc++) {
                    if (!C[cc].empty()) {
                        C_int[cc] = 0;
                        C_tot[cc] = 0;
                    }
                }
                for (int ii = 0; ii < v; ii++) {
                    for (int j = 0; j < v; j++) {
                        if (i != j && L[i][j] != -1) { // TODO: remove this for the general case
                            if (C[ii] == C[j]) {
                                C_int[ii] += L[ii][j];
                                C_tot[ii] += L[ii][j];
                            } else {
                                C_tot[ii] += L[ii][j];
                                C_tot[j] += L[ii][j];
                            }
                        }
                    }
                }
                double e_xx = 0;
                double a2_x = 0;
                for (int cc = 0; cc < C.size(); cc++) {
                    e_xx += C_int[cc];
                    a2_x += pow(C_tot[cc], 2);
                }
                // END TODO: removal
                double e_ij = C_int[c];
                double e_ii = C_int[i];

                double curGain = (e_ij - (e_ii - C_int[c])/ m)
                        + ((2 * k[i] * (C_tot[i] - L[i][c]) - 2 * k[i] * C_tot[c]) / pow(2 * m, 2));

                if (curGain > maxGain || (curGain == maxGain && l[c] < l[C_new_index])) {
                    maxGain = curGain;
                    C_new_index = c;
                }
            }
            if (maxGain > 0) {
                C[C_old_index].erase(i);
                C[C_new_index].insert(i);
            }
        }

        for (int c = 0; c < C.size(); c++) {
            if (!C[c].empty()) {
                C_int[c] = 0;
                C_tot[c] = 0;
            }
        }

        for (int i = 0; i < v; i++) {
            for (int j = 0; j < v; j++) {
                if (findInC(C, v, i) == findInC(C, v, j)) {
                    C_int[i] += L[i][j];
                    C_tot[i] += L[i][j];
                } else {
                    C_tot[i] += L[i][j];
                    C_tot[j] += L[i][j];
                }
            }
        }
        double e_xx = 0;
        double a2_x = 0;
        for (int c = 0; c < C.size(); c++) {
            if (!C[c].empty()) {
                e_xx += C_int[c];
                a2_x += pow(C_tot[c], 2);
            }
        }
        Q_c = (e_xx / m) - (a2_x / pow(2 * m, 2));

        // Print the current state of learning
        // format: "identifier,epoch,current modularity,change in modularity"
        printf("epoch,%i,%lf,%lf\n", epoch, Q_c, abs((Q_c - Q_p) / Q_p));

        if (abs((Q_c - Q_p) / Q_p) < MODULARITY_CHANGE_THRESHOLD) {
            for (int c = 0; c < C.size(); c++) {
                if (!C[c].empty()) {
                    for (auto it = C[c].begin(); it != C[c].end(); it++) {
                        //* Show the final selected communities from the algorithm
                        //* Format: "identifier,community_id,vertex_index"
                        printf("selected_communities,%i,%i\n", c, *it);
                    }
                }
            }
            break;
        } else {
            Q_p = Q_c;
        }
    }
}

double distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

int main(int argc, char *argv[]) {
    int pointsToCreate = 200;

    int rank = 0, total_tasks;
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Status stat;
    if (rank == 0) { // TODO: remove!

        srandom(rank);

        double *points = (double *) malloc(pointsToCreate * 2 * sizeof(double));

        /// Generate Points in Euclidean Space to use
        if (rank == 0) {
            for (int i = 0; i < pointsToCreate; i++) {
                points[i * 2] = (random() / (double) RAND_MAX) * 100;
                points[(i * 2) + 1] = (random() / (double) RAND_MAX) * 100;

                //* Print the points generated by the current process.
                //* Format: "identifier,i,x,y"
                //printf("generated_points,%i,%lf,%lf\n", i, points[i * 2], points[(i * 2) + 1]);
            }
        }

        struct timespec start, end;
        if (rank == 0) {
            clock_gettime(CLOCK_MONOTONIC_RAW, &start);
        }

        /// Calculate Graph values Gamma (L) and omega (w) all combined into one data structure
        double **L = (double **) malloc(pointsToCreate * sizeof(double *));
        for (int i = 0; i < pointsToCreate; i++) {
            L[i] = (double *) malloc(pointsToCreate * sizeof(double));

            for (int j = 0; j < pointsToCreate; j++) {
                double d = distance(points[i * 2], points[(i * 2) + 1], points[j * 2], points[(j * 2) + 1]);
                L[i][j] = d < DISTANCE_THRESHOLD ? d : -1;
//                L[i][j] = d < DISTANCE_THRESHOLD ? 1 : -1;


                if (L[i][j] != -1) {
                    //* Print the adjacency matrix values (and the associated weights)
                    //* Format: "identifier,i,j,w(i,j)"
                    //printf("adjacency_matrix,%i,%i,%lf\n", i, j, L[i][j]);
                }
            }
        }
        /// Create Communities data structure
        std::vector<std::set<int> > C;
        C.resize(static_cast<unsigned long>(pointsToCreate), std::set<int>());

        /// Run Parallel Louvain here.
        louvain(pointsToCreate, L, C);

        if (rank == 0) {
            clock_gettime(CLOCK_MONOTONIC_RAW, &end);
            //* Print the time take for the given number of processors used and number of cities total
            //* Format: "identifier,num processors,num cities,time(ms)"
            printf("time_taken_ms,%i,%lf\n", total_tasks,
                   (1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec) / 1e6);
        }
    } // TODO: remove

//    MPI_Finalize();
    return 0;
}
