#include <ctime>
#include <iostream>
#include <math.h>
#include<stdlib.h>
#define ncs 3

using namespace std;

int comp(const void* p1, const void* p2);
double point_dist(double * posStart, double * posEnd);
void matrix_output(double * * matrix, size_t nrs);

/*TODO
 * 1.生成一个n行乘m列的一个矩阵，这里的m列为3列，代表三维坐标，同时可以多列，即代表m维坐标
 * 2.写一个函数，用来计算两个数组之间的距离
 * 3.找出所有需要计算的点的集合，并且生成距离
 * 4.对距离进行圣墟排列
 * 5.计算每一个元素个数
 * 6.加速，我的电脑暂时用不了openmp，嘎嘎
 */

int main() {
    // init matrix
    size_t nrs = 10;

    double ** matrix = (double**) malloc(sizeof(double*) * nrs);

    for ( size_t ix = 0; ix < nrs; ++ix )
        matrix[ix] = (double*) malloc(sizeof(double) * ncs);

    for ( size_t ix = 0; ix < nrs; ++ix ){
        for(size_t ij = 0; ij < ncs; ++ij)
            matrix[ix][ij] = rand() % 32;
    }
    // init distance list
    int comb_size =  nrs * (nrs - 1) / 2;
    double * dists = (double*) malloc(sizeof(double) * comb_size);

    // compute distance
    int dist_pointer = 0;
    for (int i = 0; i < nrs; ++i) {
        for (int j = i + 1; j < nrs; ++j) {
            dists[dist_pointer] = point_dist(matrix[i], matrix[j]);
            dist_pointer += 1;
        }
    }
    // sort distance
    qsort(dists, comb_size, sizeof(double), comp);
    // print matrix
    matrix_output(matrix, nrs);
    // print dist
    cout << "Distance of each point. Ascending sorted." << endl;
    for (int i = 0; i < comb_size; ++i) {
        cout << dists[i] << endl;
    }
    return 0;
}

void matrix_output(double * * matrix, size_t nrs){
    cout << "Coordinates of each point";
    for (int i = 0; i < nrs; ++i) {
        for (int j = 0; j < ncs; ++j) {
            cout << matrix[i][j] << "\t";
        }
        cout << " " << endl;
    }
}

int comp(const void* p1, const void* p2){
    const double * a = (const double *) p1;
    const double * b = (const double *) p2;

    int value = 0;

    if(*a < *b)
        value = -1;
    else if(*a == *b)
        value = 0;
    else value = 1;

    return value;
}

double point_dist(double * posStart, double * posEnd)
{
    double dist = 0;
    for (int i = 0; i < ncs; ++i) {
        dist += pow(posStart[i] - posEnd[i], 2);
    }
    dist = sqrt(dist);
    return dist;
}




