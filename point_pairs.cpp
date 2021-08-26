/*************************************************************************
    > File Name: point_pairs.cpp
    > Author: WangChenxu
    > Mail: wangchenxu21@outlook.com
    > Created Time: Thu 26 Aug 2021 01:31:25 PM UTC
 ************************************************************************/

#include <math.h>

#include "point_pairs.h"

PointPairs::PointPairs() {
    compare_dis_x = 0;
    compare_dis_y = 0;
    abs_dis = 0.0;
    compare_mag = 0.0;
}

PointPairs::~PointPairs() {}

double PointPairs::GetCompareMag() {
    return compare_mag;
}

double PointPairs::GetAbsDis() {
    return abs_dis;
}

int PointPairs::GetCompareDisY() {
    return compare_dis_y;
}

int PointPairs::GetCompareDisX() {
    return compare_dis_x;
}

void PointPairs::SetDistance() {
    int x1, y1, x2, y2;
    int gray1, gray2;
    double mag1, mag2;
    x1 = p1.x;
    y1 = p1.y;
    gray1 = p1.gray;
    mag1 = p1.mag;
    x2 = p2.x;
    y2 = p2.y;
    gray2 = p2.gray;
    mag2 = p2.mag;

    compare_dis_x = x1 - x2;
    compare_dis_y = y1 - y2;
    abs_dis = sqrt(double(compare_dis_x * compare_dis_x + compare_dis_y * compare_dis_y));
    compare_mag = mag1 - mag2;
}
