/*************************************************************************
    > File Name: star_pairs.cpp
    > Author: WangChenxu
    > Mail: wangchenxu21@outlook.com
    > Created Time: Mon 23 Aug 2021 01:42:35 AM UTC
 ************************************************************************/

#include <math.h>

#include "star_pairs.h"

StarPairs::StarPairs() {
    compare_mag = 0.0;
    compare_gray = 0;
    abs_dis = 0.0;
    compare_dis_y = 0;
    compare_dis_x = 0;
}

StarPairs::~StarPairs() {}

void StarPairs::SetDistance() {
    int x1, y1, x2, y2;
    double mag1, mag2;
    x1 = s1.GetCoordinateX();
    y1 = s1.GetCoordinateY();
    mag1 = s1.GetStarMagnitude();
    x2 = s2.GetCoordinateX();
    y2 = s2.GetCoordinateY();
    mag2 = s2.GetStarMagnitude();

    compare_dis_x = x1 - x2;
    compare_dis_y = y1 - y2;
    abs_dis = sqrt(double(compare_dis_x * compare_dis_x + compare_dis_y * compare_dis_y));
    compare_mag = mag1 - mag2;
}

double StarPairs::GetCompareMag() {
    return compare_mag;
}

int StarPairs::GetCompareDisY() {
    return compare_dis_y;
}

int StarPairs::GetCompareDisX() {
    return compare_dis_x;
}

double StarPairs::GetAbsDis() {
    return abs_dis;
}
