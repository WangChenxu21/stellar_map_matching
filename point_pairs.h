/*************************************************************************
    > File Name: point_pairs.h
    > Author: WangChenxu
    > Mail: wangchenxu21@outlook.com
    > Created Time: Thu 26 Aug 2021 01:27:25 PM UTC
 ************************************************************************/

#ifndef _POINT_PAIRS_H_
#define _POINT_PAIRS_H_

#include "point_light.h"

class PointPairs {
public:
    PointLight p1;
    PointLight p2;

    PointPairs();
    virtual ~PointPairs();

    double GetCompareMag();
    double GetAbsDis();
    int GetCompareDisY();
    int GetCompareDisX();
    
    void SetDistance();

private:
    double compare_mag; //两点的星等差
    double abs_dis; //两点的绝对距离（欧氏距离）
    int compare_dis_y; //Y方向有向距离
    int compare_dis_x; //X方向有向距离
};

#endif
