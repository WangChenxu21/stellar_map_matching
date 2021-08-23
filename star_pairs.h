/*************************************************************************
    > File Name: star_pairs.h
    > Author: WangChenxu
    > Mail: wangchenxu21@outlook.com
    > Created Time: Mon 23 Aug 2021 01:34:04 AM UTC
 ************************************************************************/

#ifndef _STAR_PAIRS_H_
#define _STAR_PAIRS_H_

#include "star_info.h"

class StarPairs {
public:
    StarInfo s1;
    StarInfo s2;

    StarPairs();
    virtual ~StarPairs();

    void SetDistance();

    double GetCompareMag();
    int GetCompareDisY();
    int GetCompareDisX();
    double GetAbsDis();

private:
    double compare_mag; //两星星等差
    int compare_gray; //两星对应点灰度差
    double abs_dis; //两星绝对距离（欧氏距离）
    int compare_dis_y; //两星Y方向有向距离
    int compare_dis_x; //两星X方向有向距离
};

#endif
