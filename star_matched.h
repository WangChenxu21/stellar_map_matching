/*************************************************************************
    > File Name: star_matched.h
    > Author: WangChenxu
    > Mail: wangchenxu21@outlook.com
    > Created Time: Tue 24 Aug 2021 12:31:45 PM UTC
 ************************************************************************/

#ifndef _STAR_MATCHED_H_
#define _STAR_MATCHED_H_

#include "star_info.h"
#include "point_light.h"

class StarMatched {
public:
    double HD;
    PointLight star_image;
    StarInfo star_view;

    StarMatched();
    virtual ~StarMatched();
};

#endif
