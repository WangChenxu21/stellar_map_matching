/*************************************************************************
    > File Name: star_matched.cpp
    > Author: WangChenxu
    > Mail: wangchenxu21@outlook.com
    > Created Time: Tue 24 Aug 2021 12:33:34 PM UTC
 ************************************************************************/

#include "star_matched.h"

StarMatched::StarMatched() {
    HD = 10000;
    star_image.x = -1;
    star_image.y = -1;
    star_image.gray = -1;
    star_image.mag = -100;
    star_view.SetStarNumb(0);
    star_view.SetCoordinateX(-1);
    star_view.SetCoordinateY(-1);
    star_view.SetStarLongitude(-1);
    star_view.SetStarLatitude(-1);
    star_view.SetStarMagnitude(-100);
}

StarMatched::~StarMatched() {}
