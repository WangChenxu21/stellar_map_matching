/*************************************************************************
    > File Name: point_light.h
    > Author: WangChenxu
    > Mail: wangchenxu21@outlook.com
    > Created Time: Tue 24 Aug 2021 12:19:09 PM UTC
 ************************************************************************/

#ifndef _POINT_LIGHT_H_
#define _POINT_LIGHT_H_

class PointLight {
public:
    int gray;
    double mag;
    int y;
    int x;

    PointLight();
    virtual ~PointLight();
};

#endif
