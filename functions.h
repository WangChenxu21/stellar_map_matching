/*************************************************************************
    > File Name: functions.h
    > Author: WangChenxu
    > Mail: wangchenxu21@outlook.com
    > Created Time: Thu 26 Aug 2021 03:17:10 AM UTC
 ************************************************************************/

#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include <opencv2/opencv.h>

#include "point_pairs.h"
#include "star_pairs.h"
#include "star_matched.h"
#include "star_map.h"

bool StarMapMatchLTSHD(cv::Mat &lpDIBBits, long l_width, long l_height, double axes_longitude, double axes_latitude, double rotate_angle, double x_view, double y_view);

// 星点分割
void GlobalThresh(cv::Mat &lpDIBBits, long l_width, long l_height);
// 8连通域
int Connectivity8(cv::Mat &lpDIBBits, int width, int height, int* connectivity);
// 待校正图像中点和基准图像中的HD
bool HausDis(double* HD, PointLight* point_image, int size_image, StarInfo* stars_of_star_map, int size_view);

#endif
