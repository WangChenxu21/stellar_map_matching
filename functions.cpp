/*************************************************************************
    > File Name: functions.cpp
    > Author: WangChenxu
    > Mail: wangchenxu21@outlook.com
    > Created Time: Thu 26 Aug 2021 03:20:27 AM UTC
 ************************************************************************/

#include "functions.h"

bool StarMapMatchLTSHD(cv::Mat &lpDIBBits, long l_width, long l_height, double axes_longitude, double axes_latitude, double rotate_angle, double x_view, double y_view) {
    // 星点分割
    GlobalThresh(lpDIBBits, l_width, l_height);
}

void GlobalThresh(cv::Mat &lpDIBBits, long l_width, long l_height) {
    long histogram[256]; //直方图数组
    int max = 0, min = 255; //像素最大、最小值
    int threshold, threshold_new; //阈值
    unsigned int iter_times = 0; //迭代次数
    long pex_num_1 = 0, sum_1 = 0; //区域1的像素个数，灰度值总和
    long pex_num_2 = 0, sum_2 = 0; //区域2的像素个数，灰度值总和
    int mean_1, mean_2; //两区域的平均灰度

    for (int i = 0; i < l_height; ++i) {
        histogram[i] = 0;
    }

    unsigned int gray;
    // 获得直方图
    for (int i = 0; i < l_height; ++i) {
        for (int j = 0; j < l_width; ++j) {
            gray = lpDIBBits.at<uchar>(i, j);
            ++histogram[gray];

            max = gray > max ? gray : max;
            min = gray < min ? gray : min;
        }
    }

    // 初始化阈值
    threshold_new = int((max + min) / 2 + 0.5);
    threshold = 0;

    // 迭代更新阈值，最大迭代次数100次
    while (threshold_new != threshold && iter_times < 100) {
        threshold = threshold_new;

        // 被阈值分开的两区域的灰度平均值
        for (int i = 0; i < threshold; ++i) {
            pex_num_1 += histogram[i];
            sum_1 += i * histogram[i];
        }
        if (pex_num_1 == 0) {
            mean_1 = 0;
        } else {
            mean1 = int(sum_1 / pex_num_1 + 0.5);
        }

        for (int i = threshold; i < 256; ++i) {
            pex_num_2 += histogram[i];
            sum_2 += i * histogram[i];
        }
        if (pex_num_2 == 0) {
            mean_2 = 0;
        } else {
            mean_2 = int(sum_2 / pex_num_2 + 0.5);
        }

        threshold_new = int((mean_1 + mean_2) / 2 + 0.5);
    }

    // 二值化图像
    for (int i = 0; i < l_height; ++i) {
        for (int j = 0; j < l_width; ++j) {
            gray = lpDIBBits.at<uchar>(i, j);
            if (gray > threshold) {
                lpDIBBits.at<uchar>(i, j) = 255;
            } else {
                lpDIBBits.at<uchar>(i, j) = 0;
            }
        }
    }

    int* connectivity = new int [l_width * l_height];
    int numb;
    for (int i = 0; i < l_height; ++i) {
        for (int j = 0; j < l_width; ++j) {
            connectivity[i * l_width + j] = 0;
        }
    }

    numb = Connectivity8(lpDIBBits, l_width, l_height, connectivity); //8连通域
    Barycenter(lpDIBBits, l_width, l_height, connectivity, numb); //求质心

    delete [] connectivity;
}

int Connectivity8(cv::Mat &lpDIBBits, int width, int height, int* connectivity) {
    int numb = 1;
    long A = width * height;
    int* equal_mark = new int[2 * A];
    int mark_value_1, mark_value_2;
    int q = 0;
    bool used;
    unsigned int gray;

    // 左上像素
    gray = lpDIBBits.at<uchar>(0, 0);
    if (gray != 0) {
        connectivity[0] = numb;
        ++numb;
    }

    // 第一行
    for (int i = 1; i < width; ++i) {
        gray = lpDIBBits.at<uchar>(0, i);
        if (gray != 0) {
            if (connectivity[i - 1] != 0) {
                connectivity[i] = connectivity[i - 1];
            } else {
                connectivity[i] = numb;
                ++numb;
            }
        }
    }

    // 后面的所有行
//    for (int i = 1; i < height; ++i) {
//        for 
//    }
    return -1;
}

bool HausDis(double* HD, PointLight* point_image, int size_image, StarInfo* stars_of_star_map, int size_view) {
    double hd1, hd2, hd; //两点间Hausdorff距离
    int x_img, y_img, x_view, y_view;
    double dis;
    
    for (int i = 0; i < size_image; ++i) {
        for (int j = 0; j < size_view; ++j) {
            HD[i * size_view + j] = 4999;
        }
    }

    long width_img, width_view;
    width_img = min(size_image, 20);
    width_view = min(size_view, 20);

    PointPairs 
}
