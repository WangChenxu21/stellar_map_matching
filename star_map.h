/*************************************************************************
    > File Name: StarMap.h
    > Author: WangChenxu
    > Mail: wangchenxu21@outlook.com
    > Created Time: Fri 20 Aug 2021 01:10:16 AM UTC
 ************************************************************************/

#ifndef _STARMAP_H_
#define _STARMAP_H_

class StarMap {
public:
    bool CreateStarMap(StarInfo* star_of_starmap);

private:
    int julian_day_; //儒略日
    LPBYTE img_data_; //图像数据
    PBITMAPINFO img_info_; //图像信息头
    int star_counter_[13]; //星等分布统计数组
};

#endif
