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
    char* star_map_file_name; //星表路径
    char* save_file_name; //图像保存路径
    char* count_file_name; //星等分布统计文件保存路径
    double axes_latitude; //赤纬
    double axes_longitude; //赤经
    int width; //图像宽
    int height; //图像高
    double x_view; //X方向视场角
    double y_view; //Y方向视场角
    double rotate_angle; //旋转角
    double magmax; //最大星等
    double magmin; //最小星等
    double axes_error; //指向误差最大值
    int year; //年
    int month; //月
    int day; //日
    int background_value; //背景起伏最大值
    bool is_color; //产生彩色还是灰度星图
    bool is_position_emendation; //是否进行星位置修正

    StarMap();
    virtual ~StarMap();
    // todo
    void CreateStarImage(); //产生星图
    // todo
    bool CreateStarMap(StarInfo* star_of_starmap);
    // todo
    bool ShowStarMap(unsigned char *LpNewDIBBits, StarInfo *stars_of_star_map, int size, double m_x_view);

private:
    int julian_day_; //儒略日
    LPBYTE img_data_; //图像数据
    PBITMAPINFO img_info_; //图像信息头
    int star_counter_[13]; //星等分布统计数组

    void StarCount(double mag); //星等分布统计
    void SaveStarCountFile(); //保存星等分布统计文件
    // todo
    void SavePic(unsigned char *LpNewDIBBits); //保存图像
    int CoordinateTrans(StarInfo *stars_of_star_map, int size, double m_axes_longitude, double m_axes_latitude, double m_x_view, double m_y_view); //将视场中的星点赤经纬坐标转化为平面坐标

    // 位置修正函数
    // todo
    double Latitude(double x, double y, double z);
    // todo
    double Longitude(double x, double y, double z);
    // todo
    void Aberration(double &ra, double &dec, double de, double C, double D);
    // todo
    void Nutation(double &ra, double &dec, double dely_dy, double B_dB, double de);
    // todo
    void ParameterCalculation(double &dely_dy, double &B_dB, double &de, double &C, double &D);
    // todo
    void Precession(int jd, double &ra, double &dec);
    // todo
    void ProperMotion(int jd, double u_ra, double u_dec, double &ra, double &dec);
    // todo
    int JulianDay();
};

#endif
