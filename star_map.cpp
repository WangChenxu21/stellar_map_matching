/*************************************************************************
    > File Name: star_map.cpp
    > Author: WangChenxu
    > Mail: wangchenxu21@outlook.com
    > Created Time: Fri 20 Aug 2021 08:02:27 AM UTC
 ************************************************************************/

#include <fstream>

#include "star_map.h"

StarMap::StarMap() {
    axes_error = 0.0;
    year = 2007;
    month = 12;
    day = 5;
    height = 1024;
    width = 1024;
    magmax = 8;
    magmin = 0;
    x_view = 15;
    y_view = 15;
    axes_latitude = 30;
    axes_longitude = 10.1;
    rotate_angle = 0;
    is_color = false;
    is_position_emendation = false;
    star_map_file_name = "saoN3.dat";
    save_file_name = "Star.bmp";
    count_file_name = "Star.StarCounter.txt";
    for (int i = 0; i < 13; ++i) {
        star_counter_[i] = 0;
    }
}

StarMap::~StarMap() {}

void StarMap::StarCount(double mag) {
    int magcount = int(mag) + 1;
    if (magcount < 1) {
        magcount = 0;
    } else if (magcount > 11) {
        magcount = 12;
    }
    star_counter_[magcount] += 1;
}

void StarMap::SaveStarCountFile() {
    ofstream outfile(star_count_file_name, ios::app);
    for (int i = 0; i < 13; ++i) {
        if (i == 0) {
            outfile << "星等小于0的星个数为" << star_counter_[0] << std::endl;
        } else if (i == 12) {
            outfile << "星等大于11的星个数为" << star_counter_[12] << std::endl;
        } else {
            outfile << "星等从" << (i-1) << "到" << i << "的星个数为" << star_counter_[i] << std::endl;
        }
    }
    outfile.close();
}

int StarMap::CoordinateTrans(StarInfo *stars_of_star_map, int size, double m_axes_longitude, double m_axes_latitude, double m_x_view, double x_y_view) {
    int l_width = width, l_height = height;
    double m_x, m_y;
    double denominator, molecule_x, molecule_y;
    double alpha, delta;
    double alpha0 = (m_axes_longitude * 15) * pi / 180; //视轴赤经
    double delta0 = m_axes_latitude * pi / 180; //视轴赤纬
    int x_rotate, y_rotate;
    int star_count = 0; //坐标转换成功的星数
    double u_ra, u_dec;
    double dely_dy = 0, B_dB = 0, de = 0, C = 0, D = 0;

    if (is_position_emendation) {
        ParameterCalculation(dely_dy, B_dB, de, C, D);
    }

    for (int i = 0; i < size; ++i) {
        alpha = (stars_of_star_map[i].GetStarLongitude()) * 15 * pi / 180; //获得星点赤经
        delta = (stars_of_star_map[i].GetStarLatitude()) * pi / 180; //获得星点赤纬

        if (is_position_emendation) {
            u_ra = stars_of_star_map[i].GetStarProperMotionRA();
            u_dec = stars_of_star_map[i].GetStarProperMotionDEC();
            ProperMotion(julian_day_, u_ra, u_dec, alpha, delta);
            Aberration(alpha, delta, de, C, D);
            Precession(julian_day_, alpha, delta);
            Nutation(alpha, delta, dely_dy, B_dB, de);
        }

        denominator = sin(delta) * sin(delta0) + cois(delta) * cos(delta0) * cos(alpha - alpha0);
        molecule_x = cos(delta) * sin(alpha - alpha0);
        molecule_y = sin(delta) * cos(delta0) - cos(delta) * sin(delta0) * cos(alpha - alpha0);

        m_x = ((l_width / (2 * tan(m_x_view * pi / 360))) * (molecule_x / denominator) + 0.5);
        m_y = ((l_height / (2 * tan(m_y_view * pi / 360))) * (molecule_y / denominator) + 0.5);

        x_rotate = int(m_x * cos(rotate_angle_ * pi / 180) + m_y * sin(rotate_angle_ * pi / 180));
        y_rotate = int(m_y * cos(rotate_angle_ * pi / 180) - m_x * sin(rotate_angle_ * pi / 180));

        //转换到像素坐标系
        if ((x_rotate > -l_width / 2) && (x_rotate <= l_width / 2) && (y_rotate >= -l_height / 2) && (y_rotate < l_height / 2)) {
            stars_of_star_map[star_count] = stars_of_star_map[i];
            stars_of_star_map[star_count].SetCoordinateX(l_height / 2 - y_rotate);
            stars_of_star_map[star_count].SetCoordinateY(l_width / 2 + x_rotate);
            ++star_count;
        }
    }

    return star_count;
}
