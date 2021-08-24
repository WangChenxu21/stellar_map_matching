/*************************************************************************
    > File Name: star_map.cpp
    > Author: WangChenxu
    > Mail: wangchenxu21@outlook.com
    > Created Time: Fri 20 Aug 2021 08:02:27 AM UTC
 ************************************************************************/

#include <fstream>

#include "star_map.h"

#define pi 3.1415926f

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
    std::ofstream outfile(count_file_name, std::ios::app);
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

int StarMap::CoordinateTrans(StarInfo *stars_of_star_map, int size, double m_axes_longitude, double m_axes_latitude, double m_x_view, double m_y_view) {
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

        denominator = sin(delta) * sin(delta0) + cos(delta) * cos(delta0) * cos(alpha - alpha0);
        molecule_x = cos(delta) * sin(alpha - alpha0);
        molecule_y = sin(delta) * cos(delta0) - cos(delta) * sin(delta0) * cos(alpha - alpha0);

        m_x = ((l_width / (2 * tan(m_x_view * pi / 360))) * (molecule_x / denominator) + 0.5);
        m_y = ((l_height / (2 * tan(m_y_view * pi / 360))) * (molecule_y / denominator) + 0.5);

        x_rotate = int(m_x * cos(rotate_angle * pi / 180) + m_y * sin(rotate_angle * pi / 180));
        y_rotate = int(m_y * cos(rotate_angle * pi / 180) - m_x * sin(rotate_angle * pi / 180));

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

double StarMap::Latitude(double x, double y, double z) {
    double temp;
    if (fabs(z) == 1) {
        if (z > 0) {
            temp = 90;
        } else {
            temp = -90;
        }
    } else {
        temp = atan(z / sqrt(1 - pow(z, 2))) * 180 / pi;
    }

    return temp;
}

double StarMap::Longitude(double x, double y, double z) {
    double temp;
    if (x == 0 && y == 0) {
        temp = 0;
    } else if (x == 0 && y != 0) {
        if (y > 0) {
            temp = 6;
        } else {
            temp = 18;
        }
    } else {
        double deg = atan(y / x) * 180 / pi;
        if (x >= 0 && y >= 0) {
            temp = deg / 15;
        } else if (x < 0) {
            temp = (deg + 180) / 15;
        } else {
            temp = (deg + 360) / 15;
        }
    }

    return temp;
}

void StarMap::Aberration(double &ra, double &dec, double de, double C, double D) {
    double sin_ra = sin(ra);
    double cos_ra = cos(ra);
    double sec_dec = 1 / cos(dec);
    double sin_dec = sin(dec);
    double c = cos_ra * sec_dec;
    double d = sin_ra * sec_dec;
    double cc = tan(de) * cos(dec) - sin_ra * sin_dec;
    double dd = cos_ra * sin_dec;
    double del_ra = C * c + D * d;
    double del_dec = C * cc + D * dd;

    ra += del_ra;
    dec += del_dec;
}

void StarMap::Nutation(double &ra, double &dec, double dely_dy, double B_dB, double de) {
    double cos_de = cos(de);
    double sin_de = sin(de);
    double sin_ra = sin(ra);
    double cos_ra = cos(ra);
    double tan_dec = tan(dec);
    double del_ra = dely_dy * (cos_de + sin_de * sin_ra * tan_dec) + B_dB * tan_dec * cos_ra;
    double del_dec = dely_dy * sin_de * cos_ra - B_dB * sin_ra;

    ra += del_ra;
    dec += del_dec;
}

void StarMap::ParameterCalculation(double &dely_dy, double &B_dB, double &de, double &C, double &D) {
    double M = 6.240059966692;
    double ed=0.01670862;
    double E = M;
    double E1, f1, dif_f1;
    
    for (int i = 0; i < 5; ++i) {
        f1 = (E - ed * sin(E)) - M;
        dif_f1 = 1 - ed * cos(E);
        E1 = E - (f1 / dif_f1);
        int note = 0;
        if ((E1 - E) < (1 * 10 ^ (-10))) {
            if ((E1 - E) > (-1 * 10 ^ (-10))) {
                E = E1;
                note = 1;
            }
        }
        if (note == 1) {
            break;
        }
        E = E1;
    }

    double sin_E = sin(E);
    double cos_E = cos(E);
    double fd = 2 * atan(sqrt((1 + ed) / (1 - ed)) * (1 - cos_E) / sin_E);
    double ud = 4.9381882858862 + fd;
    double sin_ud = sin(ud);
    double cos_ud = cos(ud);
    double b = 2.355548393544;
    double c = 6.240035939326;
    double d = 1.627901933972;
    double e = 4.113753719200;
    double f = 2.182438624361;
    double g = 2 * (d + f);
    double h = g - e;
    
    dely_dy = (0.2062 * sin(f + f) - 17.1996 * sin(f) - 1.3187 * sin(h) + 0.1426 * sin(c) - 0.0517 * sin(h + c) + 0.0217 * sin(h - c) + 0.0129 * sin(h - f) - 0.2274 * sin(g) + 0.0712 * sin(b) - 0.0386 * sin(g - f) - 0.0301 * sin(g + b) - 0.0158 * sin(b - e) - 0.0123 * sin(b - g)) * pi /(3600 * 180);
    B_dB = (0.0895 * cos(f + f) - 9.2025 * cos(f) - 0.5736 * cos(h) - 0.0224 * cos(h + c) - 0.0977 * cos(g) - 0.0200 * cos(g - f) - 0.0129 * cos(g + b)) * pi /(3600 * 180);
    de = 0.4090928042223 - B_dB;
    
    double cos_de = cos(de);
    double sin_fd = sin(fd);
    double p = ed * sin_E / (1 - ed * cos_E);
    double q = sin_fd / sin_E;

    C = -0.00009936508495 * (p * sin_ud + q * cos_ud) * cos_de;
    D = -0.00009936508495 * (q * sin_ud - p * cos_ud);
}

void StarMap::Precession(int jd, double &ra, double &dec) {
    double del_ra, del_dec;
    double m = 0.00022361721972456;
    double n = 0.00009717173455170;
    double dmdt = 0.000000135414278898;
    double dndt = -0.000000041369151409;
    double sin_ra = sin(ra);
    double cos_ra = cos(ra);
    double tan_dec = tan(dec);
    double sin_2ra = sin(2 * ra);
    double cos_2ra = cos(2 * ra);
	double del_ra1 = m + n * sin_ra * tan_dec;
	double del_dec1 = n * cos_ra;
	double del_ra2 = dmdt + n * n * sin_2ra / 2 + (dndt * sin_ra + m * n * cos_ra) * tan_dec + n * n * sin_2ra * tan_dec * tan_dec;
	double del_dec2 = dndt * cos_ra - m * n * sin_ra - n * n * sin_ra * sin_ra * tan_dec;
	double del_ra3 = m * n / 2 + 3 * m * n * n * cos_2ra / 2 + 3 * n * dndt * sin_2ra / 2 + ((2 * n * n - m * m + 3 * n * n * cos_2ra) * n * sin_ra + (2 * m * dndt + n * dmdt) * cos_ra) * tan_dec + (3 * n * dndt * sin_2ra + 3 * m * n * n * cos_2ra) * tan_dec * tan_dec + 2 * n * n * sin_ra * (1 + 2 * cos_2ra) * tan_dec;
	double del_dec3 = -m * m * n * cos_ra - n * dmdt * sin_ra - n * n * sin_ra * sin_ra * cos_ra * (1 + 3 * n * tan_dec * tan_dec) - (3 * m * n * n * sin_2ra / 2 + 3 * n * dndt * sin_ra * sin_ra) * tan_dec;
	double del_t = jd / (365.25 * 100);
	del_ra = del_t * del_ra1 + del_t * del_t * del_ra2 / 2 + del_t * del_t * del_t * del_ra3 / 6;
	del_dec = del_t * del_dec1 + del_t * del_t * del_dec2 / 2 + del_t * del_t * del_t * del_dec3 / 6;

	ra += del_ra;
    dec += del_dec;
}

void StarMap::ProperMotion(int jd, double u_ra, double u_dec, double &ra, double &dec) {
    double del_ra;
    double del_dec;
    double del_t = jd / 365.25;
    del_ra = u_ra * del_t * 15 * pi / (3600 * 180);
    del_dec = u_dec * del_t * pi / (3600 * 180);

    ra += del_ra;
    dec += del_dec;
}

int StarMap::JulianDay() {
    int jd, day1, day2, day3, day11, day22;
    day1 = int(365.25 * (year - 2000) + 0.76);
    day11 = int(365.25 * (year - 2000));
    if (day1 == day11) {
        day22 = 60;
    } else {
        day22 = 59;
    }

    switch(month) {
		case 1: day2 = 0; break;
		case 2: day2 = 31; break;
		case 3: day2 = day22; break;
		case 4: day2 = day22 + 31; break;
		case 5: day2 = day22 + 61; break;
		case 6: day2 = day22 + 92; break;
		case 7: day2 = day22 + 122; break;
		case 8: day2 = day22 + 153; break;
		case 9: day2 = day22 + 184; break;
		case 10: day2 = day22 + 214; break;
		case 11: day2 = day22 + 244; break;
		case 12: day2 = day22 + 274; break;
		default: break; 
    }

    day3 = day - 1;
    jd = day1 + day2 + day3;

    return jd;
}

int StarMap::CreateStarMap(StarInfo* stars_of_star_map) {
    int star_count = 0;

    star_count = FindStar(stars_of_star_map, axes_longitude, axes_latitude, x_view, y_view);
    if (star_count == 0) {
        return -1;
    }

    star_count = CoordinateTrans(stars_of_star_map, star_count, axes_longitude, axes_latitude, x_view, y_view);
    if (star_count == 0) {
        return -1;
    }

    return star_count;
}

void StarMap::Dichotomy(StarInfo* stars_of_star_map_temp, double m_max_lat, double m_min_lat, double m_axes_longitude, double m_x_view) {
    char numb[6], longitude[10], latitude[11], magitude[4], sptype[2], pmra[7], pmdec[6];
    double latitude_mid = 0, latitude_next = 0, m_longitude = 0, m_latitude = 0, m_magitude = 0, m_pmRA = 0, m_pmDEC = 0;
    long m_numb = 0;
    char m_sptype = 'A';
    unsigned long long result;
    long all_star_count = 258997;
    std::fstream file;
    file.open(star_map_file_name, std::ios::in);

    long low = 1, mid = 0, high = all_star_count - 1;
    while (low <= high) {
        mid = (low + high) / 2;
        file.seekg(47 * (mid - 1) + 6, std::ios::beg);
        file.read(latitude, 11);
        latitude_mid = (atof(latitude)) * 180 / pi;
        file.seekg(47 - 11, std::ios::cur);
        file.read(latitude, 11);
        latitude_next = (atof(latitude)) * 180 / pi;

        if (latitude_next < m_min_lat) {
            low = mid + 2;
        } else if (latitude_mid > m_min_lat) {
            high = mid - 1;
        } else if (latitude_mid == m_min_lat) {
            low = all_star_count;
        } else {
            ++mid;
            low = all_star_count;
        }
    }

    file.seekg(47 * (mid - 1), std::ios::beg);
    file.read(numb, 6); //获得编号
    m_numb = atol(numb);
    file.read(latitude, 11); //获得赤纬
    m_latitude = (atof(latitude)) * 180 / pi;
    file.read(longitude, 10); //获得赤经
    m_longitude = (atof(longitude)) * 180 / (pi * 15);
    file.read(magitude, 4); //获得星等
    m_magitude = atof(magitude);
    file.read(sptype, 2); //获得光谱型
    m_sptype = sptype[0];
    file.read(pmdec, 6); //获得赤纬自行
    m_pmDEC = atof(pmdec);
    file.read(pmra, 7); //获得赤经自行
    m_pmRA = atof(pmra);

    long sign = m_numb;
    while (m_latitude <= m_max_lat && sign <= all_star_count) {
        if (m_magitude <= magmax && m_magitude >= magmin) {
            double m_x_view_temp = sqrt(2) * 2 * asin(sin(m_x_view * pi / (2 * 180)) / cos(m_latitude * pi / 180)) * 180 / pi;
            double m_max_long = m_axes_longitude + (m_x_view_temp * 0.5) / 15; //最大赤经
            double m_min_long = m_axes_longitude - (m_x_view_temp * 0.5) / 15; //最小赤经

            if ((m_longitude >= m_min_long) && (m_longitude <= m_max_long)) {
                stars_of_star_map_temp->SetStarNumb(m_numb);
                stars_of_star_map_temp->SetStarLongitude(m_longitude);
                stars_of_star_map_temp->SetStarLatitude(m_latitude);
                stars_of_star_map_temp->SetStarMagnitude(m_magitude);
                stars_of_star_map_temp->SetStarSptype(m_sptype);
                stars_of_star_map_temp->SetStarProperMotionRA(m_pmRA);
                stars_of_star_map_temp->SetStarProperMotionDEC(m_pmDEC);
                ++stars_of_star_map_temp;
            }
        }
        file.seekg(1, std::ios::cur);
        file.read(numb, 6); //获得编号
        m_numb = atol(numb);
        file.read(latitude, 11); //获得赤纬
        m_latitude = (atof(latitude)) * 180 / pi;
        file.read(longitude, 10); //获得赤经
        m_longitude = (atof(longitude)) * 180 / (pi * 15);
        file.read(magitude, 4); //获得星等
        m_magitude = atof(magitude);
        file.read(sptype, 2); //获得光谱型
        m_sptype = sptype[0];
        file.read(pmdec, 6); //获得赤纬自行
        m_pmDEC = atof(pmdec);
        file.read(pmra, 7); //获得赤经自行
        m_pmRA = atof(pmra);
        sign++;
    }

    file.close();
}

int StarMap::FindStar(StarInfo* stars_of_star_map_temp, double m_axes_longitude, double m_axes_latitude, double m_x_view, double m_y_view) {
    int k = 0;

    // 视场范围
    double max_latitude, min_latitude;
    max_latitude = m_axes_latitude + m_y_view * (sqrt(2) * sin((45 + rotate_angle) * pi / 180) * 0.5); //视场最大赤纬
    min_latitude = m_axes_latitude - m_y_view * (sqrt(2) * sin((45 + rotate_angle) * pi / 180) * 0.5); //视场最小赤纬

    // 二分法搜索星图
    Dichotomy(stars_of_star_map_temp, max_latitude, min_latitude, m_axes_longitude, m_x_view);
    
    while (stars_of_star_map_temp[k].GetStarNumb() != 0) {
        ++k;
    }
    return k;
}

void StarMap::CreateStarImage() {
    int star_count; //视场中星的个数
    julian_day_ = JulianDay();
    
//    unsigned char* lpNewDIBBits;
    cv::Mat lpNewDIBBits;
    if (is_color) {
//       lpNewDIBBits = new unsigned char[width * height * 3];
//        cv::Mat img_temp(width, height, CV_8UC3, cv::Scalar(0, 0, 0));
//        lpNewDIBBits = img_temp.clone();
        lpNewDIBBits.create(width, height, CV_8UC3);
    } else {
//        lpNewDIBBits = new unsigned char[width * height];
//        cv::Mat img_temp(width, height, CV_8UC1, cv::Scalar(0));        
//        lpNewDIBBits = img_temp.clone();
        lpNewDIBBits.create(width, height, CV_8UC1);
    }

    StarInfo stars_of_star_map[5000];
    star_count = CreateStarMap(stars_of_star_map);
    if (star_count == 0) {
        return;
    }

    if (ShowStarMap(lpNewDIBBits, stars_of_star_map, star_count, x_view)) {
        SavePic(lpNewDIBBits);
        SaveStarCountFile();
    }
}

bool StarMap::ShowStarMap(cv::Mat& lpNewDIBBits, StarInfo *stars_of_star_map, int size, double m_x_view)
{
	int l_width = width;
	int l_height = height;
	int x_orig, y_orig;
	int err_axes_x = 0, err_axes_y = 0; //光轴指向误差
	int x, y;
	double mag;
	char sptype;
	int i;
	int background;

	//指向像素的指针
//	unsigned char* lpDst;

	//图像每行的字节数
//	long lLineBytes;

	//显示星点
	srand((unsigned)time(NULL)); //模拟光轴指向误差，产生0-246的随机数
	int err = int(l_width * (axes_error / m_x_view));
	
	if (err != 0) {
		err_axes_x = rand() % err;
		err_axes_y = rand() % err;
	}

	if (is_color) {
//		lLineBytes=(((l_width * 24) + 31) / 32) * 4;

        if (background_value == 0) {
			//将图像像素全部变为0
			for (i = 0; i < l_height; ++i)
			{
				for (int j = 0; j < l_width; ++j) {
//					lpDst = unsigned char*(lpNewDIBBits + lLineBytes * (lHeight - 1 - i) + 3 * j);
//					*lpDst = 0; //blue
//					*(lpDst + 1) = 0; //green
//					*(lpDst + 2) = 0; //red
                    lpNewDIBBits.at<cv::Vec3b>(i, j)[0] = 0;
                    lpNewDIBBits.at<cv::Vec3b>(i, j)[1] = 0;
                    lpNewDIBBits.at<cv::Vec3b>(i, j)[2] = 0;
				}
			}
		} else {
            for (i = 0; i < l_height; ++i) {
				for(int j = 0; j < l_width; ++j) {
                    background = rand() % background_value;
//					lpDst = unsigned char*(lpNewDIBBits + lLineBytes * (lHeight - 1 - i) + 3 * j);
//					*lpDst = backgroud; //blue
//					*(lpDst + 1) = backgroud; //green
//					*(lpDst + 2) = backgroud; //red
					lpNewDIBBits.at<cv::Vec3b>(i, j)[0] = background;
					lpNewDIBBits.at<cv::Vec3b>(i, j)[1] = background;
					lpNewDIBBits.at<cv::Vec3b>(i, j)[2] = background;
				}
			}
		}

	   for (i = 0; i < size; ++i){
		  //星表中星的位置
		  x_orig = stars_of_star_map[i].GetCoordinateX();
		  y_orig = stars_of_star_map[i].GetCoordinateY();
		

	      //图像上星的位置（加入光轴及恒星位置误差）
		  x = int(x_orig - err / 2 + err_axes_x);
		  y = int(y_orig - err / 2 + err_axes_y);

		  mag = stars_of_star_map[i].GetStarMagnitude();

		  int starpix, ra ;
		  double zoom = 0.5;
		  starpix = int((12 - mag) * zoom + 2);
		  ra = starpix - 3;
	      double uu = 0.466 * ra;
		  sptype = stars_of_star_map[i].GetStarSptype();
          double red, green, blue;
          switch(sptype) //对恒星颜色作判断
		  {
                 case 'O': red=0.55, green=0.59, blue=1.00; break;
                 case 'B': red=0.47, green=0.70, blue=0.90; break;
		         case 'A': red=1.00, green=1.00, blue=1.00; break;
		         case 'F': red=1.00, green=1.00, blue=0.67; break;
	             case 'G': red=0.80, green=0.76, blue=0.20; break;
		         case 'K': red=0.80, green=0.41, blue=0.02; break;
		         case 'M': red=0.75, green=0.50, blue=0.20; break;
		         default : red=1.00, green=1.00, blue=1.00; break;
								
		  }
         double pex0 = 4500 * (pow(zoom, 2.5)) * (pow(1.8, (9 - mag)));
		 if ((x >= 0) && (x < l_height) && (y >= 0) && (y < l_width)) {
			StarCount(mag);
//			lpDst = unsigned char*(lpNewDIBBits + lLineBytes * (l_height - 1 - x) + 3 * y);
			for (int j = -starpix; j <= starpix; ++j) {
				if ((x + j) >= 0 && (x + j) < l_height) {
					for (int k = -starpix; k <= starpix; ++k) {
						if ((y + k) >= 0 && (y + k) < l_width) {
							double a,b,c;        //对数变换的系数          为显示需要，将计算得到的灰度值进行对比度变换
							a = 0.4;
							b = 32;
							c = 55;
								
							//斑点像素  假设服从高斯分布
							double pex;
				            float dx, dy;             
							dx = float(abs(j)); //X方向到质心距离
							dy = float(abs(k)); //Y方向到质心距离
							pex = (pex0 * exp(-(pow(dx, 2) + pow(dy, 2)) / (2 * pow(uu, 2))) / (2 * pi * pow(uu, 2)));
							if (pex > 255) {
								pex = 255;
							}
								
							double color[3];
							int colornum;
//							color[0] = (pex * blue + *(lpDst - lLineBytes * j + 3 * k + 2));
//	                        color[1] = (pex * green + *(lpDst - lLineBytes * j + 3 * k + 1));
//							color[2] = (pex * red + *(lpDst - lLineBytes * j + 3 * k));
							color[0] = pex * blue + lpNewDIBBits.at<cv::Vec3b>(x + j, y + k)[0];
							color[1] = pex * green + lpNewDIBBits.at<cv::Vec3b>(x + j, y + k)[1];
							color[2] = pex * red + lpNewDIBBits.at<cv::Vec3b>(x + j, y + k)[2];
							for (colornum = 0; colornum < 3; ++colornum) {
                            	if (color[colornum] > 255) {
									color[colornum] = 255;
                                }
//                                *(lpDst - lLineBytes * j + 3 * k + colornum) = unsigned char(color[colornum]);
                                lpNewDIBBits.at<cv::Vec3b>(x + j, y + k)[colornum] = color[colornum];
							}
						}
					}
				}
			}
		}

	   }

	} else {
//		lLineBytes = (((l_width * 8) + 31) / 32) * 4;

		if (background_value==0) {
			//将图像像素全部变为0
			for (i = 0; i < l_height; ++i) {
				for(int j = 0; j < l_width; ++j) {
//					lpDst = unsigned char*(lpNewDIBBits + lLineBytes * (l_height - 1 - i) + j);
//					*lpDst = 0;
                    lpNewDIBBits.at<uchar>(i, j) = 0;					
				}
			}
		} else {
            for (i = 0; i < l_height ; ++i) {
				for (int j = 0; j < l_width; ++j) {
                    background = rand() % background_value;
//					lpDst = unsigned char*(lpNewDIBBits + lLineBytes * (l_height - 1 - i) + j);
//					*lpDst = backgroud;
					lpNewDIBBits.at<uchar>(i, j) = background;
				}
			}
		}
	
	   for (i = 0; i < size; ++i) {
		  //星表中星的位置
		  x_orig = stars_of_star_map[i].GetCoordinateX();
		  y_orig = stars_of_star_map[i].GetCoordinateY();
		

	      //图像上星的位置（加入光轴及恒星位置误差）
		  x = int(x_orig - err / 2 + err_axes_x);
		  y = int(y_orig - err / 2 + err_axes_y);

		  mag = stars_of_star_map[i].GetStarMagnitude();
		  int starpix, ra;
		  
		  double zoom = 0.5;
		  starpix = int((12 - mag) * zoom + 2);
		  ra = starpix - 3;
	      double uu = 0.466 * ra;
          double pex0 = 4500 * (pow(zoom, 2.5)) * (pow(1.8, (9 - mag)));
		if ((x >= 0) && (x < l_height) && (y >= 0) && (y < l_width)) {
			StarCount(mag);
//			 lpDst = unsigned char*(lpNewDIBBits + lLineBytes * (l_height - 1 - x) + y);
				for (int j = -starpix; j <= starpix; ++j) {
					if ((x + j) >= 0 && (x + j) < l_height) {
						for (int k = -starpix; k <= starpix; ++k) {
							if ((y + k) >= 0 && (y + k) < l_width) {
								double a, b, c;        //对数变换的系数          为显示需要，将计算得到的灰度值进行对比度变换
								a = 0.4;
								b = 32;
								c = 55;
								
								//斑点像素  假设服从高斯分布
								double pex;
								float dx, dy;             
								dx = float(abs(j));   //X方向到质心距离
								dy = float(abs(k));   //Y方向到质心距离
								pex = (pex0 * exp(-(pow(dx, 2) + pow(dy, 2)) / (2 *pow(uu, 2))) / (2 * pi * pow(uu, 2)));
								if (pex > 255) {
									pex = 255;
								}
								if (pex < 0) {
                                   pex = 0;
								}
								
								double color;
//								 color = (*(lpDst - lLineBytes * j + k) + pex);
								color = pex + lpNewDIBBits.at<uchar>(x + j, y + k);
                                       if (color > 255)
                                       {
									            color = 255;
                                       }
//                                       *(lpDst - lLineBytes * j + k) = unsigned char(color);
                                       lpNewDIBBits.at<uchar>(x + j, y + k) = color;
				           }
						}
					}
				}

		}

	   }

	}
	
	
	return true;
}

void StarMap::SavePic(cv::Mat& lpNewDIBBits)
{
    cv::imwrite(save_file_name, lpNewDIBBits);
    /*
   if (is_color)
   {
	PBITMAPINFO pBmpInfor = m_pBMPinfo;	
	pBmpInfor = (PBITMAPINFO) new BYTE[sizeof(BITMAPINFOHEADER)+ 256 * sizeof(RGBQUAD)];
	ASSERT(pBmpInfor != NULL);
	//文件头
	PBITMAPINFOHEADER pBMPInfoHeaderr = &pBmpInfor->bmiHeader;
	pBMPInfoHeaderr->biSize = sizeof(BITMAPINFOHEADER);
	pBMPInfoHeaderr->biBitCount =24;
	pBMPInfoHeaderr->biHeight = height;
	pBMPInfoHeaderr->biWidth = width;
	pBMPInfoHeaderr->biCompression = BI_RGB;
	pBMPInfoHeaderr->biPlanes = 1;
	pBMPInfoHeaderr->biSizeImage = (long)(width*height*3);
	pBMPInfoHeaderr->biXPelsPerMeter = pBMPInfoHeaderr->biYPelsPerMeter = 0;
	pBMPInfoHeaderr->biClrUsed = 0;
	pBMPInfoHeaderr->biClrImportant = 0;

	m_pBMPinfo = pBmpInfor;
	pImgData = lpNewDIBBits;
	
	CFile filesave;

	BITMAPFILEHEADER BmpFileHeaders;
	BmpFileHeaders.bfType = 0x4D42;
	BmpFileHeaders.bfOffBits = sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER)
		+ 256 * sizeof(RGBQUAD);
	BmpFileHeaders.bfSize = BmpFileHeaders.bfOffBits 
		+ this->m_pBMPinfo->bmiHeader.biWidth*this->m_pBMPinfo->bmiHeader.biHeight
		*this->m_pBMPinfo->bmiHeader.biBitCount/8;
	
	filesave.Open(pSaveFileName, CFile::modeCreate|CFile::modeWrite   );
	filesave.Write(&BmpFileHeaders, sizeof(BmpFileHeaders));
	filesave.Write(this->m_pBMPinfo, sizeof(BITMAPINFOHEADER)+ 256*sizeof(RGBQUAD));
	filesave.WriteHuge(this->pImgData, width*height*3);
	filesave.Close();	
   }
   else
   {
	PBITMAPINFO pBmpInfor = m_pBMPinfo;	
	pBmpInfor = (PBITMAPINFO) new BYTE[sizeof(BITMAPINFOHEADER)+ 256 * sizeof(RGBQUAD)];
	ASSERT(pBmpInfor != NULL);
	//文件头
	PBITMAPINFOHEADER pBMPInfoHeaderr = &pBmpInfor->bmiHeader;
	pBMPInfoHeaderr->biSize = sizeof(BITMAPINFOHEADER);
	pBMPInfoHeaderr->biBitCount =8;
	pBMPInfoHeaderr->biHeight = height;
	pBMPInfoHeaderr->biWidth = width;
	pBMPInfoHeaderr->biCompression = BI_RGB;
	pBMPInfoHeaderr->biPlanes = 1;
	pBMPInfoHeaderr->biSizeImage = 0;
	pBMPInfoHeaderr->biXPelsPerMeter = pBMPInfoHeaderr->biYPelsPerMeter = 0;
    pBMPInfoHeaderr->biClrUsed = 256;
	pBMPInfoHeaderr->biClrImportant = 0;
	//调色板
	RGBQUAD *pGrayr = pBmpInfor->bmiColors;
	int i;
	for(i=0; i<256; i++)
	{
	pGrayr[i].rgbBlue = (BYTE)i;
	pGrayr[i].rgbGreen = (BYTE)i;
	pGrayr[i].rgbRed = (BYTE)i;
	pGrayr[i].rgbReserved = 0;
	}
		
	
	m_pBMPinfo = pBmpInfor;
	pImgData = lpNewDIBBits;
	
	CFile filesave;
    
	BITMAPFILEHEADER BmpFileHeaders;
	BmpFileHeaders.bfType = 0x4D42;
	BmpFileHeaders.bfOffBits = sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER)
		+ 256 * sizeof(RGBQUAD);
	BmpFileHeaders.bfSize = BmpFileHeaders.bfOffBits 
		+ this->m_pBMPinfo->bmiHeader.biWidth*this->m_pBMPinfo->bmiHeader.biHeight
		*this->m_pBMPinfo->bmiHeader.biBitCount/8;
	
	filesave.Open(pSaveFileName, CFile::modeCreate|CFile::modeWrite   );
	filesave.Write(&BmpFileHeaders, sizeof(BmpFileHeaders));
	filesave.Write(this->m_pBMPinfo, sizeof(BITMAPINFOHEADER)+ 256*sizeof(RGBQUAD));
	filesave.WriteHuge(this->pImgData, width*height );
	filesave.Close();	

   }
*/
}
