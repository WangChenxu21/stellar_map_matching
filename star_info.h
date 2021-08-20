/*************************************************************************
    > File Name: star_info.h
    > Author: WangChenxu
    > Mail: wangchenxu21@outlook.com
    > Created Time: Fri 20 Aug 2021 02:34:22 AM UTC
 ************************************************************************/

#ifndef _STARINFO_H_
#define _STARINFO_H_

class StarInfo {
public:
	StarInfo();
	virtual ~StarInfo();

	void SetStarMagnitude(double m_star_magnitude);
	void SetStarLatitude(double m_star_latitude);
	void SetStarLongitude(double m_star_longitude);
	void SetStarNumb(long m_star_numb);
	void SetStarSptype(char m_star_sptype);
	void SetStarProperMotionRA(double m_star_pmRA);
	void SetStarProperMotionDEC(double m_star_pmDEC);
	void SetCoordinateY(int m_coordinate_y);
	void SetCoordinateX(int m_coordinate_x);

	int GetCoordinateY();
	int GetCoordinateX();
	double GetStarMagnitude();
	double GetStarLatitude();
	double GetStarLongitude();
	long GetStarNumb();
	char GetStarSptype();
	double GetStarProperMotionRA();
	double GetStarProperMotionDEC();

private:
	int coordinate_y; //星点在平面的Y坐标
	int coordinate_x; //星点在平面的X坐标
	double star_magnitude; //星等
	double star_latitude; //赤纬
	double star_longitude; //赤经
	char star_sptype; //光谱型
	long star_numb; //编号
	double star_pmRA; //赤经自行
	double star_pmDEC; //赤纬自行
};

#endif
