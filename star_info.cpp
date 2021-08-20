/*************************************************************************
    > File Name: star_info.cpp
    > Author: WangChenxu
    > Mail: wangchenxu21@outlook.com
    > Created Time: Fri 20 Aug 2021 02:56:32 AM UTC
 ************************************************************************/

#include "star_info.h"

StarInfo::StarInfo() {
    star_numb = 0;
    star_longitude = 0.;
    star_latitude = 0.;
    star_magnitude = -10.;
    star_sptype = 'A';
    star_pmRA = 0.;
    star_pmDEC = 0.;
    coordinate_x = 0;
    coordinate_y = 0;
}

StarInfo::~StarInfo() {}

void StarInfo::SetStarMagnitude(double m_star_magnitude) {
    star_magnitude = m_star_magnitude;
}

void StarInfo::SetStarLatitude(double m_star_latitude) {
    star_latitude = m_star_latitude;
}

void StarInfo::SetStarLongitude(double m_star_longitude) {
    star_longitude = m_star_longitude;
}

void StarInfo::SetStarNumb(long m_star_numb) {
    star_numb = m_star_numb;
}

void StarInfo::SetStarSptype(char m_star_sptype) {
    star_sptype = m_star_sptype;
}

void StarInfo::SetStarProperMotionRA(double m_star_pmRA) {
    star_pmRA = m_star_pmRA;
}

void StarInfo::SetStarProperMotionDEC(double m_star_pmDEC) {
    star_pmDEC = m_star_pmDEC;
}

void StarInfo::SetCoordinateY(int m_coordinate_y) {
    coordinate_y = m_coordinate_y;
}

void StarInfo::SetCoordinateX(int m_coordinate_x) {
    coordinate_x = m_coordinate_x;
}

int StarInfo::GetCoordinateY() {
    return coordinate_y;
}

int StarInfo::GetCoordinateX() {
    return coordinate_x;
}

double StarInfo::GetStarMagnitude() {
    return star_magnitude;
}

double StarInfo::GetStarLatitude() {
    return star_latitude;
}

double StarInfo::GetStarLongitude() {
    return star_longitude;
}

long StarInfo::GetStarNumb() {
    return star_numb;
}

char StarInfo::GetStarSptype() {
    return star_sptype;
}

double StarInfo::GetStarProperMotionRA() {
    return star_pmRA;
}

double StarInfo::GetStarProperMotionDEC() {
    return star_pmDEC;
}
