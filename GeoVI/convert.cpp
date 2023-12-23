#include "convert.h"
#include <proj/crs.hpp>
#include <proj/io.hpp>
#include <proj/util.hpp>
#include <proj.h>
#include <vector>
#include <cmath>
#include <sstream>
#include <ctime>
#include <iomanip>





using namespace geovi;

static PJ* P_GIS_WGS84_TO_UTM;
static PJ* P_GIS_UTM_TO_WGS84;
static PJ* P_J;

int utm_num(double latitude,double longitude){
    return floor(floor(longitude)/6) + 31;
}

void init_pj(int64_t utm_epsg_num){
    std::string epsg_str = std::string("EPSG:").append(std::to_string(utm_epsg_num));
    P_GIS_WGS84_TO_UTM = proj_create_crs_to_crs(PJ_DEFAULT_CTX,"EPSG:4326",epsg_str.c_str(),NULL);
    if( P_GIS_WGS84_TO_UTM == 0 ){
        std::cout << "init converter <wgs84 to utm> failed" << std::endl;
    }else{
        PJ* P_for_GIS = proj_normalize_for_visualization(PJ_DEFAULT_CTX,P_GIS_WGS84_TO_UTM);
        if(P_for_GIS == 0){
            std::cout << "normalize converter < wgs84 to utm failed> " << std::endl;
            proj_destroy(P_GIS_WGS84_TO_UTM);
        }else {
            proj_destroy(P_GIS_WGS84_TO_UTM);
            P_GIS_WGS84_TO_UTM = P_for_GIS;
        }
    }

    P_GIS_UTM_TO_WGS84 = proj_create_crs_to_crs(PJ_DEFAULT_CTX,epsg_str.c_str(),"EPSG:4326",NULL);
    if(P_GIS_UTM_TO_WGS84 == 0){
        std::cout << "init converter < utm to wgs84 > failed" << std::endl;
    }else{
        PJ* P_for_GIS = proj_normalize_for_visualization(PJ_DEFAULT_CTX,P_GIS_UTM_TO_WGS84);
        if(P_for_GIS == 0){
            std::cout << "normalize converter < utm to wgs84 failed> " << std::endl;
            proj_destroy(P_GIS_UTM_TO_WGS84);
        }else {
            proj_destroy(P_GIS_UTM_TO_WGS84);
            P_GIS_UTM_TO_WGS84 = P_for_GIS;
        }
    }

}


void proj_convert_wgs84_to_utm(double lat,double lon,Point2& p){
    PJ_COORD c,c_out;
    c.lpzt.lam = lon;
    c.lpzt.phi = lat;
    c.lpzt.z = 0.0;
    c.lpzt.t = HUGE_VAL;
    c_out = proj_trans(P_GIS_WGS84_TO_UTM,PJ_FWD,c);
    p.x = c_out.xy.x;
    p.y = c_out.xy.y;
}

void  proj_convert_utm_to_wgs84(double x,double y,Point2& p){
    PJ_COORD  c,c_out;
    c.xy.x = x;
    c.xy.y = y;
    c_out = proj_trans(P_GIS_UTM_TO_WGS84,PJ_FWD,c);
    p.y = c.lp.lam;
    p.x = c.lp.phi;
}

CoordinateSystemConverter::CoordinateSystemConverter() {
    init_pj(LongitudeBands::band_14);
}

CoordinateSystemConverter::CoordinateSystemConverter(geovi::LongitudeBands bd):band(bd) {
    init_pj(bd);
}

void CoordinateSystemConverter::convert(CoordinateSystemType srcCRS,CoordinateSystemType targetCRS,Point2& point){
    if( srcCRS == CoordinateSystemType::WGS84 && targetCRS == CoordinateSystemType::UTM){
        convertBetweenWGS84AndUTM(point);
    }
    if( srcCRS == CoordinateSystemType::UTM && targetCRS == CoordinateSystemType::WGS84){
        convertBetweenUTMAndWGS84(point);
    }
}

void CoordinateSystemConverter::convertBetweenUTMAndWGS84(Point2& point){
    proj_convert_utm_to_wgs84(point.x,point.y,point);



}

void CoordinateSystemConverter::convertBetweenWGS84AndUTM(Point2& point){
    proj_convert_wgs84_to_utm(point.x,point.y,point);

}

LongitudeBands CoordinateSystemConverter::utm_identity(double longitude, double latitude) {
    return LongitudeBands(utm_num(latitude,longitude));
}

std::chrono::system_clock::time_point StringAndTimeConverter::convertStringToTime(const std::string& timeStr,const char* fmt){
    std::tm tm = {};
    std::istringstream ss(timeStr);
    ss >> std::get_time(&tm,fmt);
    std::time_t tt =std::mktime(&tm);
    return std::chrono::system_clock::from_time_t(tt);
}

std::string StringAndTimeConverter::formatTime(const std::chrono::system_clock::time_point& timePoint,const char* fmt){
    std::time_t tt = std::chrono::system_clock::to_time_t(timePoint);
    std::tm* tm = std::localtime(&tt);
    std::ostringstream oss;
    oss << std::put_time(tm,fmt);
    return oss.str();
}

unsigned long long TimeUtilHelper::get_current_millis() {
    auto now = std::chrono::system_clock::now();
    auto duration = now.time_since_epoch();
    return std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
}

double StatisticUtilHelper::mean(const std::vector<int> &data) {
    int sum = 0;
    for(auto d : data){
        sum += d;
    }
    return sum/data.size();
}

double StatisticUtilHelper::mean(const std::vector<double> &data) {
    double sum = 0.0;
    for(auto d : data){
        sum += d;
    }
    return sum/data.size();
}

double StatisticUtilHelper::variance(const std::vector<int> &data) {
    double mean = StatisticUtilHelper::mean(data);
    double sum_squared_differences = 0.0;
    for(auto d: data){
        double difference = d - mean;
        sum_squared_differences += std::pow(difference,2);
    }
    return sum_squared_differences/data.size();
}

double StatisticUtilHelper::variance(const std::vector<double> &data) {
    double mean = StatisticUtilHelper::mean(data);
    double sum_squared_differences = 0.0;
    for(auto d: data){
        double difference = d - mean;
        sum_squared_differences += std::pow(difference,2);
    }
    return sum_squared_differences/data.size();
}

int StatisticUtilHelper::max_value(const std::vector<int> &data) {
    int max = std::numeric_limits<int>::min();
    for(auto d : data){
        max = d > max ? d : max ;
    }
    return max;
}

double StatisticUtilHelper::max_value(const std::vector<double> &data) {
    double max = std::numeric_limits<int>::min();
    for(auto d : data){
        max = d > max ? d : max ;
    }
    return max;
}