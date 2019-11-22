#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <math.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "DUtils/Random.h"

Eigen::Matrix3d getR(float phi, float omiga, float kapa){
    Eigen::Matrix3d R;
    R<<cos(phi)*cos(kapa)-sin(phi)*sin(omiga)*sin(kapa),-cos(phi)*sin(kapa)-sin(phi)*sin(omiga)*cos(kapa),-sin(phi)*cos(omiga),
        cos(omiga)*sin(kapa),cos(omiga)*cos(kapa),-sin(omiga),
        sin(phi)*cos(kapa)+cos(phi)*sin(omiga)*sin(kapa),-sin(phi)*sin(kapa)+cos(phi)*sin(omiga)*cos(kapa),cos(phi)*cos(omiga);

    return R;
}