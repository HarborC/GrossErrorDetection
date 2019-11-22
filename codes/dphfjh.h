// 单片后方交会
#include "Common.h"

using namespace std;
using namespace Eigen;

void dphfjh(vector<unordered_map<int, Eigen::Vector2d>> images, unordered_map<int, Eigen::Vector3d> control_pt, int N)
{
    // 单片后方交会
    // 注意控制点是左手坐标系，所以X，Y要调换
    vector<int> cpt;
    Eigen::Vector3d mean_pt = Eigen::Vector3d(0,0,0);
    for(unordered_map<int, Eigen::Vector3d>::iterator it = control_pt.begin(); it != control_pt.end(); it++){
        if(images[2].find(it->first) != images[2].end()){
            cpt.push_back(it->first);
            mean_pt += it->second;
        }
    }
    mean_pt /= cpt.size();

    float phi3 = 0;
    float omiga3 = 0;
    float kapa3 = 0;
    float Xs3 = mean_pt(1);
    float Ys3 = mean_pt(0);
    float Zs3 = 2660;
    float f = 100.5;

    for(int it = 0; it<N; it++){
        Eigen::Matrix<double, 3, 3> R;
        R = getR(phi3, omiga3, kapa3);

        Eigen::Matrix<double, Dynamic, 6> A;
        Eigen::Matrix<double, Dynamic, 1> L;
        Eigen::Matrix<double, 6, 1> delta;
        A.resize(cpt.size()*2, 6);
        L.resize(cpt.size()*2, 1);

        for(int i = 0; i<cpt.size(); i++){
            float x = images[2][cpt[i]](0);
            float y = images[2][cpt[i]](1);
            float xw = control_pt[cpt[i]](1), yw = control_pt[cpt[i]](0), zw = control_pt[cpt[i]](2);
            
            Eigen::Vector3d Xw;
            Xw<<xw-Xs3, yw-Ys3, zw-Zs3;
            Eigen::Vector3d Xw_ = R.transpose()*Xw;

            float x_ = -f*Xw_(0)/Xw_(2);
            float y_ = -f*Xw_(1)/Xw_(2);

            A(2*i+0, 0) = y_*sin(omiga3) - (x_*(x_*cos(kapa3)-y_*sin(kapa3))/f+f*cos(kapa3))*cos(omiga3);
            A(2*i+0, 1) = -f*sin(kapa3) - x_*(x_*sin(kapa3)+y_*cos(kapa3))/f;
            A(2*i+0, 2) = y_;
            A(2*i+0, 3) = (R(0,0)*f+R(0,2)*x_)/Xw_(2);
            A(2*i+0, 4) = (R(1,0)*f+R(1,2)*x_)/Xw_(2);
            A(2*i+0, 5) = (R(2,0)*f+R(2,2)*x_)/Xw_(2);

            A(2*i+1, 0) = -x_*sin(omiga3) - (y_*(x_*cos(kapa3)-y_*sin(kapa3))/f-f*sin(kapa3))*cos(omiga3);
            A(2*i+1, 1) = -f*cos(kapa3) - y_*(x_*sin(kapa3)+y_*cos(kapa3))/f;
            A(2*i+1, 2) = -x_;
            A(2*i+1, 3) = (R(0,1)*f+R(0,2)*y_)/Xw_(2);
            A(2*i+1, 4) = (R(1,1)*f+R(1,2)*y_)/Xw_(2);
            A(2*i+1, 5) = (R(2,1)*f+R(2,2)*y_)/Xw_(2);

            L(2*i+0, 0) = x - x_;
            L(2*i+1, 0) = y - y_;
        }

        delta = (A.transpose()*A).inverse()*(A.transpose()*L);
        phi3 += delta(0);
        omiga3 += delta(1);
        kapa3 += delta(2);
        Xs3 += delta(3);
        Ys3 += delta(4);
        Zs3 += delta(5);

        cerr<<"phi, omiga, kapa, Xs3, Ys3, Zs3 :"<<delta.transpose()<<endl;
    }
    cerr<<endl;
    cerr<<"phi, omiga, kapa, Xs3, Ys3, Zs3 :"<<phi3<<" "<<omiga3<<" "<<kapa3<<" "<<Xs3<<" "<<Ys3<<" "<<Zs3<<endl;

}