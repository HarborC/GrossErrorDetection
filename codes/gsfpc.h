// 光束法区域网平差
#include "Common.h"

void gsfpc(vector<unordered_map<int, Eigen::Vector2d>> images, vector<int> idx12, vector<int> idx13, unordered_map<int, vector<int>> same_pt, int N){
    
    // 光束法区域网平差
    float f = 100.5;
    float phi[3] = {0.0, 0.0, 0.0};
    float omiga[3] = {0.0, 0.0, 0.0};
    float kapa[3] = {0.0, 0.0, 0.0};
    float Xs[3] = {0.0, 0.0, 0.0};
    float Ys[3] = {0.0, 0.0, 0.0};
    float Zs[3] = {0.0, 0.0, 0.0};

    for(int j = 0; j<idx12.size(); j++){
        Xs[1] = Xs[1] + images[0][idx12[j]](0)-images[1][idx12[j]](0);
    }
    Xs[1] = Xs[1]/idx12.size();

    for(int j = 0; j<idx13.size(); j++){
        Xs[2] = Xs[2] + images[0][idx13[j]](0)-images[2][idx13[j]](0);
    }
    Xs[2] = Xs[2]/idx13.size();

    std::vector<Eigen::Vector3d> Ps(same_pt.size());
    int nn = 0;
    for(unordered_map<int, vector<int>>::iterator it = same_pt.begin(); it != same_pt.end(); it++, nn++){
        int Pid = it->first;
        vector<int> sets = it->second;

        int cam1 = sets[0]; int cam2 = sets[1];

        Eigen::Vector2d x1 = images[cam1][Pid]; Eigen::Vector2d x2 = images[cam2][Pid];
        float Bx = Xs[cam2] - Xs[cam1];
        float Bz = Zs[cam2] - Zs[cam1];
        float N1 = (Bx*(-f) - Bz*x2(0))/(x1(0)*(-f) - x2(0)*(-f));
        float N2 = (Bx*(-f) - Bz*x1(0))/(x1(0)*(-f) - x2(0)*(-f));

        Ps[nn]=Eigen::Vector3d(Xs[cam1]+N1*x1(0),Ys[cam1]+N1*x1(1),Zs[cam1]+N1*(-f));
    }
    Eigen::Matrix<double, Dynamic, Dynamic> A;
    Eigen::Matrix<double, Dynamic, 1> L;
    A.setZero(nl*2, 11 + same_pt.size()*3);
    L.setZero(nl*2, 1);
    Eigen::Matrix<double, Dynamic, Dynamic> P;
    P.setZero(nl*2, nl*2);
    for(int i = 0; i<nl*2; i++){
        P(i,i) = 1.0;
    }

    Eigen::Matrix<double, Dynamic, 1> delta;
    delta.resize(11 + same_pt.size()*3,1);
    Eigen::Matrix<double, Dynamic, Dynamic> H;
    H.resize(11 + same_pt.size()*3,11 + same_pt.size()*3);
    Eigen::Matrix<double, Dynamic, 1> v;
    v.resize(nl*2,1);

    for(int itime = 0; itime<N; itime++){
        
        std::vector<Eigen::Matrix3d> R(3);
        for(int i = 0; i<3; i++){
            R[i] = getR(phi[i], omiga[i], kapa[i]);
        }

        A.setZero(nl*2, 11 + same_pt.size()*3);
        L.setZero(nl*2, 1);

        int npixel = 0;
        int nP = 0;
        for(unordered_map<int, vector<int>>::iterator it = same_pt.begin(); it != same_pt.end(); it++, nP++){
            int Pid = it->first;
            vector<int> sets = it->second;
            for(int ic = 0; ic<sets.size(); ic++){
                int camid = sets[ic];
                float x = images[camid][Pid](0);
                float y = images[camid][Pid](1);
                float xw = Ps[nP](0), yw = Ps[nP](1), zw = Ps[nP](2);
                
                Eigen::Vector3d Xw;
                Xw<<xw-Xs[camid], yw-Ys[camid], zw-Zs[camid];
                Eigen::Matrix3d Rc = R[camid];
                Eigen::Vector3d Xw_ = Rc.transpose()*Xw;

                float x_ = -f*Xw_(0)/Xw_(2);
                float y_ = -f*Xw_(1)/Xw_(2);

                if(camid == 2){
                    A(2*npixel+0, 5) = y_*sin(omiga[camid]) - (x_*(x_*cos(kapa[camid])-y_*sin(kapa[camid]))/f+f*cos(kapa[camid]))*cos(omiga[camid]);
                    A(2*npixel+0, 6) = -f*sin(kapa[camid]) - x_*(x_*sin(kapa[camid])+y_*cos(kapa[camid]))/f;
                    A(2*npixel+0, 7) = y_;
                    A(2*npixel+0, 8) = (Rc(0,0)*f+Rc(0,2)*x_)/Xw_(2);
                    A(2*npixel+0, 9) = (Rc(1,0)*f+Rc(1,2)*x_)/Xw_(2);
                    A(2*npixel+0, 10) = (Rc(2,0)*f+Rc(2,2)*x_)/Xw_(2);

                    A(2*npixel+1, 5) = -x_*sin(omiga[camid]) - (y_*(x_*cos(kapa[camid])-y_*sin(kapa[camid]))/f-f*sin(kapa[camid]))*cos(omiga[camid]);
                    A(2*npixel+1, 6) = -f*cos(kapa[camid]) - y_*(x_*sin(kapa[camid])+y_*cos(kapa[camid]))/f;
                    A(2*npixel+1, 7) = -x_;
                    A(2*npixel+1, 8) = (Rc(0,1)*f+Rc(0,2)*y_)/Xw_(2);
                    A(2*npixel+1, 9) = (Rc(1,1)*f+Rc(1,2)*y_)/Xw_(2);
                    A(2*npixel+1, 10) = (Rc(2,1)*f+Rc(2,2)*y_)/Xw_(2);
                }else if(camid == 1){
                    A(2*npixel+0, 0) = y_*sin(omiga[camid]) - (x_*(x_*cos(kapa[camid])-y_*sin(kapa[camid]))/f+f*cos(kapa[camid]))*cos(omiga[camid]);
                    A(2*npixel+0, 1) = -f*sin(kapa[camid]) - x_*(x_*sin(kapa[camid])+y_*cos(kapa[camid]))/f;
                    A(2*npixel+0, 2) = y_;
                    A(2*npixel+0, 3) = (Rc(1,0)*f+Rc(1,2)*x_)/Xw_(2);
                    A(2*npixel+0, 4) = (Rc(2,0)*f+Rc(2,2)*x_)/Xw_(2);

                    A(2*npixel+1, 0) = -x_*sin(omiga[camid]) - (y_*(x_*cos(kapa[camid])-y_*sin(kapa[camid]))/f-f*sin(kapa[camid]))*cos(omiga[camid]);
                    A(2*npixel+1, 1) = -f*cos(kapa[camid]) - y_*(x_*sin(kapa[camid])+y_*cos(kapa[camid]))/f;
                    A(2*npixel+1, 2) = -x_;
                    A(2*npixel+1, 3) = (Rc(1,1)*f+Rc(1,2)*y_)/Xw_(2);
                    A(2*npixel+1, 4) = (Rc(2,1)*f+Rc(2,2)*y_)/Xw_(2);
                }
                
                A(2*npixel+0, 11+3*nP+0) = f*(Rc(0,2)*Xw_(0)-Rc(0,0)*Xw_(2))/Xw_(2)/Xw_(2);
                A(2*npixel+0, 11+3*nP+1) = f*(Rc(1,2)*Xw_(0)-Rc(1,0)*Xw_(2))/Xw_(2)/Xw_(2);
                A(2*npixel+0, 11+3*nP+2) = f*(Rc(2,2)*Xw_(0)-Rc(2,0)*Xw_(2))/Xw_(2)/Xw_(2);

                A(2*npixel+1, 11+3*nP+0) = f*(Rc(0,2)*Xw_(1)-Rc(0,1)*Xw_(2))/Xw_(2)/Xw_(2);
                A(2*npixel+1, 11+3*nP+1) = f*(Rc(1,2)*Xw_(1)-Rc(1,1)*Xw_(2))/Xw_(2)/Xw_(2);
                A(2*npixel+1, 11+3*nP+2) = f*(Rc(2,2)*Xw_(1)-Rc(2,1)*Xw_(2))/Xw_(2)/Xw_(2);

                L(2*npixel+0, 0) = x - x_;
                L(2*npixel+1, 0) = y - y_;

                npixel++;
            }
        }

        //cerr<<A<<endl;

        H = (A.transpose()*P*A).inverse();
        delta = H*(A.transpose()*P*L);

        phi[1] = phi[1] + delta(0);
        omiga[1] = omiga[1] + delta(1);
        kapa[1] = kapa[1] + delta(2);
        Ys[1] = Ys[1] + delta(3);
        Zs[1] = Zs[1] + delta(4);
        phi[2] = phi[2] + delta(5);
        omiga[2] = omiga[2] + delta(6);
        kapa[2] = kapa[2] + delta(7);
        Xs[2] = Xs[2] + delta(8);
        Ys[2] = Ys[2] + delta(9);
        Zs[2] = Zs[2] + delta(10);

        for(int i = 0;i<same_pt.size();i++){
            Ps[i](0) = Ps[i](0)+delta(11+i*3+0);
            Ps[i](1) = Ps[i](1)+delta(11+i*3+1);
            Ps[i](2) = Ps[i](2)+delta(11+i*3+2);
        }

        //选权
        v = A*delta-L;
        int r = nl*2-11-same_pt.size()*3;
        double sigma0 = sqrt((v.transpose()*P*v)(0,0)/r);
        //cerr<<sigma0<<endl;
        for(int i = 0;i<nl*2;i++){
            double pi = P(i,i);
            double sigmavi = sigma0*sqrt(1/pi-(A.row(i)*H*(A.row(i).transpose()))(0,0));
            //double sigmavi = sigma0*sqrt((A.row(i)*H*(A.row(i).transpose()))(0,0));
            double Vi = abs(v(i,0)/(sigmavi+0.0000001));
            //if(Vi>3){
            //    P(i,i) = pi*exp(-Vi/3.0);
            //}
            if(Vi>=4.13){
                P(i,i) = 1.0/(Vi*Vi);
            }else{
                P(i,i) = 1.0;
            }
        }

        //cerr<<delta.transpose().head(20)<<endl;
    }

    cerr<<endl;
    int npixel = 0;
    int nP = 0;
    for(unordered_map<int, vector<int>>::iterator it = same_pt.begin(); it != same_pt.end(); it++, nP++){
        int Pid = it->first;
        vector<int> sets = it->second;
        for(int ic = 0; ic<sets.size(); ic++){

            if(v(npixel*2+1)>4*0.0028){
            //if(P(npixel*2+1)<1.0){
                cerr<<"camid = "<<sets[ic]<<" id = "<<Pid<<" v = "<<v(npixel*2+1)<<endl;
            }

            npixel++;
        }
    }

    cerr<<"phi2, omiga2, kapa2, Xs2, Ys2, Zs2 :"<<phi[1]<<" "<<omiga[1]<<" "<<kapa[1]<<" "<<Xs[1]<<" "<<Ys[1]<<" "<<Zs[1]<<endl;
    cerr<<"phi3, omiga3, kapa3, Xs3, Ys3, Zs3 :"<<phi[2]<<" "<<omiga[2]<<" "<<kapa[2]<<" "<<Xs[2]<<" "<<Ys[2]<<" "<<Zs[2]<<endl;

}