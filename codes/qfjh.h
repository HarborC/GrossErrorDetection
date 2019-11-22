// 相对定向
#include "Common.h"

void qfjh(vector<int> idx, int c1_, int c2_, vector<unordered_map<int, Eigen::Vector2d>> images, int N){
    
    // Generate sets of 8 points for each RANSAC iteration
    int mMaxIterations = 1;
    vector<int> idx_ = idx;
    int cam1 = c1_;
    int cam2 = c2_;
    int nmatches = idx_.size();
    
    /*vector<size_t> vAllIndices;
    vector<size_t> vAvailableIndices;
    for(int i=0; i<nmatches; i++)
    {
        vAllIndices.push_back(i);
    }
    vector< vector<size_t> > mvSets = vector< vector<size_t> >(mMaxIterations,vector<size_t>(5,0));

    DUtils::Random::SeedRandOnce(0);

    for(int it=0; it<mMaxIterations; it++)
    {
        vAvailableIndices = vAllIndices;

        // Select a minimum set
        for(size_t j=0; j<5; j++)
        {
            int randi = DUtils::Random::RandomInt(0,vAvailableIndices.size()-1);
            int idx = vAvailableIndices[randi];

            mvSets[it][j] = idx;

            vAvailableIndices[randi] = vAvailableIndices.back();
            vAvailableIndices.pop_back();
        }
    }*/

    for(size_t it = 0; it<mMaxIterations; it++){
        //vector<size_t> set = mvSets[it];

        float phi = 0;
        float omiga = 0;
        float kapa = 0;
        float f = 100.5;
        float u = 0;
        float v = 0;
        float Bx = 0, By, Bz;

        //for(int j = 0; j<5; j++){
            //Bx += images[cam1][idx_[set[j]]](0)-images[cam2][idx_[set[j]]](0);
        for(int j = 0; j<idx_.size();j++){
            Bx += images[cam1][idx_[j]](0)-images[cam2][idx_[j]](0);
        }
        //Bx /= 5;
        Bx /= idx_.size();
        cerr<<endl;
        cerr<<Bx<<endl;
        cerr<<endl;
        
        Eigen::Matrix<double, Dynamic, 1> L;
        Eigen::Matrix<double, Dynamic, 5> A;
        Eigen::Matrix<double, Dynamic, Dynamic> P;
        Eigen::Matrix<double, Dynamic, 1> delta;
        Eigen::Matrix<double, Dynamic, Dynamic> H;
        Eigen::Matrix<double, Dynamic, 1> re;
        A.setZero(idx_.size(),5);
        L.setZero(idx_.size(),1);
        P.setZero(idx_.size(),idx_.size());
        delta.resize(idx_.size(),1);
        H.resize(5,5);
        re.resize(idx_.size(),1);
        for(int i = 0;i<idx_.size();i++){
            P(i,i) = 1.0;
        }

        for(int tim = 0; tim < N; tim++){
            By = Bx * u;
            Bz = Bx * v;

            Eigen::Matrix<double, 3, 3> R;
            R<<cos(phi)*cos(kapa)-sin(phi)*sin(omiga)*sin(kapa),-cos(phi)*sin(kapa)-sin(phi)*sin(omiga)*cos(kapa),-sin(phi)*cos(omiga),
                cos(omiga)*sin(kapa),cos(omiga)*cos(kapa),-sin(omiga),
                sin(phi)*cos(kapa)+cos(phi)*sin(omiga)*sin(kapa),-sin(phi)*sin(kapa)+cos(phi)*sin(omiga)*cos(kapa),cos(phi)*cos(omiga);

            for(size_t i = 0; i<idx_.size(); i++){
                Eigen::Vector3d X1(images[cam1][idx_[i]](0), images[cam1][idx_[i]](1), -f);
                Eigen::Vector3d X2_(images[cam2][idx_[i]](0), images[cam2][idx_[i]](1), -f);
                Eigen::Vector3d X2 = R * X2_;

                float N2 = (Bx * X1(2) - Bz * X1(0)) / (X1(0) * X2(2) - X1(2) * X2(0));
                float N1 = (Bx * X2(2) - Bz * X2(0)) / (X1(0) * X2(2) - X1(2) * X2(0));

                A(i,0) = Bx;
                A(i,1) = -Bx * X2(1) / X2(2);
                A(i,2) = - X2(0) * X2(1) * N2 / X2(2);
                A(i,3) = - (X2(2) + X2(1) * X2(1) / X2(2)) * N2;
                A(i,4) = X2(0) * N2;
                L(i,0) = N1*X1(1)-N2*X2(1)-By;
            }

            H = (A.transpose()*P*A).inverse();
            delta = H * (A.transpose() * P * L);
            u += delta(0);
            v += delta(1);
            phi += delta(2);
            omiga += delta(3);
            kapa += delta(4);
            cerr<<"u, v , phi, omiga, kapa :"<<delta.transpose()<<endl;

            //选权
            re = A*delta-L;
            int r = idx_.size()-5;
            double sigma0 = sqrt((re.transpose()*P*re)(0,0)/r);
            //cerr<<sigma0<<endl;
            for(int i = 0;i<idx_.size();i++){
                double pi = P(i,i);
                double sigmavi = sigma0*sqrt(1/pi-(A.row(i)*H*(A.row(i).transpose()))(0,0));
                //double sigmavi = sigma0*sqrt((A.row(i)*H*(A.row(i).transpose()))(0,0));
                double Vi = abs(re(i,0)/(sigmavi+0.0000001));
                if(Vi>3){
                    P(i,i) = pi*exp(1-(Vi/3.0)*(Vi/3.0));
                }
                //if(Vi>=4.13){
                //    P(i,i) = 1.0/(Vi*Vi);
                //}else{
                //    P(i,i) = 1.0;
                //}
            }

        }

        for(size_t i = 0; i<idx_.size(); i++){
            if(re(i)>4*0.0028){
                //if(P(npixel*2+1)<1.0){
                cerr<<"id = "<<idx_[i]<<" v = "<<re(i)<<endl;
            }

        }


        cerr<<"Bx, By, Bz, phi, omiga, kapa :"<<Bx<<" "<<By<<" "<<Bz<<" "<<phi<<" "<<omiga<<" "<<kapa<<" "<<endl;
    }
}
    