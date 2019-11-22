#include "Common.h"

// Method
#include "dphfjh.h"
#include "qfjh.h"
#include "gsfpc.h"

using namespace std;
using namespace Eigen;

int main (int argc, char **argv)
{
    int N = atoi(argv[1]);

    // 同名像点、控制点、焦距
    vector<unordered_map<int, Eigen::Vector2d>> images(3);
    unordered_map<int, Eigen::Vector3d> control_pt;
    int nl = 0;

    // 读取同名像点
    ifstream in("../Data/SamePoints.DAT");
    if(!in){
        cout<<"文件打开失败！"<<endl;
        return 0;
    }

    string s;
    getline(in, s);
    for(int i = 0; i<3; i++){
        getline(in, s);
        
        while (1){
            getline(in, s);
            stringstream ss;
            ss<<s;
            int n1;
            double n2, n3;
            ss>>n1; 
            if(n1 <0){
                break;
            }
            ss>>n2; ss>>n3;
            Eigen::Vector2d pt(n2, n3);
            images[i][n1] = pt;
            nl++;
        }
    }

    in.close();

    /*for(int i = 0;i<3;i++){
        cout<<"Image No."<<i+1<<": "<<endl;
        for(unordered_map<int, Eigen::Vector2d>::iterator it = images[i].begin(); it != images[i].end(); it++){
            cout<<"id: "<<(*it).first<<" coord: x = "<<it->second(0)<<" y = "<<it->second(1)<<endl;
        }
        cout<<endl;
    }*/

    // 读取控制点
    ifstream in2("../Data/ControlPoint.DAT");
    if(!in2){
        cout<<"文件打开失败！"<<endl;
        return 0;
    }

    while(1){
        getline(in2, s);
        stringstream ss;
        ss<<s;
        int n1; double n2, n3, n4;
        ss>>n1; 
        if(n1 < 0)
            break;
        ss>>n2; ss>>n3; ss>>n4;
        Eigen::Vector3d c_pt(n2, n3, n4);
        control_pt[n1] = c_pt;
    }

    in2.close();

    /*for(unordered_map<int, Eigen::Vector3d>::iterator it = control_pt.begin(); it != control_pt.end(); it++){
        cout<<"id: "<<it->first<<" Control Point: X = "<<it->second(0)<<" Y = "<<it->second(1)<<" Z  = "<<it->second(2)<<endl;
    }
    cout<<endl;*/

    // 同名点与像片的关系
    vector<int> idx12, idx23, idx13;
    unordered_map<int, vector<int>> same_pt;
    
    for(unordered_map<int, Eigen::Vector2d>::iterator it = images[0].begin(); it != images[0].end(); it++){
        int index = it->first;
        if(images[1].find(index) != images[1].end()){
            idx12.push_back(index);
        }
        if(images[2].find(index) != images[2].end()){
            idx13.push_back(index);
        }
        if(same_pt.find(index) == same_pt.end()){
            vector<int> v;
            v.push_back(0);
            same_pt[index] = v;
        }
    }
    for(unordered_map<int, Eigen::Vector2d>::iterator it = images[1].begin(); it != images[1].end(); it++){
        int index = it->first;
        if(images[2].find(index) != images[2].end()){
            idx23.push_back(index);
        }
        if(same_pt.find(index) == same_pt.end()){
            vector<int> v;
            v.push_back(1);
            same_pt[index] = v;
        }else{
            same_pt[index].push_back(1);
        }
    }
    for(unordered_map<int, Eigen::Vector2d>::iterator it = images[2].begin(); it != images[2].end(); it++){
        int index = it->first;
        if(same_pt.find(index) == same_pt.end()){
            vector<int> v;
            v.push_back(2);
            same_pt[index] = v;
        }else{
            same_pt[index].push_back(2);
        }
    }
    
    int n1, n2;
    n1 = n2 = 0;
    for(unordered_map<int, vector<int>>::iterator it = same_pt.begin(); it != same_pt.end(); it++){
        if(it->second.size() == 2){
            n1++;
        }else if (it->second.size() == 3){
            n2++;
        }
    }

    /*for(unordered_map<int, vector<int>>::iterator it = same_pt.begin(); it != same_pt.end(); it++){
        cout<<"id = "<<it->first<<" images = ";
        for(int j = 0;j<it->second.size(); j++){
            cout<<it->second[j]<<" ";
        }
        cout<<endl;
    }*/

    // 单片后方交会
    dphfjh(images, control_pt, N);

    // 前方交会
    qfjh(idx12, 0, 1, images, N);

    // 光束法平差
    gsfpc(images, idx12, idx13, same_pt, N);

    return 0;
}