#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

double dist(const double *p1, const double *p2){
    return ((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
}
void triangle_normal(const double *a, const double *b, const double *c, double *normal)
{
    
    double ab[3], ac[3];

    ab[0] = b[0] - a[0];
    ab[1] = b[1] - a[1];
    ab[2] = b[2] - a[2];

    ac[0] = c[0] - a[0];
    ac[1] = c[1] - a[1];
    ac[2] = c[2] - a[2];

    normal[0] = ab[1] * ac[2] - ab[2] * ac[1];
    normal[1] = ab[2] * ac[0] - ab[0] * ac[2];
    normal[2] = ab[0] * ac[1] - ab[1] * ac[0];

    double norm = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

    normal[0] /= norm;
    normal[1] /= norm;
    normal[2] /= norm;
}

vector<double> *out_normal(vector<double> &points,  vector<size_t> &triangles){
    size_t num_points = points.size()/3;
    size_t num_triangles = triangles.size()/3;
    cout<<points.size()<<endl;
    cout<<triangles.size()<<endl;
    vector<double> *out_normal = new vector<double>(num_points*3);
    for (size_t i=0; i<num_triangles; i++){
        auto p1_index = triangles[i*3+0];
        auto p2_index = triangles[i*3+1];
        auto p3_index = triangles[i*3+2];
        double *p1 = &(points[p1_index*3]);
        double *p2 = &(points[p2_index*3]);
        double *p3 = &(points[p3_index*3]);
        double local_normal[3];
        triangle_normal(p1,p2,p3,local_normal);
        // calculate weight.
        double pc[3];
        pc[0] = (p1[0]+p2[0]+p3[0])/3.0;
        pc[1] = (p1[1]+p2[1]+p3[1])/3.0;
        pc[2] = (p1[2]+p2[2]+p3[2])/3.0;
        //cout <<"center point:"<< pc[0] << pc[1] << pc[2]<< endl;
        //cout <<"out normal:"<< local_normal[0] << local_normal[1] << local_normal[2]<< endl;
        double d1 = dist(p1,pc);
        double d2 = dist(p2,pc);
        double d3 = dist(p3,pc);
        (*out_normal)[p1_index*3+0] += local_normal[0]/d1;
        (*out_normal)[p1_index*3+1] += local_normal[1]/d1;
        (*out_normal)[p1_index*3+2] += local_normal[2]/d1;
        (*out_normal)[p2_index*3+0] += local_normal[0]/d2;
        (*out_normal)[p2_index*3+1] += local_normal[1]/d2;
        (*out_normal)[p2_index*3+2] += local_normal[2]/d2;
        (*out_normal)[p3_index*3+0] += local_normal[0]/d3;
        (*out_normal)[p3_index*3+1] += local_normal[1]/d3;
        (*out_normal)[p3_index*3+2] += local_normal[2]/d3;
    }
            
    for (size_t i=0; i<num_points; i++){
        double *normalize = &(*out_normal)[i*3];
        double norm = sqrt(normalize[0] * normalize[0] + normalize[1] * normalize[1] + normalize[2] * normalize[2]);
        (*out_normal)[i*3+0] /= norm;
        (*out_normal)[i*3+1] /= norm;
        (*out_normal)[i*3+2] /= norm;
        double sum=0.0;
        sum += (points[i * 3 + 0] * (*out_normal)[i * 3 + 0]) ;
        sum += (points[i * 3 + 1] * (*out_normal)[i * 3 + 1]) ;
        sum += (points[i * 3 + 2] * (*out_normal)[i * 3 + 2]) ;
        double length = 0.0;
        length += points[i*3+0]*points[3*i+0]; 
        length += points[i*3+1]*points[3*i+1]; 
        length += points[i*3+2]*points[3*i+2]; 
        cout<<"error norm:"<<length + sum<<endl;
    }
    return out_normal;
}

void getNodeCoord_TriInd( vector<double> &points,  vector<size_t> &triangles)
{
    double double_temp;
    size_t size_t_temp;
	string node="points0.txt";
    string triangle="triangles0.txt";
	ifstream file;
	file.open(node,  ios::in);
    while(! file.eof())
    {
        file >> double_temp;
        points.push_back(double_temp);
    }
    points.pop_back();
	file.close();
	file.open(triangle,  ios::in);
    while(! file.eof())
    {
        file >> size_t_temp;
        triangles.push_back(size_t_temp);
    }
    triangles.pop_back();
	file.close();
}
int main(){
    vector<double> points;
    vector<size_t> triangles;

    getNodeCoord_TriInd(points, triangles);
    auto a = out_normal(points, triangles);

    cout << points.size()<< endl;
    cout << triangles.size()<< endl;
    return 0;
}

