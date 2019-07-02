#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#define PI 3.1415926

using namespace std;

const int nodeNum = 402;//300;                 //网格节点个数
const int triNum = 800;//596;                  //三角形的个数
const int start = 0;

double nodeCoord[nodeNum][3];            //节点的坐标
int triInd[triNum][3];                   //三角形节点的编号

double nodeCoord1[nodeNum];
int triInd1[triNum];

double n[nodeNum][3];                    //节点的单位法向量
double triNor[triNum][3];                //每个三角形面上的单位法向量
double KH[nodeNum];                      //平均曲率

//子函数的声明
void getNodeCoord_TriInd();              //获取节点坐标和三角形的节点编号
void unitNorVec(int num, double a[3], double b[3], double c[3]);    //计算给定三角形的法向量
void calculateTriNor();                  //计算每个三角形的外单位法向量
double calculateomega(double a[3], double b[3], double c[3]);   //计算法向量中的权系数
void calculateN();                       //计算每个节点的外单位法向量
bool triType(double A[3], double B[3], double C[3]);     //判断三角形是锐角还是钝角
double acuteAngArea(double A[3], double B[3], double C[3]);      //计算锐角三角形的面积
double obtuseAngArea(double A[3], double B[3], double C[3]);      //计算钝角三角形的面积
double calculateAm(int num);            //计算节点周围小邻域的面积
double calcalateAngle(int num);                 //计算曲率后半部分表达式
void calculateKH();                    //计算得到所有曲率


py::array_t<double> curvature(py::array_t<double> input1, py::array_t<int> input2)
{
	py::array_t<double> output = py::array_t<double>(nodeNum);
	
	py::buffer_info buf1 = input1.request();
    py::buffer_info buf2 = input2.request();
	py::buffer_info buf3 = output.request();

    double* nodes = (double*)buf1.ptr;
    int*    trias =	   (int*)buf2.ptr;
	double* curva = (double*)buf3.ptr;

	for (int i = 0; i < nodeNum; i++){
		nodeCoord1[i] =i+2;
		nodeCoord[i][0] = nodes[i*3];
		nodeCoord[i][1] = nodes[i*3+1];
		nodeCoord[i][2] = nodes[i*3+2];
	}
	
	for (int i = 0; i < triNum; i++){
		triInd[i][0] = trias[i*3];
		triInd[i][1] = trias[i*3+1];
		triInd[i][2] = trias[i*3+2];
	}
	calculateTriNor();
	calculateN();
	calculateKH();
	for (int i = 0; i < nodeNum; i++){
		curva[i] = KH[i];
	}
	return output;
}
PYBIND11_MODULE(curvature,m){
    m.def("curvature", &curvature, "haven't reckon about the explanation yet");
}
// c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` curvature.cpp -o curvature`python3-config --extension-suffix`

void unitNorVec(int num, double a[3], double b[3], double c[3])    //计算给定三角形的法向量
{
	double ab[3], ac[3], normal[3];
	double mold;

	ab[0] = b[0] - a[0];
	ab[1] = b[1] - a[1];
	ab[2] = b[2] - a[2];

	ac[0] = c[0] - a[0];
	ac[1] = c[1] - a[1];
	ac[2] = c[2] - a[2];

	normal[0] = ab[1] * ac[2] - ab[2] * ac[1];
	normal[1] = ab[2] * ac[0] - ab[0] * ac[2];
	normal[2] = ab[0] * ac[1] - ab[1] * ac[0];

	mold = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

	triNor[num][0] = normal[0] / mold;
	triNor[num][1] = normal[1] / mold;
	triNor[num][2] = normal[2] / mold;
}

void calculateTriNor()                  //计算每个三角形的外单位法向量
{
	int i;
	for (i = 0; i < triNum; i++)
	{
		unitNorVec(i, nodeCoord[triInd[i][0] - start],
			nodeCoord[triInd[i][1] - start], nodeCoord[triInd[i][2] - start]);
	}
}

double calculateomega(double a[3], double b[3], double c[3])   //计算法向量中的权系数
{
	double coefficient, sum[3];                       //a,b,c中a[3]为所求顶点的的坐标

	sum[0] = (a[0] + b[0] + c[0]) / 3.0;
	sum[1] = (a[1] + b[1] + c[1]) / 3.0;
	sum[2] = (a[2] + b[2] + c[2]) / 3.0;

	sum[0] = sum[0] - a[0];
	sum[1] = sum[1] - a[1];
	sum[2] = sum[2] - a[2];

	coefficient = 1.0 / (sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]);

	return coefficient;
}

void calculateN()                     //计算每个节点的外单位法向量
{
	int i, node;
	double nMold, omega;
	for (node = 0; node < nodeNum; node++)
	{
		n[node][0] = 0;
		n[node][1] = 0;
		n[node][2] = 0;
		for (i = 0; i < triNum; i++)              //判断三角形的第一个节点编号
		{
			if ((node + start) == triInd[i][0])
			{
				omega = calculateomega(nodeCoord[triInd[i][0] - start],
					nodeCoord[triInd[i][1] - start], nodeCoord[triInd[i][2] - start]);
				n[node][0] += omega * triNor[i][0];
				n[node][1] += omega * triNor[i][1];
				n[node][2] += omega * triNor[i][2];
			}

			if ((node + start) == triInd[i][1])
			{
				omega = calculateomega(nodeCoord[triInd[i][1] - start],
					nodeCoord[triInd[i][2] - start], nodeCoord[triInd[i][0] - start]);
				n[node][0] += omega * triNor[i][0];
				n[node][1] += omega * triNor[i][1];
				n[node][2] += omega * triNor[i][2];
			}

			if ((node + start) == triInd[i][2])
			{
				omega = calculateomega(nodeCoord[triInd[i][2] - start],
					nodeCoord[triInd[i][0] - start], nodeCoord[triInd[i][1] - start]);
				n[node][0] += omega * triNor[i][0];
				n[node][1] += omega * triNor[i][1];
				n[node][2] += omega * triNor[i][2];
			}

		}
	}

	for (node = 0; node < nodeNum; node++)
	{
		nMold = sqrt(n[node][0] * n[node][0] + n[node][1] * n[node][1] +
			n[node][2] * n[node][2]);

		n[node][0] /= nMold;
		n[node][1] /= nMold;
		n[node][2] /= nMold;
	}
}


bool triType(double A[3], double B[3], double C[3])     //判断三角形是锐角还是钝角
{
	double AB[3], AC[3];
	double BA[3], BC[3];
	double ABAC, BABC;

	AB[0] = B[0] - A[0];
	AB[1] = B[1] - A[1];
	AB[2] = B[2] - A[2];

	AC[0] = C[0] - A[0];
	AC[1] = C[1] - A[1];
	AC[2] = C[2] - A[2];

	ABAC = (AB[0] * AC[0] + AB[1] * AC[1] + AB[2] * AC[2]) / (sqrt(AB[0] * AB[0] + AB[1] * AB[1]
		+ AB[2] * AB[2]) * sqrt(AC[0] * AC[0] + AC[1] * AC[1] + AC[2] * AC[2]));

	if ((180 * acos(ABAC) / PI) >= 90)             //钝角和直角返回2，锐角返回1
		return true;

	BA[0] = A[0] - B[0];
	BA[1] = A[1] - B[1];
	BA[2] = A[2] - B[2];

	BC[0] = C[0] - B[0];
	BC[1] = C[1] - B[1];
	BC[2] = C[2] - B[2];

	BABC = (BA[0] * BC[0] + BA[1] * BC[1] + BA[2] * BC[2]) / (sqrt(BA[0] * BA[0] + BA[1] * BA[1]
		+ BA[2] * BA[2]) * sqrt(BC[0] * BC[0] + BC[1] * BC[1] + BC[2] * BC[2]));

	if ((180 * acos(BABC) / PI) >= 90)
		return true;

	if (((180 * acos(ABAC) / PI) + (180 * acos(BABC) / PI)) < 90)
		return true;

	return false;
}

double acuteAngArea(double A[3], double B[3], double C[3])      //计算锐角三角形的面积
{
	int i;
	double tArea;
	double angle1, angle2, b, c;
	double ab, ac;
	double BA[3], BC[3];
	double CA[3], CB[3];

	for (i = 0; i < 3; i++)
	{
		BA[i] = A[i] - B[i];
		BC[i] = C[i] - B[i];
		CA[i] = A[i] - C[i];
		CB[i] = B[i] - C[i];
	}

	angle1 = (BA[0] * BC[0] + BA[1] * BC[1] + BA[2] * BC[2]) /
		(sqrt(BA[0] * BA[0] + BA[1] * BA[1] + BA[2] * BA[2]) *
			sqrt(BC[0] * BC[0] + BC[1] * BC[1] + BC[2] * BC[2]));
	b = angle1 / (sin(acos(angle1)));
	ab = sqrt(BA[0] * BA[0] + BA[1] * BA[1] + BA[2] * BA[2]);

	angle2 = (CA[0] * CB[0] + CA[1] * CB[1] + CA[2] * CB[2]) /
		(sqrt(CA[0] * CA[0] + CA[1] * CA[1] + CA[2] * CA[2]) *
			sqrt(CB[0] * CB[0] + CB[1] * CB[1] + CB[2] * CB[2]));
	c = angle2 / (sin(acos(angle2)));
	ac = sqrt(CA[0] * CA[0] + CA[1] * CA[1] + CA[2] * CA[2]);

	tArea = (ab * ab * c + ac * ac * b) / 8.0;

	return tArea;
}

double obtuseAngArea(double A[3], double B[3], double C[3])      //计算钝角三角形的面积
{
	int i;
	double D[3], E[3], F[3];
	double AD[3], AE[3], AF[3];
	double ADAF[3], AEAF[3];
	double tArea;

	for (i = 0; i < 3; i++)
	{
		D[i] = (A[i] + B[i]) / 2.0;
		E[i] = (A[i] + C[i]) / 2.0;
		F[i] = (B[i] + C[i]) / 2.0;
	}

	for (i = 0; i < 3; i++)
	{
		AD[i] = D[i] - A[i];
		AE[i] = E[i] - A[i];
		AF[i] = F[i] - A[i];
	}

	ADAF[0] = AD[1] * AF[2] - AF[1] * AD[2];
	ADAF[1] = AD[2] * AF[0] - AF[2] * AD[0];
	ADAF[2] = AD[0] * AF[1] - AF[0] * AD[1];
	tArea = sqrt(ADAF[0] * ADAF[0] + ADAF[1] * ADAF[1] + ADAF[2] * ADAF[2]) / 2.0;

	AEAF[0] = AE[1] * AF[2] - AF[1] * AE[2];
	AEAF[1] = AE[2] * AF[0] - AF[2] * AE[0];
	AEAF[2] = AE[0] * AF[1] - AF[0] * AE[1];
	tArea += sqrt(AEAF[0] * AEAF[0] + AEAF[1] * AEAF[1] + AEAF[2] * AEAF[2]) / 2.0;

	return tArea;
}

double calculateAm(int num)            //计算节点周围小邻域的面积
{
	int i;
	bool type;
	double Am;

	Am = 0.0;
	for (i = 0; i < triNum; i++)
	{
		if ((num + start) == triInd[i][0])            //判断三角形的第一个节点编号
		{
			type = triType(nodeCoord[triInd[i][0] - start], nodeCoord[triInd[i][1] - start],
				nodeCoord[triInd[i][2] - start]);
			if (type)                                //判断是锐角还是钝角
			{
				Am += obtuseAngArea(nodeCoord[triInd[i][0] - start], nodeCoord[triInd[i][1] - start],
					nodeCoord[triInd[i][2] - start]);
			}
			else
			{
				Am += acuteAngArea(nodeCoord[triInd[i][0] - start], nodeCoord[triInd[i][1] - start],
					nodeCoord[triInd[i][2] - start]);
			}
		}

		if ((num + start) == triInd[i][1])            //判断三角形的第二个节点编号
		{
			type = triType(nodeCoord[triInd[i][1] - start], nodeCoord[triInd[i][2] - start],
				nodeCoord[triInd[i][0] - start]);
			if (type)                           //判断是锐角还是钝角
			{
				Am += obtuseAngArea(nodeCoord[triInd[i][1] - start], nodeCoord[triInd[i][2] - start],
					nodeCoord[triInd[i][0] - start]);
			}
			else
			{
				Am += acuteAngArea(nodeCoord[triInd[i][1] - start], nodeCoord[triInd[i][2] - start],
					nodeCoord[triInd[i][0] - start]);
			}
		}

		if ((num + start) == triInd[i][2])            //判断三角形的第三个节点编号
		{
			type = triType(nodeCoord[triInd[i][2] - start], nodeCoord[triInd[i][0] - start],
				nodeCoord[triInd[i][1] - start]);
			if (type)                  //判断是锐角还是钝角
			{
				Am += obtuseAngArea(nodeCoord[triInd[i][2] - start], nodeCoord[triInd[i][0] - start],
					nodeCoord[triInd[i][1] - start]);
			}
			else
			{
				Am += acuteAngArea(nodeCoord[triInd[i][2] - start], nodeCoord[triInd[i][0] - start],
					nodeCoord[triInd[i][1] - start]);
			}
		}
	}

	return Am;
}

double calcalateAngle(int num)                 //计算曲率后半部分表达式
{
	int i, j;
	int index[20][2];                     //存放包含num节点的三角形的节点编号
	int sum = 0;                          //统计num节点周围有多少个三角形
	int p1[20];                           //节点编号
	bool p2[20];
	double A[20][3];                      //向量Pj - Pi
	double B1[20][3];
	double B2[20][3];
	double C1[20][3];
	double C2[20][3];
	double arfa[20][2];
	double beta[20][2];
	double latter;

	for (i = 0; i < triNum; i++)        //判断节点num周围有哪些节点
	{
		if ((num + start) == triInd[i][0])
		{
			index[sum][0] = triInd[i][1];
			index[sum][1] = triInd[i][2];
			sum++;
		}

		if ((num + start) == triInd[i][1])
		{
			index[sum][0] = triInd[i][2];
			index[sum][1] = triInd[i][0];
			sum++;
		}

		if ((num + start) == triInd[i][2])
		{
			index[sum][0] = triInd[i][0];
			index[sum][1] = triInd[i][1];
			sum++;
		}
	}

	for (i = 0; i < sum; i++)
		p2[i] = true;
	p1[0] = index[0][0];
	p1[1] = index[0][1];

	for (i = 1; i < sum - 1; i++)
	{
		for (j = 1; j < sum; j++)
		{
			if (p1[i] == index[j][0] && p2[j])
			{
				p1[i + 1] = index[j][1];
				p2[j] = false;
				break;
			}

			if (p1[i] == index[j][1] && p2[j])
			{
				p1[i + 1] = index[j][0];
				p2[j] = false;
				break;
			}
		}
	}
	for (i = 0; i < 3; i++)               //计算第0个节点
	{
		A[0][i] = nodeCoord[num][i] - nodeCoord[p1[0] - start][i];
		B1[0][i] = nodeCoord[num][i] - nodeCoord[p1[sum - 1] - start][i];
		B2[0][i] = nodeCoord[p1[0] - start][i] - nodeCoord[p1[sum - 1] - start][i];
		C1[0][i] = nodeCoord[num][i] - nodeCoord[p1[1] - start][i];
		C2[0][i] = nodeCoord[p1[0] - start][i] - nodeCoord[p1[1] - start][i];
	}
	arfa[0][0] = (B1[0][0] * B2[0][0] + B1[0][1] * B2[0][1] + B1[0][2] * B2[0][2]) /
		(sqrt(B1[0][0] * B1[0][0] + B1[0][1] * B1[0][1] + B1[0][2] * B1[0][2]) *
			sqrt(B2[0][0] * B2[0][0] + B2[0][1] * B2[0][1] + B2[0][2] * B2[0][2]));
	arfa[0][1] = arfa[0][0] / (sin(acos(arfa[0][0])));

	beta[0][0] = (C1[0][0] * C2[0][0] + C1[0][1] * C2[0][1] + C1[0][2] * C2[0][2]) /
		(sqrt(C1[0][0] * C1[0][0] + C1[0][1] * C1[0][1] + C1[0][2] * C1[0][2]) *
			sqrt(C2[0][0] * C2[0][0] + C2[0][1] * C2[0][1] + C2[0][2] * C2[0][2]));
	beta[0][1] = beta[0][0] / (sin(acos(beta[0][0])));

	for (i = 1; i < sum - 1; i++)       //计算第一个到第sum - 2个节点
	{
		for (j = 0; j < 3; j++)
		{
			A[i][j] = nodeCoord[num][j] - nodeCoord[p1[i] - start][j];
			B1[i][j] = nodeCoord[num][j] - nodeCoord[p1[i - 1] - start][j];
			B2[i][j] = nodeCoord[p1[i] - start][j] - nodeCoord[p1[i - 1] - start][j];
			C1[i][j] = nodeCoord[num][j] - nodeCoord[p1[i + 1] - start][j];
			C2[i][j] = nodeCoord[p1[i] - start][j] - nodeCoord[p1[i + 1] - start][j];
		}
		arfa[i][0] = (B1[i][0] * B2[i][0] + B1[i][1] * B2[i][1] + B1[i][2] * B2[i][2]) /
			(sqrt(B1[i][0] * B1[i][0] + B1[i][1] * B1[i][1] + B1[i][2] * B1[i][2]) *
				sqrt(B2[i][0] * B2[i][0] + B2[i][1] * B2[i][1] + B2[i][2] * B2[i][2]));
		arfa[i][1] = arfa[i][0] / (sin(acos(arfa[i][0])));

		beta[i][0] = (C1[i][0] * C2[i][0] + C1[i][1] * C2[i][1] + C1[i][2] * C2[i][2]) /
			(sqrt(C1[i][0] * C1[i][0] + C1[i][1] * C1[i][1] + C1[i][2] * C1[i][2]) *
				sqrt(C2[i][0] * C2[i][0] + C2[i][1] * C2[i][1] + C2[i][2] * C2[i][2]));
		beta[i][1] = beta[i][0] / (sin(acos(beta[i][0])));
	}

	for (i = 0; i < 3; i++)               //计算第sum - 1个节点
	{
		A[sum - 1][i] = nodeCoord[num][i] - nodeCoord[p1[sum - 1] - start][i];
		B1[sum - 1][i] = nodeCoord[num][i] - nodeCoord[p1[sum - 2] - start][i];
		B2[sum - 1][i] = nodeCoord[p1[sum - 1] - start][i] - nodeCoord[p1[sum - 2] - start][i];
		C1[sum - 1][i] = nodeCoord[num][i] - nodeCoord[p1[0] - start][i];
		C2[sum - 1][i] = nodeCoord[p1[sum - 1] - start][i] - nodeCoord[p1[0] - start][i];
	}
	arfa[sum - 1][0] = (B1[sum - 1][0] * B2[sum - 1][0] + B1[sum - 1][1] * B2[sum - 1][1] + B1[sum - 1][2] * B2[sum - 1][2]) /
		(sqrt(B1[sum - 1][0] * B1[sum - 1][0] + B1[sum - 1][1] * B1[sum - 1][1] + B1[sum - 1][2] * B1[sum - 1][2]) *
			sqrt(B2[sum - 1][0] * B2[sum - 1][0] + B2[sum - 1][1] * B2[sum - 1][1] + B2[sum - 1][2] * B2[sum - 1][2]));
	arfa[sum - 1][1] = arfa[sum - 1][0] / (sin(acos(arfa[sum - 1][0])));

	beta[sum - 1][0] = (C1[sum - 1][0] * C2[sum - 1][0] + C1[sum - 1][1] * C2[sum - 1][1] + C1[sum - 1][2] * C2[sum - 1][2]) /
		(sqrt(C1[sum - 1][0] * C1[sum - 1][0] + C1[sum - 1][1] * C1[sum - 1][1] + C1[sum - 1][2] * C1[sum - 1][2]) *
			sqrt(C2[sum - 1][0] * C2[sum - 1][0] + C2[sum - 1][1] * C2[sum - 1][1] + C2[sum - 1][2] * C2[sum - 1][2]));
	beta[sum - 1][1] = beta[sum - 1][0] / (sin(acos(beta[sum - 1][0])));

	latter = 0.0;
	for (i = 0; i < sum; i++)
	{
		latter += (arfa[i][1] + beta[i][1]) * (A[i][0] * n[num][0] + A[i][1] * n[num][1]
			+ A[i][2] * n[num][2]);
	}

	return latter;
}

void calculateKH()           //计算得到所有曲率
{
	int i;
	for (i = 0; i < nodeNum; i++)
	{
		KH[i] = (0.25 / calculateAm(i)) * calcalateAngle(i);
	}
}

