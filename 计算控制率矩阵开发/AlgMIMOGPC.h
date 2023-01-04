/*
多变量广义预测控制
*/
#pragma once
#include "AlgICSBase.h"
#include <malloc.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
#ifdef DCE_LINUX
#define MAX_PATH          260
#define max(a,b)            (((a) > (b)) ? (a) : (b))
typedef __gnu_cxx::vector<double> doublevector;
#else
	#if defined(_MSC_VER) && (_MSC_VER >= 1500)
	typedef std::vector<double> doublevector;
	#else
	typedef stdext::vector<double> doublevector;
	#endif
#endif

const int MIMOGPC_MODEL_NA_MAX=9;
const int MIMOGPC_MODEL_NB_MAX=9;
const int MIMOGPC_MODELNUM_MAX=16;
const int M1MIMOGPC_MODEL_NA_MAX=9;//三路外扰模型阶数及个数MR
const int M1MIMOGPC_MODEL_NB_MAX=9;
const int M1MIMOGPC_MODELNUM_MAX=16;
const int M2MIMOGPC_MODEL_NA_MAX=9;
const int M2MIMOGPC_MODEL_NB_MAX=9;
const int M2MIMOGPC_MODELNUM_MAX=16;
const int M3MIMOGPC_MODEL_NA_MAX=9;
const int M3MIMOGPC_MODEL_NB_MAX=9;
const int M3MIMOGPC_MODELNUM_MAX=16;


class AlgICSMIMOGPC :public AlgICSBase
{
public:
	AlgICSMIMOGPC(AlgIdType AlgId);
	virtual ~AlgICSMIMOGPC();

	virtual int Init(const char* pTunePara, AlgInOut *pInVar, AlgInOut *pOutVar, bool* pIsConnect);
	virtual int Reinit(const char* pTunePara, bool* pIsConnect);
	virtual int Run(unsigned int PeriodMs, AlgInOut * pInVar, AlgInOut *pOutVar);
	virtual int Trace(AlgInOut * pInVar, AlgInOut *pOutVar);
	virtual int SetTunePara(int ParaPos, const char* pParaVal);
	void SendSyncData(int* Len,char* pDataBuf);	
	void RecSyncData(int Len,char* pDataBuf);
	virtual void cleanup();
	virtual void AlgCmd(CmdIdType CmdID,float Data);

private:
	bool ValidatePara();//检修常数有效性
	bool ReadParay1p1FromFile();//读控制量1对过程量1模型
	bool ReadParay1p2FromFile();//读控制量1对过程量2模型
	bool ReadParay2p1FromFile();//读控制量2对过程量1模型
	bool ReadParay2p2FromFile();//读控制量2对过程量2模型
	bool ReadParav1p1FromFile();//读外扰1对过程量1模型
	bool ReadParav1p2FromFile();//读外扰1对过程量2模型
	bool ReadParav2p1FromFile();//读外扰2对过程量1模型
	bool ReadParav2p2FromFile();//读外扰2对过程量2模型
	bool ReadParav3p1FromFile();//读外扰3对过程量1模型
	bool ReadParav3p2FromFile();//读外扰3对过程量2模型
	bool malloc_nm(int n,int m,int size_1,int size_2,double*** a);//动态分配一个n行m列数组
	void free_nm(int n,double** a);//释放动态分配的一个n行m列数组
	void conv(double a[],double b[],int na,int nb,double c[]);//求卷积
	void multidiophantine(double a[],int na,double b[],int nb,int N,double** E,double** F,double** G);//多步丢番图方程求解
	//void multidiophantine2(double a[],int na,double b1[],int nb1,double b2[],int nb2,int N,double** E,double** F1,double** F2,double** G);
	void matrix_init(double** a,int n,int m,double a0);//将n行m列矩阵a每个值初始化
	//void gpcmatrix(double a[],int na,double b[],int nb,int d,int N,int Nu,double** G,double** G1,double** F);//计算GPC(广义预测控制)算法所需矩阵
	void gpcmatrix(double a[],int na,double b[],int nb,int d_OR,int d2,int N,int Nu,double** G,double** G1,double** F);//计算GPC(广义预测控制)算法所需矩阵
	void gpcmatrix_2(double a[],int na,double b1[],int nb1,int d1_OR,double b2[],int nb2,int d2_OR,int d_another_min,int N,int Nu,double** G,double** G1,double** F);
	bool GpcMatrixCalculate();//动态分配内存，根据广义预测控制原理计算相关矩阵的值
	void GpcMatrixClear();//释放动态分配的内存
	int Calculate_max(int a,int b);
	int Calculate_min(int a,int b);
	
	unsigned char GetInputQuality();//返回输入品质（最差），跟踪输入除外
	void myc2d(double* a,double* b,int na,int nb,double Ts,double* new_a,double* new_b);//连续系统离散化程序，采用双线性变换法（tustin），
	void myc2d2(double* a,double* b1,double* b2,int na,int nb1,int nb2,double Ts,double* new_a,double* new_b1,double* new_b2);//连续系统离散化程序，采用双线性变换法（tustin），
	
	void CRateCalculate(double **A ,double **B , vector<double>& result1 ,vector<double>& result2,int GMA ,int GMB ,int GN,float Gama1,float Gama2,float eta1,float eta2);//计算控制率
	void DeletePnt();
	void LUP_Descomposition(double *A,double *L,double *U,int *P,int N);//LUP分解
	void LUP_Solve(double *L,double *U,int *P,double *b,int N,double *x);//LUP求解方程
	int getNext(int i, int m, int n);/* 后继 */
	int getPre(int i, int m, int n);/* 前驱 */
	void movedata(double *mtx, int i, int m, int n);/* 处理以下标i为起点的环 */
	void transpose(double *mtx, int m, int n);/* 转置，即循环处理所有环 */
	void LUP_solve_inverse(double *A , int N , double *inv_A);//LUP求逆(将每列b求出的各列x进行组装)
	void MatrixCombine(double** G1,double** G2,int G1N,int G1M ,int G2N,int G2M ,double** GF);

private:
	//const int MAX_MODEL_NUM=4;
	//const int VAR_NUM_PER_MODEL=3;
	enum{IN_NUM=44,				// 输入变量个数 MR
		OUT_NUM=5,				// 输出变量个数 MR
		MIN_CONST_POS=0,		// 常数位置起始值（约定从1开始）
		MAX_CONST_POS=59,		// 常数位置最大值
		POS_I_SETLEFTTEM=0,		// 输入变量位置
		POS_I_REALLEFTTEM,
		POS_I_VALPOSOFSPRAY,
		POS_I_DELOFSPRAY1,		// 输入变量位置：第一组模型的第一个输入
		POS_I_MAXOFYVALUE=15,
		POS_I_MINOFYVALUE,
		POS_I_ALGSTIN=0,
		POS_O_OUT=0
	};

	// 按设计的常数类型、顺序声明一个包含全部常数的结构。如果4字节、2字节变量（double/short/int）混用，
	// 需要指定按2字节对齐
	struct sumTunePara
	{
		short Ts;//离散周期
		short N;//预测时域
		short Nu1;//控制时域1
		short Nu2;//控制时域2
		short TransPeriod;//切换过渡周期数
		short TransferModel;//切换erro校正模式和非erro校正模式
		float bak;
		float alpha1;//设定值1柔化系数
		float alpha2;//设定值2柔化系数
		float gama1;//控制权重系数1
		float gama2;//控制权重系数2
		float eta1;//过程量权重系数1
		float eta2;//过程量权重系数2
		
		//int GK;//控制律柔化系数
		
	}m_TunePara;
	
	//需要从文件中读取的参数
	struct ParaFromFile
	{	
		int totalmodel;//等于模型个数减1，从0计数
		int na;
		int nb;
		int nb2;
		double a[MIMOGPC_MODELNUM_MAX][MIMOGPC_MODEL_NA_MAX];
		double aa[MIMOGPC_MODELNUM_MAX][MIMOGPC_MODEL_NA_MAX+1];//aa=a*[1 -1];
		double b[MIMOGPC_MODELNUM_MAX][MIMOGPC_MODEL_NB_MAX];
		double b2[MIMOGPC_MODELNUM_MAX][MIMOGPC_MODEL_NB_MAX];
		int d[MIMOGPC_MODELNUM_MAX];
		int d2[MIMOGPC_MODELNUM_MAX];
		float T;
		int d_max;
		int d2_max;
		int flag[7];
		int NLTDflag;
		int NLTDna;
		int NLTDnb;
		double NLTD_a[MIMOGPC_MODELNUM_MAX][MIMOGPC_MODEL_NA_MAX];
		double NLTD_aa[MIMOGPC_MODELNUM_MAX][MIMOGPC_MODEL_NA_MAX+1];//aa=a*[1 -1];
		double NLTD_b[MIMOGPC_MODELNUM_MAX][MIMOGPC_MODEL_NB_MAX];
		int NLTD_d[MIMOGPC_MODELNUM_MAX];
		int NLTD_d_max;
	}m_y1p1_ModelPara,c_y1p1_ModelPara,m_y1p2_ModelPara,c_y1p2_ModelPara,
	m_y2p1_ModelPara,c_y2p1_ModelPara,m_y2p2_ModelPara,c_y2p2_ModelPara,
	m_v1p1_ModelPara,c_v1p1_ModelPara,m_v1p2_ModelPara,c_v1p2_ModelPara,
	m_v2p1_ModelPara,c_v2p1_ModelPara,m_v2p2_ModelPara,c_v2p2_ModelPara,
	m_v3p1_ModelPara,c_v3p1_ModelPara,m_v3p2_ModelPara,c_v3p2_ModelPara;//MR
	/*c_ModelPara为连续系统模型，m_ModelPara为离散系统模型
	  两路控制量定义为：m_y1p1_ModelPara,c_y1p1_ModelPara,m_y1p2_ModelPara,c_y1p2_ModelPara；m_y2p1_ModelPara,c_y2p1_ModelPara,m_y2p2_ModelPara,c_y2p2_ModelPara
	  三路外扰定义为：m_v1p1_ModelPara,c_v1p1_ModelPara,m_v1p2_ModelPara,c_v1p2_ModelPara；m_v2p1_ModelPara,c_v2p1_ModelPara,m_v2p2_ModelPara,c_v2p2_ModelPara；m_v3p1_ModelPara,c_v3p1_ModelPara,m_v3p2_ModelPara,c_v3p2_ModelPara
	*/
	
	//GPC算法所需矩阵，init()函数动态分配内存,cleanup()函数释放
	struct GpcMatrix
	{
		double** G;
		double** G1;
		double** F;
		//double* G2;
		double ym[10];//保存过去时刻预测输出值,数组长度大于被控对象参数a最大阶次;
		//doublevector GR;//控制率
		bool IsAlloc;//表示是否已经动态分配，参数错时不分配
	}m_y1p1_GpcMatrix[MIMOGPC_MODELNUM_MAX],m_y1p2_GpcMatrix[MIMOGPC_MODELNUM_MAX],
	m_y2p1_GpcMatrix[MIMOGPC_MODELNUM_MAX],m_y2p2_GpcMatrix[MIMOGPC_MODELNUM_MAX],
	m_v1p1_GpcMatrix[MIMOGPC_MODELNUM_MAX],m_v1p2_GpcMatrix[MIMOGPC_MODELNUM_MAX],
	m_v2p1_GpcMatrix[MIMOGPC_MODELNUM_MAX],m_v2p2_GpcMatrix[MIMOGPC_MODELNUM_MAX],
	m_v3p1_GpcMatrix[MIMOGPC_MODELNUM_MAX],m_v3p2_GpcMatrix[MIMOGPC_MODELNUM_MAX];
	

	struct GpcOCRMatrix
	{
		double** OCR;
		doublevector GR;//控制率
		bool IsAlloc;//表示是否已经动态分配，参数错时不分配
	}m_p1_OCR[MIMOGPC_MODELNUM_MAX],m_p2_OCR[MIMOGPC_MODELNUM_MAX];

	double yk1_lasttime[10];//保存过去时刻输出值,数组长度大于被控对象参数a最大阶次;
	//double yk1p2_lasttime[10];//保存过去时刻输出值,数组长度大于被控对象参数a最大阶次;
	double yk2_lasttime[10];//保存过去时刻输出值,数组长度大于被控对象参数a最大阶次;
	//double yk2p2_lasttime[10];//保存过去时刻输出值,数组长度大于被控对象参数a最大阶次;
	double du1p1;
	double du2p1;
	double du1p2;
	double du2p2;
	double du1;
	double du2;

	doublevector dUk1;//保存过去时刻控制器输出变化量,数组长度为b阶次加d; 
	doublevector dUk2;//保存过去时刻控制器输出变化量,数组长度为b阶次加d; 
	doublevector dv1k;//保存过去时刻外扰1输出变化量,数组长度为b阶次加d; MR
	doublevector dv2k;//保存过去时刻外扰2输出变化量,数组长度为b阶次加d; MR
	doublevector dv3k;//保存过去时刻外扰3输出变化量,数组长度为b阶次加d; MR
	//doublevector dunfk1;//保存过去时刻前馈输出变化量,数组长度为b阶次加d; MR
	//doublevector dunfk2;//保存过去时刻前馈输出变化量,数组长度为b阶次加d; MR
	doublevector dv1p1k;//临时修改dv1k的数组 长度是nb+d;
	doublevector dv1p2k;//临时修改dv1k的数组 长度是nb+d;
	doublevector dv2p1k;
	doublevector dv2p2k;
	doublevector dv3p1k;
	doublevector dv3p2k;
	doublevector dv1p1k_nltd;//临时修改dv1k的数组 长度是nb+d;
	doublevector dv1p2k_nltd;//临时修改dv1k的数组 长度是nb+d;
	doublevector dv2p1k_nltd;
	doublevector dv2p2k_nltd;
	doublevector dv3p1k_nltd;
	doublevector dv3p2k_nltd;
	//doublevector dnfk1;
	//doublevector dnfk2;
	doublevector Yr1;
	doublevector Yr2;
	//doublevector Y1MOUT;
	//int Y1MFLAG;
	//doublevector Y2MOUT;
	//doublevector Y3MOUT;
	double uk1;//保存上上时刻控制器1输出值
	double uk2;//保存上上时刻控制器2输出值
	double v1k;//保存上上时刻外扰1输出值 MR
	double v2k;//保存上上时刻外扰2输出值 MR
	double v3k;//保存上上时刻外扰3输出值 MR
	//double nufk1;//保存上上时刻前馈输出值 MR
	//double nufk2;//保存上上时刻前馈输出值 MR
	double u1_lasttime;//保存上一时刻控制器输出值
	double u2_lasttime;//保存上一时刻控制器输出值
	double uo_y1;//加主动前馈的控制器输出值
	double uo_y2;//加主动前馈的控制器输出值
	double v1;//保存上一时刻外扰1输出值 MR
	double v2;//保存上一时刻外扰2输出值 MR
	double v3;//保存上一时刻外扰2输出值 MR
	//double nuf1;//被动前馈量1 MR
	//double nuf2;//被动前馈量2 MR
	//double puf1;//主动前馈量1 MR
	//double puf2;//主动前馈量2 MR
	double w1;//衰减因子1 MR //注意读入顺序
	double w2;//衰减因子2 MR 
	double w3;//衰减因子3 MR
	double w4;//衰减因子4 MR
	double w5;//衰减因子5 MR
	//double w6;//衰减因子6 MR
	//double w7;//衰减因子7 MR
	double ny1v2;//NLTD微分输入1 MR
	double ny2v2;//NLTD微分输入2 MR
	double nv1v2;//NLTD微分输入3 MR
	double nv2v2;//NLTD微分输入4 MR
	double nv3v2;//NLTD微分输入5 MR
	//double nuf1v2;//NLTD微分输入6 MR
	//double nuf2v2;//NLTD微分输入7 MR
	double ScaleTop1;//控制量1输出上限 MR
	double ScaleBot1;//控制量1输出下限 MR
	double DeltaLH1;//增量1输出上限 MR
	double DeltaLL1;//增量1输出下限 MR
	double ScaleTop2;//控制量2输出上限 MR
	double ScaleBot2;//控制量2输出下限 MR
	double DeltaLH2;//增量2输出上限 MR
	double DeltaLL2;//增量2输出下限 MR
	bool m_IsParaError;
	int TimeCount;//控制区时间计数，达到预测控制采样时间，run一次==========改===============
	double pv1_now;//当前时刻过程量采样值
	double pv2_now;//当前时刻过程量采样值

	int Num_y1p1;//预测控制器采用的预测模型编号
	int sel_y1p1_now;//当前时刻选择的模型编号
	int sel_y1p1_last;//上一时刻选择的模型编号
	int sel_y1p1_transtime;//模型切换累计时间==========改===============
	int Num_y1p2;//预测控制器采用的预测模型编号
	int sel_y1p2_now;//当前时刻选择的模型编号
	int sel_y1p2_last;//上一时刻选择的模型编号
	int sel_y1p2_transtime;//模型切换累计时间==========改===============
	
	int Num_y2p1;//预测控制器采用的预测模型编号
	int sel_y2p1_now;//当前时刻选择的模型编号
	int sel_y2p1_last;//上一时刻选择的模型编号
	int sel_y2p1_transtime;//模型切换累计时间==========改===============
	int Num_y2p2;//预测控制器采用的预测模型编号
	int sel_y2p2_now;//当前时刻选择的模型编号
	int sel_y2p2_last;//上一时刻选择的模型编号
	int sel_y2p2_transtime;//模型切换累计时间==========改===============

	int Num_v1p1;//外扰采用的预测模型编号 ============以下是三路外扰 MR============
	int sel_v1p1_now;//外扰当前时刻选择的模型编号
	int sel_v1p1_last;//外扰上一时刻选择的模型编号
	int sel_v1p1_transtime;//外扰模型切换累计时间
	int Num_v1p2;//外扰采用的预测模型编号 ============以下是三路外扰 MR============
	int sel_v1p2_now;//外扰当前时刻选择的模型编号
	int sel_v1p2_last;//外扰上一时刻选择的模型编号
	int sel_v1p2_transtime;//外扰模型切换累计时间

	int Num_v2p1;//外扰采用的预测模型编号
	int sel_v2p1_now;//外扰当前时刻选择的模型编号
	int sel_v2p1_last;//外扰上一时刻选择的模型编号
	int sel_v2p1_transtime;//外扰模型切换累计时间
	int Num_v2p2;//外扰采用的预测模型编号
	int sel_v2p2_now;//外扰当前时刻选择的模型编号
	int sel_v2p2_last;//外扰上一时刻选择的模型编号
	int sel_v2p2_transtime;//外扰模型切换累计时间

	int Num_v3p1;//外扰采用的预测模型编号
	int sel_v3p1_now;//外扰当前时刻选择的模型编号
	int sel_v3p1_last;//外扰上一时刻选择的模型编号
	int sel_v3p1_transtime;//外扰模型切换累计时间
	int Num_v3p2;//外扰采用的预测模型编号
	int sel_v3p2_now;//外扰当前时刻选择的模型编号
	int sel_v3p2_last;//外扰上一时刻选择的模型编号
	int sel_v3p2_transtime;//外扰模型切换累计时间

	unsigned int AlgState;
	int DataLen;//发送数据长度
	int initflag;//初始化标志
    double erro_p1_now;
	double erro_p2_now;
	int transferflag;//模型重读标志位（置1重读）
	bool ModelUpdate;//模型重读引脚
	double GpcMatrix_y1mnow;//MGPC+当前控制量模型输出
	double GpcMatrix_y2mnow;//MGPC+当前控制量模型输出
	int MainDropFlag;

};
