/*
���������Ԥ�����
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
const int M1MIMOGPC_MODEL_NA_MAX=9;//��·����ģ�ͽ���������MR
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
	bool ValidatePara();//���޳�����Ч��
	bool ReadParay1p1FromFile();//��������1�Թ�����1ģ��
	bool ReadParay1p2FromFile();//��������1�Թ�����2ģ��
	bool ReadParay2p1FromFile();//��������2�Թ�����1ģ��
	bool ReadParay2p2FromFile();//��������2�Թ�����2ģ��
	bool ReadParav1p1FromFile();//������1�Թ�����1ģ��
	bool ReadParav1p2FromFile();//������1�Թ�����2ģ��
	bool ReadParav2p1FromFile();//������2�Թ�����1ģ��
	bool ReadParav2p2FromFile();//������2�Թ�����2ģ��
	bool ReadParav3p1FromFile();//������3�Թ�����1ģ��
	bool ReadParav3p2FromFile();//������3�Թ�����2ģ��
	bool malloc_nm(int n,int m,int size_1,int size_2,double*** a);//��̬����һ��n��m������
	void free_nm(int n,double** a);//�ͷŶ�̬�����һ��n��m������
	void conv(double a[],double b[],int na,int nb,double c[]);//����
	void multidiophantine(double a[],int na,double b[],int nb,int N,double** E,double** F,double** G);//�ಽ����ͼ�������
	//void multidiophantine2(double a[],int na,double b1[],int nb1,double b2[],int nb2,int N,double** E,double** F1,double** F2,double** G);
	void matrix_init(double** a,int n,int m,double a0);//��n��m�о���aÿ��ֵ��ʼ��
	//void gpcmatrix(double a[],int na,double b[],int nb,int d,int N,int Nu,double** G,double** G1,double** F);//����GPC(����Ԥ�����)�㷨�������
	void gpcmatrix(double a[],int na,double b[],int nb,int d_OR,int d2,int N,int Nu,double** G,double** G1,double** F);//����GPC(����Ԥ�����)�㷨�������
	void gpcmatrix_2(double a[],int na,double b1[],int nb1,int d1_OR,double b2[],int nb2,int d2_OR,int d_another_min,int N,int Nu,double** G,double** G1,double** F);
	bool GpcMatrixCalculate();//��̬�����ڴ棬���ݹ���Ԥ�����ԭ�������ؾ����ֵ
	void GpcMatrixClear();//�ͷŶ�̬������ڴ�
	int Calculate_max(int a,int b);
	int Calculate_min(int a,int b);
	
	unsigned char GetInputQuality();//��������Ʒ�ʣ����������������
	void myc2d(double* a,double* b,int na,int nb,double Ts,double* new_a,double* new_b);//����ϵͳ��ɢ�����򣬲���˫���Ա任����tustin����
	void myc2d2(double* a,double* b1,double* b2,int na,int nb1,int nb2,double Ts,double* new_a,double* new_b1,double* new_b2);//����ϵͳ��ɢ�����򣬲���˫���Ա任����tustin����
	
	void CRateCalculate(double **A ,double **B , vector<double>& result1 ,vector<double>& result2,int GMA ,int GMB ,int GN,float Gama1,float Gama2,float eta1,float eta2);//���������
	void DeletePnt();
	void LUP_Descomposition(double *A,double *L,double *U,int *P,int N);//LUP�ֽ�
	void LUP_Solve(double *L,double *U,int *P,double *b,int N,double *x);//LUP��ⷽ��
	int getNext(int i, int m, int n);/* ��� */
	int getPre(int i, int m, int n);/* ǰ�� */
	void movedata(double *mtx, int i, int m, int n);/* �������±�iΪ���Ļ� */
	void transpose(double *mtx, int m, int n);/* ת�ã���ѭ���������л� */
	void LUP_solve_inverse(double *A , int N , double *inv_A);//LUP����(��ÿ��b����ĸ���x������װ)
	void MatrixCombine(double** G1,double** G2,int G1N,int G1M ,int G2N,int G2M ,double** GF);

private:
	//const int MAX_MODEL_NUM=4;
	//const int VAR_NUM_PER_MODEL=3;
	enum{IN_NUM=44,				// ����������� MR
		OUT_NUM=5,				// ����������� MR
		MIN_CONST_POS=0,		// ����λ����ʼֵ��Լ����1��ʼ��
		MAX_CONST_POS=59,		// ����λ�����ֵ
		POS_I_SETLEFTTEM=0,		// �������λ��
		POS_I_REALLEFTTEM,
		POS_I_VALPOSOFSPRAY,
		POS_I_DELOFSPRAY1,		// �������λ�ã���һ��ģ�͵ĵ�һ������
		POS_I_MAXOFYVALUE=15,
		POS_I_MINOFYVALUE,
		POS_I_ALGSTIN=0,
		POS_O_OUT=0
	};

	// ����Ƶĳ������͡�˳������һ������ȫ�������Ľṹ�����4�ֽڡ�2�ֽڱ�����double/short/int�����ã�
	// ��Ҫָ����2�ֽڶ���
	struct sumTunePara
	{
		short Ts;//��ɢ����
		short N;//Ԥ��ʱ��
		short Nu1;//����ʱ��1
		short Nu2;//����ʱ��2
		short TransPeriod;//�л�����������
		short TransferModel;//�л�erroУ��ģʽ�ͷ�erroУ��ģʽ
		float bak;
		float alpha1;//�趨ֵ1�ữϵ��
		float alpha2;//�趨ֵ2�ữϵ��
		float gama1;//����Ȩ��ϵ��1
		float gama2;//����Ȩ��ϵ��2
		float eta1;//������Ȩ��ϵ��1
		float eta2;//������Ȩ��ϵ��2
		
		//int GK;//�������ữϵ��
		
	}m_TunePara;
	
	//��Ҫ���ļ��ж�ȡ�Ĳ���
	struct ParaFromFile
	{	
		int totalmodel;//����ģ�͸�����1����0����
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
	/*c_ModelParaΪ����ϵͳģ�ͣ�m_ModelParaΪ��ɢϵͳģ��
	  ��·����������Ϊ��m_y1p1_ModelPara,c_y1p1_ModelPara,m_y1p2_ModelPara,c_y1p2_ModelPara��m_y2p1_ModelPara,c_y2p1_ModelPara,m_y2p2_ModelPara,c_y2p2_ModelPara
	  ��·���Ŷ���Ϊ��m_v1p1_ModelPara,c_v1p1_ModelPara,m_v1p2_ModelPara,c_v1p2_ModelPara��m_v2p1_ModelPara,c_v2p1_ModelPara,m_v2p2_ModelPara,c_v2p2_ModelPara��m_v3p1_ModelPara,c_v3p1_ModelPara,m_v3p2_ModelPara,c_v3p2_ModelPara
	*/
	
	//GPC�㷨�������init()������̬�����ڴ�,cleanup()�����ͷ�
	struct GpcMatrix
	{
		double** G;
		double** G1;
		double** F;
		//double* G2;
		double ym[10];//�����ȥʱ��Ԥ�����ֵ,���鳤�ȴ��ڱ��ض������a���״�;
		//doublevector GR;//������
		bool IsAlloc;//��ʾ�Ƿ��Ѿ���̬���䣬������ʱ������
	}m_y1p1_GpcMatrix[MIMOGPC_MODELNUM_MAX],m_y1p2_GpcMatrix[MIMOGPC_MODELNUM_MAX],
	m_y2p1_GpcMatrix[MIMOGPC_MODELNUM_MAX],m_y2p2_GpcMatrix[MIMOGPC_MODELNUM_MAX],
	m_v1p1_GpcMatrix[MIMOGPC_MODELNUM_MAX],m_v1p2_GpcMatrix[MIMOGPC_MODELNUM_MAX],
	m_v2p1_GpcMatrix[MIMOGPC_MODELNUM_MAX],m_v2p2_GpcMatrix[MIMOGPC_MODELNUM_MAX],
	m_v3p1_GpcMatrix[MIMOGPC_MODELNUM_MAX],m_v3p2_GpcMatrix[MIMOGPC_MODELNUM_MAX];
	

	struct GpcOCRMatrix
	{
		double** OCR;
		doublevector GR;//������
		bool IsAlloc;//��ʾ�Ƿ��Ѿ���̬���䣬������ʱ������
	}m_p1_OCR[MIMOGPC_MODELNUM_MAX],m_p2_OCR[MIMOGPC_MODELNUM_MAX];

	double yk1_lasttime[10];//�����ȥʱ�����ֵ,���鳤�ȴ��ڱ��ض������a���״�;
	//double yk1p2_lasttime[10];//�����ȥʱ�����ֵ,���鳤�ȴ��ڱ��ض������a���״�;
	double yk2_lasttime[10];//�����ȥʱ�����ֵ,���鳤�ȴ��ڱ��ض������a���״�;
	//double yk2p2_lasttime[10];//�����ȥʱ�����ֵ,���鳤�ȴ��ڱ��ض������a���״�;
	double du1p1;
	double du2p1;
	double du1p2;
	double du2p2;
	double du1;
	double du2;

	doublevector dUk1;//�����ȥʱ�̿���������仯��,���鳤��Ϊb�״μ�d; 
	doublevector dUk2;//�����ȥʱ�̿���������仯��,���鳤��Ϊb�״μ�d; 
	doublevector dv1k;//�����ȥʱ������1����仯��,���鳤��Ϊb�״μ�d; MR
	doublevector dv2k;//�����ȥʱ������2����仯��,���鳤��Ϊb�״μ�d; MR
	doublevector dv3k;//�����ȥʱ������3����仯��,���鳤��Ϊb�״μ�d; MR
	//doublevector dunfk1;//�����ȥʱ��ǰ������仯��,���鳤��Ϊb�״μ�d; MR
	//doublevector dunfk2;//�����ȥʱ��ǰ������仯��,���鳤��Ϊb�״μ�d; MR
	doublevector dv1p1k;//��ʱ�޸�dv1k������ ������nb+d;
	doublevector dv1p2k;//��ʱ�޸�dv1k������ ������nb+d;
	doublevector dv2p1k;
	doublevector dv2p2k;
	doublevector dv3p1k;
	doublevector dv3p2k;
	doublevector dv1p1k_nltd;//��ʱ�޸�dv1k������ ������nb+d;
	doublevector dv1p2k_nltd;//��ʱ�޸�dv1k������ ������nb+d;
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
	double uk1;//��������ʱ�̿�����1���ֵ
	double uk2;//��������ʱ�̿�����2���ֵ
	double v1k;//��������ʱ������1���ֵ MR
	double v2k;//��������ʱ������2���ֵ MR
	double v3k;//��������ʱ������3���ֵ MR
	//double nufk1;//��������ʱ��ǰ�����ֵ MR
	//double nufk2;//��������ʱ��ǰ�����ֵ MR
	double u1_lasttime;//������һʱ�̿��������ֵ
	double u2_lasttime;//������һʱ�̿��������ֵ
	double uo_y1;//������ǰ���Ŀ��������ֵ
	double uo_y2;//������ǰ���Ŀ��������ֵ
	double v1;//������һʱ������1���ֵ MR
	double v2;//������һʱ������2���ֵ MR
	double v3;//������һʱ������2���ֵ MR
	//double nuf1;//����ǰ����1 MR
	//double nuf2;//����ǰ����2 MR
	//double puf1;//����ǰ����1 MR
	//double puf2;//����ǰ����2 MR
	double w1;//˥������1 MR //ע�����˳��
	double w2;//˥������2 MR 
	double w3;//˥������3 MR
	double w4;//˥������4 MR
	double w5;//˥������5 MR
	//double w6;//˥������6 MR
	//double w7;//˥������7 MR
	double ny1v2;//NLTD΢������1 MR
	double ny2v2;//NLTD΢������2 MR
	double nv1v2;//NLTD΢������3 MR
	double nv2v2;//NLTD΢������4 MR
	double nv3v2;//NLTD΢������5 MR
	//double nuf1v2;//NLTD΢������6 MR
	//double nuf2v2;//NLTD΢������7 MR
	double ScaleTop1;//������1������� MR
	double ScaleBot1;//������1������� MR
	double DeltaLH1;//����1������� MR
	double DeltaLL1;//����1������� MR
	double ScaleTop2;//������2������� MR
	double ScaleBot2;//������2������� MR
	double DeltaLH2;//����2������� MR
	double DeltaLL2;//����2������� MR
	bool m_IsParaError;
	int TimeCount;//������ʱ��������ﵽԤ����Ʋ���ʱ�䣬runһ��==========��===============
	double pv1_now;//��ǰʱ�̹���������ֵ
	double pv2_now;//��ǰʱ�̹���������ֵ

	int Num_y1p1;//Ԥ����������õ�Ԥ��ģ�ͱ��
	int sel_y1p1_now;//��ǰʱ��ѡ���ģ�ͱ��
	int sel_y1p1_last;//��һʱ��ѡ���ģ�ͱ��
	int sel_y1p1_transtime;//ģ���л��ۼ�ʱ��==========��===============
	int Num_y1p2;//Ԥ����������õ�Ԥ��ģ�ͱ��
	int sel_y1p2_now;//��ǰʱ��ѡ���ģ�ͱ��
	int sel_y1p2_last;//��һʱ��ѡ���ģ�ͱ��
	int sel_y1p2_transtime;//ģ���л��ۼ�ʱ��==========��===============
	
	int Num_y2p1;//Ԥ����������õ�Ԥ��ģ�ͱ��
	int sel_y2p1_now;//��ǰʱ��ѡ���ģ�ͱ��
	int sel_y2p1_last;//��һʱ��ѡ���ģ�ͱ��
	int sel_y2p1_transtime;//ģ���л��ۼ�ʱ��==========��===============
	int Num_y2p2;//Ԥ����������õ�Ԥ��ģ�ͱ��
	int sel_y2p2_now;//��ǰʱ��ѡ���ģ�ͱ��
	int sel_y2p2_last;//��һʱ��ѡ���ģ�ͱ��
	int sel_y2p2_transtime;//ģ���л��ۼ�ʱ��==========��===============

	int Num_v1p1;//���Ų��õ�Ԥ��ģ�ͱ�� ============��������·���� MR============
	int sel_v1p1_now;//���ŵ�ǰʱ��ѡ���ģ�ͱ��
	int sel_v1p1_last;//������һʱ��ѡ���ģ�ͱ��
	int sel_v1p1_transtime;//����ģ���л��ۼ�ʱ��
	int Num_v1p2;//���Ų��õ�Ԥ��ģ�ͱ�� ============��������·���� MR============
	int sel_v1p2_now;//���ŵ�ǰʱ��ѡ���ģ�ͱ��
	int sel_v1p2_last;//������һʱ��ѡ���ģ�ͱ��
	int sel_v1p2_transtime;//����ģ���л��ۼ�ʱ��

	int Num_v2p1;//���Ų��õ�Ԥ��ģ�ͱ��
	int sel_v2p1_now;//���ŵ�ǰʱ��ѡ���ģ�ͱ��
	int sel_v2p1_last;//������һʱ��ѡ���ģ�ͱ��
	int sel_v2p1_transtime;//����ģ���л��ۼ�ʱ��
	int Num_v2p2;//���Ų��õ�Ԥ��ģ�ͱ��
	int sel_v2p2_now;//���ŵ�ǰʱ��ѡ���ģ�ͱ��
	int sel_v2p2_last;//������һʱ��ѡ���ģ�ͱ��
	int sel_v2p2_transtime;//����ģ���л��ۼ�ʱ��

	int Num_v3p1;//���Ų��õ�Ԥ��ģ�ͱ��
	int sel_v3p1_now;//���ŵ�ǰʱ��ѡ���ģ�ͱ��
	int sel_v3p1_last;//������һʱ��ѡ���ģ�ͱ��
	int sel_v3p1_transtime;//����ģ���л��ۼ�ʱ��
	int Num_v3p2;//���Ų��õ�Ԥ��ģ�ͱ��
	int sel_v3p2_now;//���ŵ�ǰʱ��ѡ���ģ�ͱ��
	int sel_v3p2_last;//������һʱ��ѡ���ģ�ͱ��
	int sel_v3p2_transtime;//����ģ���л��ۼ�ʱ��

	unsigned int AlgState;
	int DataLen;//�������ݳ���
	int initflag;//��ʼ����־
    double erro_p1_now;
	double erro_p2_now;
	int transferflag;//ģ���ض���־λ����1�ض���
	bool ModelUpdate;//ģ���ض�����
	double GpcMatrix_y1mnow;//MGPC+��ǰ������ģ�����
	double GpcMatrix_y2mnow;//MGPC+��ǰ������ģ�����
	int MainDropFlag;

};
