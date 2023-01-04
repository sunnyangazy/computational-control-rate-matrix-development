/*============================================================================================================
�㷨�������Ԥ������㷨
���ܣ���ȡ�����ļ��еı��ض���Ķ������ģ�ͣ�������ɢʱ����ɢ����
	  Ӧ�ý���ʽ����Ԥ������㷨����������������������ѡ�����ģ��
���룺AlgStIn���������룩,Ulast����һʱ�̿�������,SP���趨ֵ��,PV������ֵ��,ConditionIn����ǰ���й���������
�����AlgStOut���㷨״̬�����,U�������������,Udelta�����������������,CondNumSel����ǰѡ��Ĺ����ţ�
���������ڲ������㷨ģ����룬ģ�Ͳ������ļ�����
===============================================================================================================*/
#include "stdafx.h"
#include "AlgMIMOGPC.h"
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include "time.h"
#include "stdlib.h"
#define debug 0
#define PRINTTIME 0


using namespace std;


bool AlgICSMIMOGPC::GpcMatrixCalculate()
{
	int ModelNum;	
	int i,j,flag=1;	
	GpcMatrixClear();//�����ڴ�֮ǰ���ͷ�֮ǰ�ѷ����ڴ�
	for (ModelNum=0;ModelNum<=m_y1p1_ModelPara.totalmodel;ModelNum++)
	{
		int d_selsctmin;
		int d_selsctmax;
		int dan_selsctmin;
		int dan_selsctmax;
		d_selsctmin = m_y1p1_ModelPara.d[ModelNum];
		d_selsctmax = m_y1p1_ModelPara.d[ModelNum];
		dan_selsctmin=m_y2p1_ModelPara.d[ModelNum];
		dan_selsctmax=m_y2p1_ModelPara.d[ModelNum];
		if(c_y1p1_ModelPara.nb2!=0)
		{	
			d_selsctmin = min(m_y1p1_ModelPara.d[ModelNum],m_y1p1_ModelPara.d2[ModelNum]);
			d_selsctmax = max(m_y1p1_ModelPara.d[ModelNum],m_y1p1_ModelPara.d2[ModelNum]);
		}
		if(c_y2p1_ModelPara.nb2!=0)
		{	
			dan_selsctmin = min(m_y2p1_ModelPara.d[ModelNum],m_y2p1_ModelPara.d2[ModelNum]);
			dan_selsctmax = max(m_y2p1_ModelPara.d[ModelNum],m_y2p1_ModelPara.d2[ModelNum]);
		}
		if(malloc_nm(m_TunePara.N-d_selsctmin+1,m_TunePara.Nu1,sizeof(double),sizeof(double*),(double***)&m_y1p1_GpcMatrix[ModelNum].G))
			matrix_init(m_y1p1_GpcMatrix[ModelNum].G,m_TunePara.N-d_selsctmin+1,m_TunePara.Nu1,0);
		else
		{
			flag=0;
			break;
		}
		if(malloc_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y1p1_ModelPara.nb+d_selsctmax-2,sizeof(double),sizeof(double*),(double***)&m_y1p1_GpcMatrix[ModelNum].G1))
			matrix_init(m_y1p1_GpcMatrix[ModelNum].G1,m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y1p1_ModelPara.nb+d_selsctmax-2,0);
		else
		{
			flag=0;
			break;
		}
		if(malloc_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y1p1_ModelPara.na,sizeof(double),sizeof(double*),(double***)&m_y1p1_GpcMatrix[ModelNum].F))
			matrix_init(m_y1p1_GpcMatrix[ModelNum].F,m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y1p1_ModelPara.na,0);
		else
		{
			flag=0;
			break;
		}
		m_y1p1_GpcMatrix[ModelNum].IsAlloc=true;//��̬�����ڴ�ɹ���
		//�������ֵ
		//gpcmatrix(m_y1p1_ModelPara.a[ModelNum],m_y1p1_ModelPara.na,m_y1p1_ModelPara.b[ModelNum],m_y1p1_ModelPara.nb,m_y1p1_ModelPara.d[ModelNum],m_TunePara.N,m_TunePara.Nu1,m_y1p1_GpcMatrix[ModelNum].G,m_y1p1_GpcMatrix[ModelNum].G1,m_y1p1_GpcMatrix[ModelNum].F);
		if(c_y1p1_ModelPara.nb2!=0)
		{
			gpcmatrix_2(m_y1p1_ModelPara.a[ModelNum],m_y1p1_ModelPara.na,m_y1p1_ModelPara.b[ModelNum],m_y1p1_ModelPara.nb,m_y1p1_ModelPara.d[ModelNum],m_y1p1_ModelPara.b2[ModelNum],m_y1p1_ModelPara.nb2,m_y1p1_ModelPara.d2[ModelNum],dan_selsctmin,m_TunePara.N,m_TunePara.Nu1,m_y1p1_GpcMatrix[ModelNum].G,m_y1p1_GpcMatrix[ModelNum].G1,m_y1p1_GpcMatrix[ModelNum].F);
		}
		else
		{
			gpcmatrix(m_y1p1_ModelPara.a[ModelNum],m_y1p1_ModelPara.na,m_y1p1_ModelPara.b[ModelNum],m_y1p1_ModelPara.nb,m_y1p1_ModelPara.d[ModelNum],m_y2p1_ModelPara.d[ModelNum],m_TunePara.N,m_TunePara.Nu1,m_y1p1_GpcMatrix[ModelNum].G,m_y1p1_GpcMatrix[ModelNum].G1,m_y1p1_GpcMatrix[ModelNum].F);
		}
		//cout << "end..." << endl;
	}


	for (ModelNum=0;ModelNum<=m_y1p2_ModelPara.totalmodel;ModelNum++)
	{
		int d_selsctmin;
		int d_selsctmax;
		int dan_selsctmin;
		int dan_selsctmax;
		d_selsctmin = m_y1p2_ModelPara.d[ModelNum];
		d_selsctmax = m_y1p2_ModelPara.d[ModelNum];
		dan_selsctmin=m_y2p2_ModelPara.d[ModelNum];
		dan_selsctmax=m_y2p2_ModelPara.d[ModelNum];
		if(c_y1p2_ModelPara.nb2!=0)
		{	
			d_selsctmin = min(m_y1p2_ModelPara.d[ModelNum],m_y1p2_ModelPara.d2[ModelNum]);
			d_selsctmax = max(m_y1p2_ModelPara.d[ModelNum],m_y1p2_ModelPara.d2[ModelNum]);
		}
		if(c_y2p2_ModelPara.nb2!=0)
		{	
			dan_selsctmin = min(m_y2p2_ModelPara.d[ModelNum],m_y2p2_ModelPara.d2[ModelNum]);
			dan_selsctmax = max(m_y2p2_ModelPara.d[ModelNum],m_y2p2_ModelPara.d2[ModelNum]);
		}
		if(malloc_nm(m_TunePara.N-d_selsctmin+1,m_TunePara.Nu1,sizeof(double),sizeof(double*),(double***)&m_y1p2_GpcMatrix[ModelNum].G))
			matrix_init(m_y1p2_GpcMatrix[ModelNum].G,m_TunePara.N-d_selsctmin+1,m_TunePara.Nu1,0);
		else
		{
			flag=0;
			break;
		}
		if(malloc_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y1p2_ModelPara.nb+d_selsctmax-2,sizeof(double),sizeof(double*),(double***)&m_y1p2_GpcMatrix[ModelNum].G1))
			matrix_init(m_y1p2_GpcMatrix[ModelNum].G1,m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y1p2_ModelPara.nb+d_selsctmax-2,0);
		else
		{
			flag=0;
			break;
		}
		if(malloc_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y1p2_ModelPara.na,sizeof(double),sizeof(double*),(double***)&m_y1p2_GpcMatrix[ModelNum].F))
			matrix_init(m_y1p2_GpcMatrix[ModelNum].F,m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y1p2_ModelPara.na,0);
		else
		{
			flag=0;
			break;
		}
		m_y1p2_GpcMatrix[ModelNum].IsAlloc=true;//��̬�����ڴ�ɹ���
		//�������ֵ
		//gpcmatrix(m_y1p2_ModelPara.a[ModelNum],m_y1p2_ModelPara.na,m_y1p2_ModelPara.b[ModelNum],m_y1p2_ModelPara.nb,m_y1p2_ModelPara.d[ModelNum],m_TunePara.N,m_TunePara.Nu1,m_y1p2_GpcMatrix[ModelNum].G,m_y1p2_GpcMatrix[ModelNum].G1,m_y1p2_GpcMatrix[ModelNum].F);
		if(c_y1p2_ModelPara.nb2!=0)
		{
			gpcmatrix_2(m_y1p2_ModelPara.a[ModelNum],m_y1p2_ModelPara.na,m_y1p2_ModelPara.b[ModelNum],m_y1p2_ModelPara.nb,m_y1p2_ModelPara.d[ModelNum],m_y1p2_ModelPara.b2[ModelNum],m_y1p2_ModelPara.nb2,m_y1p2_ModelPara.d2[ModelNum],dan_selsctmin,m_TunePara.N,m_TunePara.Nu1,m_y1p2_GpcMatrix[ModelNum].G,m_y1p2_GpcMatrix[ModelNum].G1,m_y1p2_GpcMatrix[ModelNum].F);
		}
		else
		{
			gpcmatrix(m_y1p2_ModelPara.a[ModelNum],m_y1p2_ModelPara.na,m_y1p2_ModelPara.b[ModelNum],m_y1p2_ModelPara.nb,m_y1p2_ModelPara.d[ModelNum],m_y2p2_ModelPara.d[ModelNum],m_TunePara.N,m_TunePara.Nu1,m_y1p2_GpcMatrix[ModelNum].G,m_y1p2_GpcMatrix[ModelNum].G1,m_y1p2_GpcMatrix[ModelNum].F);
		}
	}


	for(ModelNum=0;ModelNum<=m_y2p1_ModelPara.totalmodel;ModelNum++)
	{
		int d_selsctmin;
		int d_selsctmax;
		int dan_selsctmin;
		int dan_selsctmax;
		d_selsctmin = m_y2p1_ModelPara.d[ModelNum];
		d_selsctmax = m_y2p1_ModelPara.d[ModelNum];
		dan_selsctmin=m_y1p1_ModelPara.d[ModelNum];
		dan_selsctmax=m_y1p1_ModelPara.d[ModelNum];
		if(c_y2p1_ModelPara.nb2!=0)
		{	
			d_selsctmin = min(m_y2p1_ModelPara.d[ModelNum],m_y2p1_ModelPara.d2[ModelNum]);
			d_selsctmax = max(m_y2p1_ModelPara.d[ModelNum],m_y2p1_ModelPara.d2[ModelNum]);
		}
		if(c_y1p1_ModelPara.nb2!=0)
		{	
			dan_selsctmin = min(m_y1p1_ModelPara.d[ModelNum],m_y1p1_ModelPara.d2[ModelNum]);
			dan_selsctmax = max(m_y1p1_ModelPara.d[ModelNum],m_y1p1_ModelPara.d2[ModelNum]);
		}
		if(malloc_nm(m_TunePara.N-d_selsctmin+1,m_TunePara.Nu2,sizeof(double),sizeof(double*),(double***)&m_y2p1_GpcMatrix[ModelNum].G))
			matrix_init(m_y2p1_GpcMatrix[ModelNum].G,m_TunePara.N-d_selsctmin+1,m_TunePara.Nu2,0);
		else
		{
			flag=0;
			break;
		}
		if(malloc_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y2p1_ModelPara.nb+d_selsctmax-2,sizeof(double),sizeof(double*),(double***)&m_y2p1_GpcMatrix[ModelNum].G1))
			matrix_init(m_y2p1_GpcMatrix[ModelNum].G1,m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y2p1_ModelPara.nb+d_selsctmax-2,0);
		else
		{
			flag=0;
			break;
		}
		if(malloc_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y2p1_ModelPara.na,sizeof(double),sizeof(double*),(double***)&m_y2p1_GpcMatrix[ModelNum].F))
			matrix_init(m_y2p1_GpcMatrix[ModelNum].F,m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y2p1_ModelPara.na,0);
		else
		{
			flag=0;
			break;
		}
		m_y2p1_GpcMatrix[ModelNum].IsAlloc=true;//��̬�����ڴ�ɹ���
		//�������ֵ
		//gpcmatrix(m_y2p1_ModelPara.a[ModelNum],m_y2p1_ModelPara.na,m_y2p1_ModelPara.b[ModelNum],m_y2p1_ModelPara.nb,m_y2p1_ModelPara.d[ModelNum],m_TunePara.N,m_TunePara.Nu2,m_y2p1_GpcMatrix[ModelNum].G,m_y2p1_GpcMatrix[ModelNum].G1,m_y2p1_GpcMatrix[ModelNum].F);
		if(c_y2p1_ModelPara.nb2!=0)
		{
			gpcmatrix_2(m_y2p1_ModelPara.a[ModelNum],m_y2p1_ModelPara.na,m_y2p1_ModelPara.b[ModelNum],m_y2p1_ModelPara.nb,m_y2p1_ModelPara.d[ModelNum],m_y2p1_ModelPara.b2[ModelNum],m_y2p1_ModelPara.nb2,m_y2p1_ModelPara.d2[ModelNum],dan_selsctmin,m_TunePara.N,m_TunePara.Nu2,m_y2p1_GpcMatrix[ModelNum].G,m_y2p1_GpcMatrix[ModelNum].G1,m_y2p1_GpcMatrix[ModelNum].F);
		}
		else
		{
			gpcmatrix(m_y2p1_ModelPara.a[ModelNum],m_y2p1_ModelPara.na,m_y2p1_ModelPara.b[ModelNum],m_y2p1_ModelPara.nb,m_y2p1_ModelPara.d[ModelNum],m_y1p1_ModelPara.d[ModelNum],m_TunePara.N,m_TunePara.Nu2,m_y2p1_GpcMatrix[ModelNum].G,m_y2p1_GpcMatrix[ModelNum].G1,m_y2p1_GpcMatrix[ModelNum].F);
		}
	}

	for (ModelNum=0;ModelNum<=m_y2p2_ModelPara.totalmodel;ModelNum++)
	{
		int d_selsctmin;
		int d_selsctmax;
		int dan_selsctmin;
		int dan_selsctmax;
		d_selsctmin = m_y2p2_ModelPara.d[ModelNum];
		d_selsctmax = m_y2p2_ModelPara.d[ModelNum];
		dan_selsctmin=m_y1p2_ModelPara.d[ModelNum];
		dan_selsctmax=m_y1p2_ModelPara.d[ModelNum];
		if(c_y2p2_ModelPara.nb2!=0)
		{	
			d_selsctmin = min(m_y2p2_ModelPara.d[ModelNum],m_y2p2_ModelPara.d2[ModelNum]);
			d_selsctmax = max(m_y2p2_ModelPara.d[ModelNum],m_y2p2_ModelPara.d2[ModelNum]);
		}
		if(c_y1p2_ModelPara.nb2!=0)
		{	
			dan_selsctmin = min(m_y1p2_ModelPara.d[ModelNum],m_y1p2_ModelPara.d2[ModelNum]);
			dan_selsctmax = max(m_y1p2_ModelPara.d[ModelNum],m_y1p2_ModelPara.d2[ModelNum]);
		}
		
		if(malloc_nm(m_TunePara.N-d_selsctmin+1,m_TunePara.Nu2,sizeof(double),sizeof(double*),(double***)&m_y2p2_GpcMatrix[ModelNum].G))
			matrix_init(m_y2p2_GpcMatrix[ModelNum].G,m_TunePara.N-d_selsctmin+1,m_TunePara.Nu2,0);
		else
		{
			flag=0;
			break;
		}
		if(malloc_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y2p2_ModelPara.nb+d_selsctmax-2,sizeof(double),sizeof(double*),(double***)&m_y2p2_GpcMatrix[ModelNum].G1))
			matrix_init(m_y2p2_GpcMatrix[ModelNum].G1,m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y2p2_ModelPara.nb+d_selsctmax-2,0);
		else
		{
			flag=0;
			break;
		}
		if(malloc_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y2p2_ModelPara.na,sizeof(double),sizeof(double*),(double***)&m_y2p2_GpcMatrix[ModelNum].F))
			matrix_init(m_y2p2_GpcMatrix[ModelNum].F,m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,m_y2p2_ModelPara.na,0);
		else
		{
			flag=0;
			break;
		}
		m_y2p2_GpcMatrix[ModelNum].IsAlloc=true;//��̬�����ڴ�ɹ���
		//�������ֵ
		//gpcmatrix(m_y2p2_ModelPara.a[ModelNum],m_y2p2_ModelPara.na,m_y2p2_ModelPara.b[ModelNum],m_y2p2_ModelPara.nb,m_y2p2_ModelPara.d[ModelNum],m_TunePara.N,m_TunePara.Nu2,m_y2p2_GpcMatrix[ModelNum].G,m_y2p2_GpcMatrix[ModelNum].G1,m_y2p2_GpcMatrix[ModelNum].F);
		if(c_y2p2_ModelPara.nb2!=0)
		{
			gpcmatrix_2(m_y2p2_ModelPara.a[ModelNum],m_y2p2_ModelPara.na,m_y2p2_ModelPara.b[ModelNum],m_y2p2_ModelPara.nb,m_y2p2_ModelPara.d[ModelNum],m_y2p2_ModelPara.b2[ModelNum],m_y2p2_ModelPara.nb2,m_y2p2_ModelPara.d2[ModelNum],dan_selsctmin,m_TunePara.N,m_TunePara.Nu2,m_y2p2_GpcMatrix[ModelNum].G,m_y2p2_GpcMatrix[ModelNum].G1,m_y2p2_GpcMatrix[ModelNum].F);
		}
		else
		{
			gpcmatrix(m_y2p2_ModelPara.a[ModelNum],m_y2p2_ModelPara.na,m_y2p2_ModelPara.b[ModelNum],m_y2p2_ModelPara.nb,m_y2p2_ModelPara.d[ModelNum],m_y1p2_ModelPara.d[ModelNum],m_TunePara.N,m_TunePara.Nu2,m_y2p2_GpcMatrix[ModelNum].G,m_y2p2_GpcMatrix[ModelNum].G1,m_y2p2_GpcMatrix[ModelNum].F);
		}
	}


	for(ModelNum=0;ModelNum<=Calculate_max(m_y1p1_ModelPara.totalmodel,m_y2p1_ModelPara.totalmodel);ModelNum++)
	{
		int d_selsctmin;
		d_selsctmin = m_y1p1_ModelPara.d[ModelNum];
		int dan_selsctmin;
		dan_selsctmin=m_y2p1_ModelPara.d[ModelNum];
		if(c_y1p1_ModelPara.nb2!=0)
		{	
			d_selsctmin = min(m_y1p1_ModelPara.d[ModelNum],m_y1p1_ModelPara.d2[ModelNum]);
		}
		if(c_y2p1_ModelPara.nb2!=0)
		{	
			dan_selsctmin = min(m_y2p1_ModelPara.d[ModelNum],m_y2p1_ModelPara.d2[ModelNum]);
		}
		if(malloc_nm(m_TunePara.N-Calculate_min(d_selsctmin,dan_selsctmin)+1,m_TunePara.Nu1+m_TunePara.Nu2,sizeof(double),sizeof(double*),(double***)&m_p1_OCR[ModelNum].OCR))
			matrix_init(m_p1_OCR[ModelNum].OCR,m_TunePara.N-Calculate_min(d_selsctmin,dan_selsctmin)+1,m_TunePara.Nu1+m_TunePara.Nu2,0);
		else
		{
			flag=0;
			break;
		}
		m_p1_OCR[ModelNum].IsAlloc=true;//��̬�����ڴ�ɹ���
	}
	for(ModelNum=0;ModelNum<=Calculate_max(m_y1p2_ModelPara.totalmodel,m_y2p2_ModelPara.totalmodel);ModelNum++)
	{
		int d_selsctmin;
		d_selsctmin = m_y1p2_ModelPara.d[ModelNum];
		int dan_selsctmin;
		dan_selsctmin=m_y2p2_ModelPara.d[ModelNum];
		if(c_y1p2_ModelPara.nb2!=0)
		{	
			d_selsctmin = min(m_y1p2_ModelPara.d[ModelNum],m_y1p2_ModelPara.d2[ModelNum]);
		}
		if(c_y2p2_ModelPara.nb2!=0)
		{	
			dan_selsctmin = min(m_y2p2_ModelPara.d[ModelNum],m_y2p2_ModelPara.d2[ModelNum]);
		}
		if(malloc_nm(m_TunePara.N-Calculate_min(d_selsctmin,dan_selsctmin)+1,m_TunePara.Nu1+m_TunePara.Nu2,sizeof(double),sizeof(double*),(double***)&m_p2_OCR[ModelNum].OCR))
			matrix_init(m_p2_OCR[ModelNum].OCR,m_TunePara.N-Calculate_min(d_selsctmin,dan_selsctmin)+1,m_TunePara.Nu1+m_TunePara.Nu2,0);
		else
		{
			flag=0;
			break;
		}
		m_p2_OCR[ModelNum].IsAlloc=true;//��̬�����ڴ�ɹ���
	}

/*	double** A0,**B0,**C0;	
	malloc_nm(4,2,sizeof(double),sizeof(double*),(double***)&A0);
	malloc_nm(4,3,sizeof(double),sizeof(double*),(double***)&B0);
	malloc_nm(4,5,sizeof(double),sizeof(double*),(double***)&C0);
	int ia=0;
	int ib=0;
	int ic=0;
	int ja=0;
	int jb=0;
	int jc=0;
	for(ia=0;ia<4;ia++)
	{
		for(ja=0;ja<2;ja++)
		{	
			A0[ia][ja]=1;
			printf("%f ",A0[ia][ja]);
		}
		printf("\n");
	}
	printf("\n");
	printf("\n");
	for(ib=0;ib<4;ib++)
	{
		for(jb=0;jb<3;jb++)
		{	
			B0[ib][jb]=1;
			printf("%f ",B0[ib][jb]);
		}
		printf("\n");
	}
	printf("\n");
	printf("\n");
	for(ic=0;ic<4;ic++)
	{
		for(jc=0;jc<5;jc++)
		{	
			C0[ic][jc]=0;
			printf("%f ",C0[ic][jc]);
		}
		printf("\n");
	}
	printf("\n");
	printf("\n");
	MatrixCombine(B0,A0,4,3,4,2,C0);
	for(ic=0;ic<4;ic++)
	{
		for(jc=0;jc<5;jc++)
		{	
			printf("%f ",C0[ic][jc]);
		}
		printf("\n");
	}
*/	
	for(int ModNum = 0;ModNum<=Calculate_max(m_y1p1_ModelPara.totalmodel,m_y2p1_ModelPara.totalmodel);ModNum++)
	{
		int d_selsctmin;
		d_selsctmin = m_y1p1_ModelPara.d[ModNum];
		int dan_selsctmin;
		dan_selsctmin=m_y2p1_ModelPara.d[ModNum];
		if(c_y1p1_ModelPara.nb2!=0)
		{	
			d_selsctmin = min(m_y1p1_ModelPara.d[ModNum],m_y1p1_ModelPara.d2[ModNum]);
		}
		if(c_y2p1_ModelPara.nb2!=0)
		{	
			dan_selsctmin = min(m_y2p1_ModelPara.d[ModNum],m_y2p1_ModelPara.d2[ModNum]);
		}
		MatrixCombine(m_y1p1_GpcMatrix[ModNum].G,m_y2p1_GpcMatrix[ModNum].G,m_TunePara.N-d_selsctmin+1,m_TunePara.Nu1 ,m_TunePara.N-dan_selsctmin+1,m_TunePara.Nu2 ,m_p1_OCR[ModNum].OCR);
	}
	for(int ModNum = 0;ModNum<=Calculate_max(m_y1p2_ModelPara.totalmodel,m_y2p2_ModelPara.totalmodel);ModNum++)
	{
		int d_selsctmin;
		d_selsctmin = m_y1p2_ModelPara.d[ModNum];
		int dan_selsctmin;
		dan_selsctmin=m_y2p2_ModelPara.d[ModNum];
		if(c_y1p2_ModelPara.nb2!=0)
		{	
			d_selsctmin = min(m_y1p2_ModelPara.d[ModNum],m_y1p2_ModelPara.d2[ModNum]);
		}
		if(c_y2p2_ModelPara.nb2!=0)
		{	
			dan_selsctmin = min(m_y2p2_ModelPara.d[ModNum],m_y2p2_ModelPara.d2[ModNum]);
		}
		MatrixCombine(m_y1p2_GpcMatrix[ModNum].G,m_y2p2_GpcMatrix[ModNum].G,m_TunePara.N-d_selsctmin+1,m_TunePara.Nu1 ,m_TunePara.N-dan_selsctmin+1,m_TunePara.Nu2 ,m_p2_OCR[ModNum].OCR);
	}


	for(int ModNum = 0;ModNum<=Calculate_max(m_y1p1_ModelPara.totalmodel,m_y2p1_ModelPara.totalmodel);ModNum++)
	{
		int d_selsctmin;
		d_selsctmin = m_y1p1_ModelPara.d[ModNum];
		int dan_selsctmin;
		dan_selsctmin=m_y2p1_ModelPara.d[ModNum];
		if(c_y1p1_ModelPara.nb2!=0)
		{	
			d_selsctmin = min(m_y1p1_ModelPara.d[ModNum],m_y1p1_ModelPara.d2[ModNum]);
		}
		if(c_y2p1_ModelPara.nb2!=0)
		{	
			dan_selsctmin = min(m_y2p1_ModelPara.d[ModNum],m_y2p1_ModelPara.d2[ModNum]);
		}
		m_p1_OCR[ModNum].GR.resize((m_TunePara.Nu1+m_TunePara.Nu2)*(m_TunePara.N-Calculate_min(d_selsctmin,dan_selsctmin)+1),0);	
	}
	for(int ModNum = 0;ModNum<=Calculate_max(m_y1p2_ModelPara.totalmodel,m_y2p2_ModelPara.totalmodel);ModNum++)
	{
		int d_selsctmin;
		d_selsctmin = m_y1p2_ModelPara.d[ModNum];
		int dan_selsctmin;
		dan_selsctmin=m_y2p2_ModelPara.d[ModNum];
		if(c_y1p2_ModelPara.nb2!=0)
		{	
			d_selsctmin = min(m_y1p2_ModelPara.d[ModNum],m_y1p2_ModelPara.d2[ModNum]);
		}
		if(c_y2p2_ModelPara.nb2!=0)
		{	
			dan_selsctmin = min(m_y2p2_ModelPara.d[ModNum],m_y2p2_ModelPara.d2[ModNum]);
		}
		m_p2_OCR[ModNum].GR.resize((m_TunePara.Nu1+m_TunePara.Nu2)*(m_TunePara.N-Calculate_min(d_selsctmin,dan_selsctmin)+1),0);
			
	}
	for(int ModNum = 0;ModNum<=m_y1p1_ModelPara.totalmodel;ModNum++)
	{
		int d_y1p1_selsctmin;
		d_y1p1_selsctmin = m_y1p1_ModelPara.d[ModNum];
		int d_y1p2_selsctmin;
		d_y1p2_selsctmin = m_y1p2_ModelPara.d[ModNum];
		int d_y2p1_selsctmin;
		d_y2p1_selsctmin = m_y2p1_ModelPara.d[ModNum];
		int d_y2p2_selsctmin;
		d_y2p2_selsctmin = m_y2p2_ModelPara.d[ModNum];
		if(c_y1p1_ModelPara.nb2!=0)
		{	
			d_y1p1_selsctmin = min(m_y1p1_ModelPara.d[ModNum],m_y1p1_ModelPara.d2[ModNum]);
		}
		if(c_y1p2_ModelPara.nb2!=0)
		{	
			d_y1p2_selsctmin = min(m_y1p2_ModelPara.d[ModNum],m_y1p2_ModelPara.d2[ModNum]);
		}
		if(c_y2p1_ModelPara.nb2!=0)
		{	
			d_y2p1_selsctmin = min(m_y2p1_ModelPara.d[ModNum],m_y2p1_ModelPara.d2[ModNum]);
		}
		if(c_y2p2_ModelPara.nb2!=0)
		{	
			d_y2p2_selsctmin = min(m_y2p2_ModelPara.d[ModNum],m_y2p2_ModelPara.d2[ModNum]);
		}

		CRateCalculate(m_p1_OCR[ModNum].OCR,m_p2_OCR[ModNum].OCR,m_p1_OCR[ModNum].GR,m_p2_OCR[ModNum].GR,(m_TunePara.N-Calculate_min(d_y1p1_selsctmin,d_y2p1_selsctmin)+1),(m_TunePara.N-Calculate_min(d_y1p2_selsctmin,d_y2p2_selsctmin)+1),m_TunePara.Nu1+m_TunePara.Nu2,m_TunePara.gama1,m_TunePara.gama2,m_TunePara.eta1,m_TunePara.eta2);
	}


	if(flag==1)	
		return true;
	else
		return false;
}


//�ͷŶ�̬������ڴ�
void AlgICSMIMOGPC::GpcMatrixClear()
{
	int ModelNum;
	for (ModelNum=0;ModelNum<MIMOGPC_MODELNUM_MAX;ModelNum++)
	{
		//cout << "begin to free :" << endl;
		//cout << "begin GpcMatrixClear()!" << endl;
		if(m_y1p1_GpcMatrix[ModelNum].IsAlloc==true)
		{
			int d_selsctmin;
			int dan_selsctmin;
			d_selsctmin = m_y1p1_ModelPara.d[ModelNum];
			dan_selsctmin=m_y2p1_ModelPara.d[ModelNum];
			if(c_y1p1_ModelPara.nb2!=0)
			{	
				d_selsctmin = min(m_y1p1_ModelPara.d[ModelNum],m_y1p1_ModelPara.d2[ModelNum]);
			}
			if(c_y2p1_ModelPara.nb2!=0)
			{	
				dan_selsctmin = min(m_y2p1_ModelPara.d[ModelNum],m_y2p1_ModelPara.d2[ModelNum]);
			}
			free_nm(m_TunePara.N-d_selsctmin+1,(double**)m_y1p1_GpcMatrix[ModelNum].G);
			free_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,(double**)m_y1p1_GpcMatrix[ModelNum].G1);
			free_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,(double**)m_y1p1_GpcMatrix[ModelNum].F);
		}
		m_y1p1_GpcMatrix[ModelNum].IsAlloc=false;

		if(m_y1p2_GpcMatrix[ModelNum].IsAlloc==true)
		{
			int d_selsctmin;
			int dan_selsctmin;
			d_selsctmin = m_y1p2_ModelPara.d[ModelNum];
			dan_selsctmin=m_y2p2_ModelPara.d[ModelNum];
			if(c_y1p2_ModelPara.nb2!=0)
			{	
				d_selsctmin = min(m_y1p2_ModelPara.d[ModelNum],m_y1p2_ModelPara.d2[ModelNum]);
			}
			if(c_y2p2_ModelPara.nb2!=0)
			{	
				dan_selsctmin = min(m_y2p2_ModelPara.d[ModelNum],m_y2p2_ModelPara.d2[ModelNum]);
			}
			free_nm(m_TunePara.N-d_selsctmin+1,(double**)m_y1p2_GpcMatrix[ModelNum].G);
			free_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,(double**)m_y1p2_GpcMatrix[ModelNum].G1);
			free_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,(double**)m_y1p2_GpcMatrix[ModelNum].F);
		}
		m_y1p2_GpcMatrix[ModelNum].IsAlloc=false;

		if(m_y2p1_GpcMatrix[ModelNum].IsAlloc==true)
		{
			int d_selsctmin;
			int dan_selsctmin;
			d_selsctmin = m_y2p1_ModelPara.d[ModelNum];
			dan_selsctmin=m_y1p1_ModelPara.d[ModelNum];
			if(c_y2p1_ModelPara.nb2!=0)
			{	
				d_selsctmin = min(m_y2p1_ModelPara.d[ModelNum],m_y2p1_ModelPara.d2[ModelNum]);
			}
			if(c_y1p1_ModelPara.nb2!=0)
			{	
				dan_selsctmin = min(m_y1p1_ModelPara.d[ModelNum],m_y1p1_ModelPara.d2[ModelNum]);
			}
			free_nm(m_TunePara.N-d_selsctmin+1,(double**)m_y2p1_GpcMatrix[ModelNum].G);
			free_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,(double**)m_y2p1_GpcMatrix[ModelNum].G1);
			free_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,(double**)m_y2p1_GpcMatrix[ModelNum].F);
		}
		m_y2p1_GpcMatrix[ModelNum].IsAlloc=false;

		if(m_y2p2_GpcMatrix[ModelNum].IsAlloc==true)
		{
			int d_selsctmin;
			int dan_selsctmin;
			d_selsctmin = m_y2p2_ModelPara.d[ModelNum];
			dan_selsctmin=m_y1p2_ModelPara.d[ModelNum];
			if(c_y2p2_ModelPara.nb2!=0)
			{	
				d_selsctmin = min(m_y2p2_ModelPara.d[ModelNum],m_y2p2_ModelPara.d2[ModelNum]);
			}
			if(c_y1p2_ModelPara.nb2!=0)
			{	
				dan_selsctmin = min(m_y1p2_ModelPara.d[ModelNum],m_y1p2_ModelPara.d2[ModelNum]);
			}
			free_nm(m_TunePara.N-d_selsctmin+1,(double**)m_y2p2_GpcMatrix[ModelNum].G);
			free_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,(double**)m_y2p2_GpcMatrix[ModelNum].G1);
			free_nm(m_TunePara.N-min(d_selsctmin,dan_selsctmin)+1,(double**)m_y2p2_GpcMatrix[ModelNum].F);
		}
		m_y2p2_GpcMatrix[ModelNum].IsAlloc=false;
	
		if(m_p1_OCR[ModelNum].IsAlloc==true)
		{
			int d_selsctmin;
			int dan_selsctmin;
			d_selsctmin = m_y1p1_ModelPara.d[ModelNum];
			dan_selsctmin=m_y2p1_ModelPara.d[ModelNum];
			if(c_y1p1_ModelPara.nb2!=0)
			{	
				d_selsctmin = min(m_y1p1_ModelPara.d[ModelNum],m_y1p1_ModelPara.d2[ModelNum]);
			}
			if(c_y2p1_ModelPara.nb2!=0)
			{	
				dan_selsctmin = min(m_y2p1_ModelPara.d[ModelNum],m_y2p1_ModelPara.d2[ModelNum]);
			}
			free_nm(m_TunePara.N-Calculate_min(d_selsctmin,dan_selsctmin)+1,(double**)m_p1_OCR[ModelNum].OCR);
		}
		m_p1_OCR[ModelNum].IsAlloc=false;
		
		if(m_p2_OCR[ModelNum].IsAlloc==true)
		{
			int d_selsctmin;
			int dan_selsctmin;
			d_selsctmin = m_y1p2_ModelPara.d[ModelNum];
			dan_selsctmin=m_y2p2_ModelPara.d[ModelNum];
			if(c_y1p2_ModelPara.nb2!=0)
			{	
				d_selsctmin = min(m_y1p2_ModelPara.d[ModelNum],m_y1p2_ModelPara.d2[ModelNum]);
			}
			if(c_y2p2_ModelPara.nb2!=0)
			{	
				dan_selsctmin = min(m_y2p2_ModelPara.d[ModelNum],m_y2p2_ModelPara.d2[ModelNum]);
			}
			free_nm(m_TunePara.N-Calculate_min(d_selsctmin,dan_selsctmin)+1,(double**)m_p2_OCR[ModelNum].OCR);
		}
		m_p2_OCR[ModelNum].IsAlloc=false;
	
	}
	//cout << "GpcMatrixClear() finally complete!" << endl;
}

/*��̬����һ��n��m�����飬size_1=sizeof(type),size_2=sizeof(type*),aΪ��̬����Ķ�ά�����ַ*/
bool AlgICSMIMOGPC::malloc_nm(int n,int m,int size_1,int size_2,double*** a)
{
	if ((*a=(double**)malloc(n*size_2))==NULL)
		return false;
	int i;
	for (i=0;i<n;i++)
		if (((*a)[i]=(double*)malloc(m*size_1))==NULL)
			return false;
		else
			memset((*a)[i],'0',m*size_1);//�����ַ�
	return true;
}

/*�ͷŶ�̬�����һ��n��m������,aΪ��̬����Ķ�ά�����ַ*/
void AlgICSMIMOGPC::free_nm(int n,double** a)
{
	int i;
	for (i=0;i<n;i++)
		free((a)[i]);
	free(a);
	a=0;
}

/*����*/
void AlgICSMIMOGPC::conv(double a[],double b[],int na,int nb,double c[])
{
	int i,j;
	int k=na+nb-1;

	for (i=0;i<k;i++)
		c[i]=0;

	for (i=0;i<k;i++){
		for (j=max(0,i+1-nb);j<=min(i,na-1);j++)
			c[i]+=a[j]*b[i-j];
	}
}

/*%****************************************************************************************************
  %���ܣ��ಽDiophanine���̵����
  %���ø�ʽ��multidiophantine(double a[],int na,double b[],int nb,int N,double*** E,double*** F,double*** G)
  %�������������ʽa��bϵ���Լ��״Σ�Ԥ�ⲽ��
  %���������E��F��GΪDiophanine���̵Ľ⣬����ָ�룬E��N��N�У�F��N��nb-1+N�У�G��N��na-1��
  %���У�na��nbΪ����a��b�ĳ���
%******************************************************************************************************/
void AlgICSMIMOGPC::multidiophantine(double a[],int na,double b[],int nb,int N,double** E,double** F,double** G)
{
	int i,j;

	na=na-1;nb=nb-1;
	/*����ֵ*/
	E[0][0]=1;
	conv(b,E[0],nb+1,N,F[0]);
	for (i=0;i<na;i++)
		//G[0][i+na-1]=-a[i+1];
		G[0][i]=-a[i+1];
	
	for (j=0;j<N-1;j++){
		for (i=0;i<=j;i++)//for (i=0;i<=j;i++)
			E[j+1][i]=E[j][i];
		E[j+1][j+1]=G[j][0];

		for(i=1;i<na;i++)
			G[j+1][i-1]=G[j][i]-G[j][0]*a[i];
	
		G[j+1][na-1]=-G[j][0]*a[na];
		conv(b,E[j+1],nb+1,N,F[j+1]);
	}
}

//��n��m�о���aÿ��ֵ��ʼ��Ϊa0
void AlgICSMIMOGPC::matrix_init(double** a,int n,int m,double a0)
{
	int i,j;
	for (i=0;i<n;i++){
		for (j=0;j<m;j++)
			a[i][j]=a0;
	}	
}

/*%********************************************************************************************************
  %���ܣ�����GPC(����Ԥ�����)�㷨�������(���ض���ģ�ͣ�a*y=b*u+e/��)
  %���ø�ʽ��gpcmatrix(double a[],int na,double b[],int nb,int N,int Nu,int d,double** G,double** G1,double** F)
  %�����������ɢϵͳģ��a��bϵ���������������ȣ�Ԥ�ⲽ��N,���Ʋ���Nu(d>1,N>Nu)
  %���������G(N-d+1,Nu)��G1(Nu-d+1,nb+d-2)��F(Nu-d+1,na-1),
  %��Y(k)=G��U(k)+G1��U(k-1)+FY(k-1),Y(k)Ϊδ��ʱ������Y(k-1)Ϊ��ȥʱ����֪����
%**********************************************************************************************************/
void AlgICSMIMOGPC::gpcmatrix(double a[],int na,double b[],int nb,int d_OR,int d2,int N,int Nu,double** G,double** G1,double** F)
{
	na=na-1;
	int i,j,naa,nbb;
	double *bb,*aa;
	double ta[2]={1,-1};
	
	if (d_OR<1)
	{
		d_OR=1;
	}
	if (d2<1)
	{
		d2=1;
	}

	bb=(double*)malloc((d_OR-1+nb)*sizeof(double));//�ӳٴ���1����bǰ���d-1��0
	for (i=0;i<(d_OR-1);i++){
		bb[i]=0;
	}
	for (i=0;i<(nb);i++){
		bb[d_OR-1+i]=b[i];
	}

	aa=(double*)malloc((na+2)*sizeof(double));
	
	naa=na+1;nbb=nb+d_OR-1;
	conv(a,ta,na+1,2,aa);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
	
	double** E0,**F0,**G0;//������ʱ����
	
	malloc_nm(N,N,sizeof(double),sizeof(double*),(double***)&E0);
	malloc_nm(N,N+nbb-1,sizeof(double),sizeof(double*),(double***)&F0);
	malloc_nm(N,naa,sizeof(double),sizeof(double*),(double***)&G0);
	multidiophantine(aa,naa+1,bb,nbb,N,E0,F0,G0);//��ⶪ��ͼ����
	
	for (i=0;i<N-min(d_OR,d2)+1;i++)
	{//���F����
		for(j=0;j<=na;j++)
		{
			F[i][j]=G0[i+min(d_OR,d2)-1][j];
		}
	}
/*
	for(i=1;i<=N-d_OR+1;i++)
	{
		for (j=1;j<=min(i,Nu);j++)
			G[i-1][j-1]=F0[i+d_OR-2][i+d_OR-j-1];//���G����
		for (j=1;j<=nbb-1;j++)
			G1[i-1][j-1]=F0[i+d_OR-2][i+d_OR-2+j];//���G1����
	}
*/
	for(i=1;i<=N-d_OR+1;i++)
	{
		for (j=1;j<=min(i,Nu);j++)
			G[i-1][j-1]=F0[i+d_OR-2][i+d_OR-j-1];//���G����
	}
	for(i=1;i<=N-min(d_OR,d2)+1;i++)
	{
		for (j=1;j<=nbb-1;j++)
			G1[i-1][j-1]=F0[i+min(d_OR,d2)-2][i+min(d_OR,d2)-2+j];//���G1����
	}

/*
	if(d_OR>d2)
	{
		
		for(i=1;i<=N-d_OR+1;i++)
		{
			for (j=1;j<=nbb-1;j++)
			{
				G1[i+d_OR-d2-1][j-1]=F0[i+d_OR-2][i+d_OR-2+j];//���G1����
			}
				
		}
		for(i=1;i<=d_OR-d2;i++)
		{
			for (j=1;j<=nbb-1;j++)
			{
				G1[d_OR-d2-i][j-1]=F0[d_OR-i-1][d_OR-1-i+j];//���G1����
			}
				
		}
	}
*/

	//cout << "gpcmatrix begin to free!" << endl;
	free_nm(N,(double**)E0);//�ͷ���ʱ����Ķ�̬����
	free_nm(N,(double**)F0);
	free_nm(N,(double**)G0);
	
	free(bb);
	free(aa);
	//cout << "gpcmatrix complete to free!" << endl;
}

void AlgICSMIMOGPC::gpcmatrix_2(double a[],int na,double b1[],int nb1,int d1_OR,double b2[],int nb2,int d2_OR,int d_another_min,int N,int Nu,double** G,double** G1,double** F)
{
	na=na-1;
	int i,j,naa,nbb;
	double *bb,*aa;
	double ta[2]={1,-1};
	
	if (d1_OR<1){
		d1_OR=1;
	}
	if (d2_OR<1){
		d2_OR=1;
	}

	bb=(double*)malloc((max(d1_OR,d2_OR)-1+max(nb1,nb2))*sizeof(double));//�ӳٴ���1����bǰ���d-1��0
	for (i=0;i<(max(d1_OR,d2_OR)-1+max(nb1,nb2));i++)
	{
		bb[i]=0;
	}
	for (i=0;i<(nb1);i++)
	{
		bb[d1_OR-1+i]=b1[i];
	}
	for (i=0;i<(nb2);i++)
	{
		bb[d2_OR-1+i]=bb[d2_OR-1+i]+b2[i];
	}


	aa=(double*)malloc((na+2)*sizeof(double));
	
	naa=na+1;
	nbb=max(d1_OR,d2_OR)-1+max(nb1,nb2);

	conv(a,ta,na+1,2,aa);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
	
	double** E0,**F0,**G0;//������ʱ����
	
	malloc_nm(N,N,sizeof(double),sizeof(double*),(double***)&E0);
	malloc_nm(N,N+nbb-1,sizeof(double),sizeof(double*),(double***)&F0);	
	malloc_nm(N,naa,sizeof(double),sizeof(double*),(double***)&G0);
	multidiophantine(aa,naa+1,bb,nbb,N,E0,F0,G0);
	for (i=0;i<N-min(min(d1_OR,d2_OR),d_another_min)+1;i++)
	{//���F����
		for(j=0;j<=na;j++)
			F[i][j]=G0[i+min(min(d1_OR,d2_OR),d_another_min)-1][j];
	}
	for(i=1;i<=N-min(d1_OR,d2_OR)+1;i++)
	{
		for (j=1;j<=min(i,Nu);j++)
			G[i-1][j-1]=F0[i+min(d1_OR,d2_OR)-2][i+min(d1_OR,d2_OR)-j-1];//���G����
	}
	for(i=1;i<=N-min(min(d1_OR,d2_OR),d_another_min)+1;i++)
	{
		for (j=1;j<=nbb-1;j++)
			G1[i-1][j-1]=F0[i+min(min(d1_OR,d2_OR),d_another_min)-2][i+min(min(d1_OR,d2_OR),d_another_min)-2+j];//���G1����
	}
	//cout << "gpcmatrix begin to free!" << endl;
	free_nm(N,(double**)E0);//�ͷ���ʱ����Ķ�̬����
	free_nm(N,(double**)F0);
	free_nm(N,(double**)G0);
	
	free(bb);
	free(aa);
	//cout << "gpcmatrix complete to free!" << endl;
}




//�������ļ���ȡϵͳģ��
//�ļ���Ϊ��MgpcModel+���+.para;���Ϊ6λ���֣�ǰ��λΪҳ�ţ�������λ��0������λΪID�ţ�������λ��0;
//���磺��5ҳ���㷨IDΪ3�����Ϊ��005003�������ļ���Ϊ��MgpcModel005003.para;
//�ļ�Ŀ¼������Ŀ¼\domain013\data\ctrl\obj\drop001
bool AlgICSMIMOGPC::ReadParay1p1FromFile()//Ҫ�� MR ?
{	
	char filedir[]="\/edpf\/prj\/running\/ctrl\/obj";//���ļ�Ŀ¼
 //   char filedir2[]="E:\\project\\DPU59running\\ctrl\\obj";//������DPU�ļ�Ŀ¼
	char filedir2[]="C:\\Program Files\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
    char filedir3[]="C:\\Program Files (x86)\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
	//char filedir2[]="C:\\Projects\\edpf\\MGPCMatrix";//������DPU�ļ�Ŀ¼
	AlgIdType AlgId,SheetNo,m_RecId,FinalId;
	AlgId=GetId();//��ȡ�㷨Id
	SheetNo=AlgId>>16;
	m_RecId=AlgId&0x000000FF;
	SheetNo=AlgId>>16;

	short int	St;
	char	filename[MAX_PATH]="";
	char	FilePath[MAX_PATH]="";
	char	StrBuf[256];
	size_t	Len=MAX_PATH;
	short int	rc;
	short int	j;

	FILE *fp=NULL;
	char str[MAX_PATH];
	char *p;
	char *delims={" []=\n"};
	// �����ļ���
	sprintf(filename,"Mgpcc1p1Model%03d%03d.para",SheetNo,m_RecId);
	sprintf(FilePath,"%s/%s",filedir,filename);
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(FilePath,"%s\\%s",filedir2,filename);
		if((fp=fopen(FilePath,"r"))==NULL)
		{
			sprintf(FilePath,"%s\\%s",filedir3,filename);
		    if((fp=fopen(FilePath,"r"))==NULL)
		    {
			    sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
			    printf("%s",StrBuf);
			    AlgState=AlgState|0x0001;
			    return false;
		    }
		}
	}

	AlgState=AlgState&0xFFFE;

	for(int i=0;i<7;i++)
		c_y1p1_ModelPara.flag[i]=0;
	
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
		printf("%s",StrBuf);
		AlgState=AlgState|0x0001;
		return false;
	}

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"totalmodel=",11)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y1p1_ModelPara.totalmodel=strtod(p,NULL)-1;
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y1p1_ModelPara.totalmodel>(MIMOGPC_MODELNUM_MAX-1)||c_y1p1_ModelPara.totalmodel<0){
		c_y1p1_ModelPara.flag[0]=0;
		rewind(fp);
	}
	else{
		c_y1p1_ModelPara.flag[0]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"na=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y1p1_ModelPara.na=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y1p1_ModelPara.na>MIMOGPC_MODEL_NA_MAX||c_y1p1_ModelPara.na<2)
	{
		c_y1p1_ModelPara.flag[1]=0;
		rewind(fp);
	}
	else
	{
		c_y1p1_ModelPara.flag[1]=1;
	}	
	
	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y1p1_ModelPara.nb=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y1p1_ModelPara.nb>MIMOGPC_MODEL_NB_MAX||c_y1p1_ModelPara.nb<1||c_y1p1_ModelPara.nb>c_y1p1_ModelPara.na)
	{
		c_y1p1_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_y1p1_ModelPara.flag[2]=1;
	}
	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb2=",4)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y1p1_ModelPara.nb2=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y1p1_ModelPara.nb2>MIMOGPC_MODEL_NB_MAX||c_y1p1_ModelPara.nb2>c_y1p1_ModelPara.na)
	{
		c_y1p1_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_y1p1_ModelPara.flag[2]=1;
	}



	AlgState=AlgState&0xFFFE;

	if (c_y1p1_ModelPara.flag[0]&&c_y1p1_ModelPara.flag[1]&&c_y1p1_ModelPara.flag[2])//must read "totalmodel" and "na" before
	{
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(str[0]=='A')
			{
				c_y1p1_ModelPara.flag[4]=1;
				for(int i=0;i<=c_y1p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_y1p1_ModelPara.na))
						  {
							c_y1p1_ModelPara.a[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_y1p1_ModelPara.na)//check a[i][0]=1;
							c_y1p1_ModelPara.flag[4]=0;
					}
					else
					{
						c_y1p1_ModelPara.flag[4]=0;
					}

				}
			}
			if(strncmp(str,"B=",2)==0)
			{
				c_y1p1_ModelPara.flag[5]=1;
				for(int i=0;i<=c_y1p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_y1p1_ModelPara.nb))
						  {
							c_y1p1_ModelPara.b[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_y1p1_ModelPara.nb)
							c_y1p1_ModelPara.flag[5]=0;
					}
					else
					{
						c_y1p1_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"B2=",3)==0)
			{
				c_y1p1_ModelPara.flag[5]=1;
				for(int i=0;i<=c_y1p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<c_y1p1_ModelPara.nb2))
						{
							c_y1p1_ModelPara.b2[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j!=c_y1p1_ModelPara.nb2)//check a[i][0]=1;
							c_y1p1_ModelPara.flag[5]=0;
					}
					else
					{
						c_y1p1_ModelPara.flag[5]=0;
					}
				}
			}

			if(strncmp(str,"D=",2)==0)
			{
				c_y1p1_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<=c_y1p1_ModelPara.totalmodel))
						{
							c_y1p1_ModelPara.d[j]=max(m_TunePara.Ts,strtod(p,NULL));
							if (c_y1p1_ModelPara.d[j]<1)
							{
								c_y1p1_ModelPara.d[j]=1;
							//	c_ModelPara.flag[6]=0;
							}
							if(c_y1p1_ModelPara.d_max<c_y1p1_ModelPara.d[j])
								c_y1p1_ModelPara.d_max=c_y1p1_ModelPara.d[j];
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j<c_y1p1_ModelPara.totalmodel)
							c_y1p1_ModelPara.flag[6]=0;
				}
				else
				{
					c_y1p1_ModelPara.flag[6]=0;
				}
			}
			if(strncmp(str,"D2=",3)==0)
			{
				c_y1p1_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
					p=strtok(str,delims);
					j=0;
					while((p!=NULL)&&(j<=c_y1p1_ModelPara.totalmodel))
					{
						c_y1p1_ModelPara.d2[j]=max(m_TunePara.Ts,strtod(p,NULL));
						if (c_y1p1_ModelPara.d2[j]<1)
						{
							c_y1p1_ModelPara.d2[j]=1;
							//	c_ModelPara.flag[6]=0;
						}
						if(c_y1p1_ModelPara.d2_max<c_y1p1_ModelPara.d2[j])
							c_y1p1_ModelPara.d2_max=c_y1p1_ModelPara.d2[j];
						j++;
						p=strtok(NULL,delims); 
					} 
					if (j<c_y1p1_ModelPara.totalmodel)
						c_y1p1_ModelPara.flag[6]=0;
					break;
				}
				else
				{
					c_y1p1_ModelPara.flag[6]=0;
				}
			}


			//fputs(str,stdout);
		}
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Load OK.\n",FilePath);
		printf("%s",StrBuf);
		c_y1p1_ModelPara.T=m_TunePara.Ts;
		c_y1p1_ModelPara.flag[3]=1;//T���������ļ�����
		j=0;//AlgState=AlgState&0xFFFE;AlgState=AlgState|0x0001;		
		for (int i=0;i<7;i++){
			m_y1p1_ModelPara.flag[i]=c_y1p1_ModelPara.flag[i];
			if (1==c_y1p1_ModelPara.flag[i])
			{
				j++;
				AlgState=AlgState&(~(0x0001<<(i+1)));
			}
			else
				AlgState=AlgState|0x0001<<(i+1);
		}

		if (7==j){
			//����ϵͳ��ɢ��na>nb
			m_y1p1_ModelPara.totalmodel=c_y1p1_ModelPara.totalmodel;
			m_y1p1_ModelPara.na=c_y1p1_ModelPara.na;
			m_y1p1_ModelPara.nb=c_y1p1_ModelPara.na;
			m_y1p1_ModelPara.nb2=c_y1p1_ModelPara.na;
			m_y1p1_ModelPara.T=c_y1p1_ModelPara.T*1000;//����Ϊ��λ
			m_y1p1_ModelPara.d_max=1;
			m_y1p1_ModelPara.d2_max=1;

			for(int i=0;i<=m_y1p1_ModelPara.totalmodel;i++)
			{
				if(c_y1p1_ModelPara.nb2!=0)
				{
					myc2d2(c_y1p1_ModelPara.a[i],c_y1p1_ModelPara.b[i],c_y1p1_ModelPara.b2[i],c_y1p1_ModelPara.na,c_y1p1_ModelPara.nb,c_y1p1_ModelPara.nb2,m_y1p1_ModelPara.T/1000,m_y1p1_ModelPara.a[i],m_y1p1_ModelPara.b[i],m_y1p1_ModelPara.b2[i]);
				}
				else
				{
					myc2d(c_y1p1_ModelPara.a[i],c_y1p1_ModelPara.b[i],c_y1p1_ModelPara.na,c_y1p1_ModelPara.nb,m_y1p1_ModelPara.T/1000,m_y1p1_ModelPara.a[i],m_y1p1_ModelPara.b[i]);
				}
				m_y1p1_ModelPara.d[i]=c_y1p1_ModelPara.d[i]/(m_y1p1_ModelPara.T/1000);
				if(c_y1p1_ModelPara.nb2!=0)
				{
					m_y1p1_ModelPara.d2[i]=c_y1p1_ModelPara.d2[i]/(m_y1p1_ModelPara.T/1000);
				}
				if(m_y1p1_ModelPara.d_max<m_y1p1_ModelPara.d[i])
					m_y1p1_ModelPara.d_max=m_y1p1_ModelPara.d[i];
				if(c_y1p1_ModelPara.nb2!=0)
				{
					if(m_y1p1_ModelPara.d2_max<m_y1p1_ModelPara.d2[i])
						m_y1p1_ModelPara.d2_max=m_y1p1_ModelPara.d2[i];
					m_y1p1_ModelPara.d_max = max(m_y1p1_ModelPara.d_max,m_y1p1_ModelPara.d2_max);
				}
				double ta[2]={1,-1};
				conv(m_y1p1_ModelPara.a[i],ta,m_y1p1_ModelPara.na,2,m_y1p1_ModelPara.aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
			}
				return true;
		}
		else{
			//cout<< "input wrong!" << endl;
			return false;		
		}
	}
	else
	{
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
		printf("%s",StrBuf);
		return false;
	}
}

//�������ļ���ȡϵͳģ��
//�ļ���Ϊ��MgpcModel+���+.para;���Ϊ6λ���֣�ǰ��λΪҳ�ţ�������λ��0������λΪID�ţ�������λ��0;
//���磺��5ҳ���㷨IDΪ3�����Ϊ��005003�������ļ���Ϊ��MgpcModel005003.para;
//�ļ�Ŀ¼������Ŀ¼\domain013\data\ctrl\obj\drop001
bool AlgICSMIMOGPC::ReadParay1p2FromFile()//Ҫ�� MR ?
{	
	char filedir[]="\/edpf\/prj\/running\/ctrl\/obj";//���ļ�Ŀ¼
 //   char filedir2[]="E:\\project\\DPU59running\\ctrl\\obj";//������DPU�ļ�Ŀ¼
	char filedir2[]="C:\\Program Files\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
    char filedir3[]="C:\\Program Files (x86)\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
	//char filedir2[]="C:\\Projects\\edpf\\MGPCMatrix";//������DPU�ļ�Ŀ¼
	AlgIdType AlgId,SheetNo,m_RecId,FinalId;
	AlgId=GetId();//��ȡ�㷨Id
	SheetNo=AlgId>>16;
	m_RecId=AlgId&0x000000FF;
	SheetNo=AlgId>>16;

	short int	St;
	char	filename[MAX_PATH]="";
	char	FilePath[MAX_PATH]="";
	char	StrBuf[256];
	size_t	Len=MAX_PATH;
	short int	rc;
	short int	j;

	FILE *fp=NULL;
	char str[MAX_PATH];
	char *p;
	char *delims={" []=\n"};
	// �����ļ���
	sprintf(filename,"Mgpcc1p2Model%03d%03d.para",SheetNo,m_RecId);
	sprintf(FilePath,"%s/%s",filedir,filename);
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(FilePath,"%s\\%s",filedir2,filename);
		if((fp=fopen(FilePath,"r"))==NULL)
		{
			sprintf(FilePath,"%s\\%s",filedir3,filename);
		    if((fp=fopen(FilePath,"r"))==NULL)
		    {
			    sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
			    printf("%s",StrBuf);
			    AlgState=AlgState|0x0001;
			    return false;
		    }
		}
	}

	AlgState=AlgState&0xFFFE;

	for(int i=0;i<7;i++)
		c_y1p2_ModelPara.flag[i]=0;
	
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
		printf("%s",StrBuf);
		AlgState=AlgState|0x0001;
		return false;
	}

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"totalmodel=",11)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y1p2_ModelPara.totalmodel=strtod(p,NULL)-1;
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y1p2_ModelPara.totalmodel>(MIMOGPC_MODELNUM_MAX-1)||c_y1p2_ModelPara.totalmodel<0){
		c_y1p2_ModelPara.flag[0]=0;
		rewind(fp);
	}
	else{
		c_y1p2_ModelPara.flag[0]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"na=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y1p2_ModelPara.na=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y1p2_ModelPara.na>MIMOGPC_MODEL_NA_MAX||c_y1p2_ModelPara.na<2)
	{
		c_y1p2_ModelPara.flag[1]=0;
		rewind(fp);
	}
	else
	{
		c_y1p2_ModelPara.flag[1]=1;
	}	
	
	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y1p2_ModelPara.nb=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y1p2_ModelPara.nb>MIMOGPC_MODEL_NB_MAX||c_y1p2_ModelPara.nb<1||c_y1p2_ModelPara.nb>c_y1p2_ModelPara.na)
	{
		c_y1p2_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_y1p2_ModelPara.flag[2]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb2=",4)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y1p2_ModelPara.nb2=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y1p2_ModelPara.nb2>MIMOGPC_MODEL_NB_MAX||c_y1p2_ModelPara.nb2>c_y1p2_ModelPara.na)
	{
		c_y1p2_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_y1p2_ModelPara.flag[2]=1;
	}


	AlgState=AlgState&0xFFFE;

	if (c_y1p2_ModelPara.flag[0]&&c_y1p2_ModelPara.flag[1]&&c_y1p2_ModelPara.flag[2])//must read "totalmodel" and "na" before
	{
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(str[0]=='A')
			{
				c_y1p2_ModelPara.flag[4]=1;
				for(int i=0;i<=c_y1p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_y1p2_ModelPara.na))
						  {
							c_y1p2_ModelPara.a[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_y1p2_ModelPara.na)//check a[i][0]=1;
							c_y1p2_ModelPara.flag[4]=0;
					}
					else
					{
						c_y1p2_ModelPara.flag[4]=0;
					}

				}
			}
			if(strncmp(str,"B=",2)==0)
			{
				c_y1p2_ModelPara.flag[5]=1;
				for(int i=0;i<=c_y1p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_y1p2_ModelPara.nb))
						  {
							c_y1p2_ModelPara.b[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_y1p2_ModelPara.nb)
							c_y1p2_ModelPara.flag[5]=0;
					}
					else
					{
						c_y1p2_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"B2=",3)==0)
			{
				c_y1p2_ModelPara.flag[5]=1;
				for(int i=0;i<=c_y1p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<c_y1p2_ModelPara.nb2))
						{
							c_y1p2_ModelPara.b2[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j!=c_y1p2_ModelPara.nb2)//check a[i][0]=1;
							c_y1p2_ModelPara.flag[5]=0;
					}
					else
					{
						c_y1p2_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"D=",2)==0)
			{
				c_y1p2_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<=c_y1p2_ModelPara.totalmodel))
						{
							c_y1p2_ModelPara.d[j]=max(m_TunePara.Ts,strtod(p,NULL));
							if (c_y1p2_ModelPara.d[j]<1)
							{
								c_y1p2_ModelPara.d[j]=1;
							//	c_ModelPara.flag[6]=0;
							}
							if(c_y1p2_ModelPara.d_max<c_y1p2_ModelPara.d[j])
								c_y1p2_ModelPara.d_max=c_y1p2_ModelPara.d[j];
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j<c_y1p2_ModelPara.totalmodel)
							c_y1p2_ModelPara.flag[6]=0;
				}
				else
				{
					c_y1p2_ModelPara.flag[6]=0;
				}
			}
			if(strncmp(str,"D2=",3)==0)
			{
				c_y1p2_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
					p=strtok(str,delims);
					j=0;
					while((p!=NULL)&&(j<=c_y1p2_ModelPara.totalmodel))
					{
						c_y1p2_ModelPara.d2[j]=max(m_TunePara.Ts,strtod(p,NULL));
						if (c_y1p2_ModelPara.d2[j]<1)
						{
							c_y1p2_ModelPara.d2[j]=1;
							//	c_ModelPara.flag[6]=0;
						}
						if(c_y1p2_ModelPara.d2_max<c_y1p2_ModelPara.d2[j])
							c_y1p2_ModelPara.d2_max=c_y1p2_ModelPara.d2[j];
						j++;
						p=strtok(NULL,delims); 
					} 
					if (j<c_y1p2_ModelPara.totalmodel)
						c_y1p2_ModelPara.flag[6]=0;
					break;
				}
				else
				{
					c_y1p2_ModelPara.flag[6]=0;
				}
			}
			//fputs(str,stdout);
		}
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Load OK.\n",FilePath);
		printf("%s",StrBuf);
		c_y1p2_ModelPara.T=m_TunePara.Ts;
		c_y1p2_ModelPara.flag[3]=1;//T���������ļ�����
		j=0;//AlgState=AlgState&0xFFFE;AlgState=AlgState|0x0001;		
		for (int i=0;i<7;i++){
			m_y1p2_ModelPara.flag[i]=c_y1p2_ModelPara.flag[i];
			if (1==c_y1p2_ModelPara.flag[i])
			{
				j++;
				AlgState=AlgState&(~(0x0001<<(i+1)));
			}
			else
				AlgState=AlgState|0x0001<<(i+1);
		}

		if (7==j){
			//����ϵͳ��ɢ��na>nb
			m_y1p2_ModelPara.totalmodel=c_y1p2_ModelPara.totalmodel;
			m_y1p2_ModelPara.na=c_y1p2_ModelPara.na;
			m_y1p2_ModelPara.nb=c_y1p2_ModelPara.na;
			m_y1p2_ModelPara.nb2=c_y1p2_ModelPara.na;
			m_y1p2_ModelPara.T=c_y1p2_ModelPara.T*1000;//����Ϊ��λ
			m_y1p2_ModelPara.d_max=1;
			m_y1p2_ModelPara.d2_max=1;

			for(int i=0;i<=m_y1p2_ModelPara.totalmodel;i++)
			{
				if(c_y1p2_ModelPara.nb2!=0)
				{
					myc2d2(c_y1p2_ModelPara.a[i],c_y1p2_ModelPara.b[i],c_y1p2_ModelPara.b2[i],c_y1p2_ModelPara.na,c_y1p2_ModelPara.nb,c_y1p2_ModelPara.nb2,m_y1p2_ModelPara.T/1000,m_y1p2_ModelPara.a[i],m_y1p2_ModelPara.b[i],m_y1p2_ModelPara.b2[i]);
				}
				else
				{
					myc2d(c_y1p2_ModelPara.a[i],c_y1p2_ModelPara.b[i],c_y1p2_ModelPara.na,c_y1p2_ModelPara.nb,m_y1p2_ModelPara.T/1000,m_y1p2_ModelPara.a[i],m_y1p2_ModelPara.b[i]);
				}
				m_y1p2_ModelPara.d[i]=c_y1p2_ModelPara.d[i]/(m_y1p2_ModelPara.T/1000);
				if(c_y1p2_ModelPara.nb2!=0)
				{
					m_y1p2_ModelPara.d2[i]=c_y1p2_ModelPara.d2[i]/(m_y1p2_ModelPara.T/1000);
				}
				if(m_y1p2_ModelPara.d_max<m_y1p2_ModelPara.d[i])
					m_y1p2_ModelPara.d_max=m_y1p2_ModelPara.d[i];
				if(c_y1p2_ModelPara.nb2!=0)
				{
					if(m_y1p2_ModelPara.d2_max<m_y1p2_ModelPara.d2[i])
						m_y1p2_ModelPara.d2_max=m_y1p2_ModelPara.d2[i];
					m_y1p2_ModelPara.d_max = max(m_y1p2_ModelPara.d_max,m_y1p2_ModelPara.d2_max);
				}
				double ta[2]={1,-1};
				conv(m_y1p2_ModelPara.a[i],ta,m_y1p2_ModelPara.na,2,m_y1p2_ModelPara.aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
			}
				return true;
		}
		else{
			//cout<< "input wrong!" << endl;
			return false;		
		}
	}
	else
	{
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
		printf("%s",StrBuf);
		return false;
	}
}

//�������ļ���ȡϵͳģ��
//�ļ���Ϊ��MgpcModel+���+.para;���Ϊ6λ���֣�ǰ��λΪҳ�ţ�������λ��0������λΪID�ţ�������λ��0;
//���磺��5ҳ���㷨IDΪ3�����Ϊ��005003�������ļ���Ϊ��MgpcModel005003.para;
//�ļ�Ŀ¼������Ŀ¼\domain013\data\ctrl\obj\drop001
bool AlgICSMIMOGPC::ReadParay2p1FromFile()//Ҫ�� MR ?
{	
	char filedir[]="\/edpf\/prj\/running\/ctrl\/obj";//���ļ�Ŀ¼
 //   char filedir2[]="E:\\project\\DPU59running\\ctrl\\obj";//������DPU�ļ�Ŀ¼
	char filedir2[]="C:\\Program Files\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
    char filedir3[]="C:\\Program Files (x86)\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
	//char filedir2[]="C:\\Projects\\edpf\\MGPCMatrix";//������DPU�ļ�Ŀ¼
	AlgIdType AlgId,SheetNo,m_RecId,FinalId;
	AlgId=GetId();//��ȡ�㷨Id
	SheetNo=AlgId>>16;
	m_RecId=AlgId&0x000000FF;
	SheetNo=AlgId>>16;

	short int	St;
	char	filename[MAX_PATH]="";
	char	FilePath[MAX_PATH]="";
	char	StrBuf[256];
	size_t	Len=MAX_PATH;
	short int	rc;
	short int	j;

	FILE *fp=NULL;
	char str[MAX_PATH];
	char *p;
	char *delims={" []=\n"};
	// �����ļ���
	sprintf(filename,"Mgpcc2p1Model%03d%03d.para",SheetNo,m_RecId);
	sprintf(FilePath,"%s/%s",filedir,filename);
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(FilePath,"%s\\%s",filedir2,filename);
		if((fp=fopen(FilePath,"r"))==NULL)
		{
			sprintf(FilePath,"%s\\%s",filedir3,filename);
		    if((fp=fopen(FilePath,"r"))==NULL)
		    {
			    sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
			    printf("%s",StrBuf);
			    AlgState=AlgState|0x0001;
			    return false;
		    }
		}
	}

	AlgState=AlgState&0xFFFE;

	for(int i=0;i<7;i++)
		c_y2p1_ModelPara.flag[i]=0;
	
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
		printf("%s",StrBuf);
		AlgState=AlgState|0x0001;
		return false;
	}

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"totalmodel=",11)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y2p1_ModelPara.totalmodel=strtod(p,NULL)-1;
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y2p1_ModelPara.totalmodel>(MIMOGPC_MODELNUM_MAX-1)||c_y2p1_ModelPara.totalmodel<0){
		c_y2p1_ModelPara.flag[0]=0;
		rewind(fp);
	}
	else{
		c_y2p1_ModelPara.flag[0]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"na=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y2p1_ModelPara.na=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y2p1_ModelPara.na>MIMOGPC_MODEL_NA_MAX||c_y2p1_ModelPara.na<2)
	{
		c_y2p1_ModelPara.flag[1]=0;
		rewind(fp);
	}
	else
	{
		c_y2p1_ModelPara.flag[1]=1;
	}	
	
	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y2p1_ModelPara.nb=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y2p1_ModelPara.nb>MIMOGPC_MODEL_NB_MAX||c_y2p1_ModelPara.nb<1||c_y2p1_ModelPara.nb>c_y2p1_ModelPara.na)
	{
		c_y2p1_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_y2p1_ModelPara.flag[2]=1;
	}
	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb2=",4)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y2p1_ModelPara.nb2=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y2p1_ModelPara.nb2>MIMOGPC_MODEL_NB_MAX||c_y2p1_ModelPara.nb2>c_y2p1_ModelPara.na)
	{
		c_y2p1_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_y2p1_ModelPara.flag[2]=1;
	}

	AlgState=AlgState&0xFFFE;

	if (c_y2p1_ModelPara.flag[0]&&c_y2p1_ModelPara.flag[1]&&c_y2p1_ModelPara.flag[2])//must read "totalmodel" and "na" before
	{
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(str[0]=='A')
			{
				c_y2p1_ModelPara.flag[4]=1;
				for(int i=0;i<=c_y2p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_y2p1_ModelPara.na))
						  {
							c_y2p1_ModelPara.a[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_y2p1_ModelPara.na)//check a[i][0]=1;
							c_y2p1_ModelPara.flag[4]=0;
					}
					else
					{
						c_y2p1_ModelPara.flag[4]=0;
					}

				}
			}
			if(strncmp(str,"B=",2)==0)
			{
				c_y2p1_ModelPara.flag[5]=1;
				for(int i=0;i<=c_y2p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_y2p1_ModelPara.nb))
						  {
							c_y2p1_ModelPara.b[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_y2p1_ModelPara.nb)
							c_y2p1_ModelPara.flag[5]=0;
					}
					else
					{
						c_y2p1_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"B2=",3)==0)
			{
				c_y2p1_ModelPara.flag[5]=1;
				for(int i=0;i<=c_y2p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<c_y2p1_ModelPara.nb2))
						{
							c_y2p1_ModelPara.b2[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j!=c_y2p1_ModelPara.nb2)//check a[i][0]=1;
							c_y2p1_ModelPara.flag[5]=0;
					}
					else
					{
						c_y2p1_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"D=",2)==0)
			{
				c_y2p1_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<=c_y2p1_ModelPara.totalmodel))
						{
							c_y2p1_ModelPara.d[j]=max(m_TunePara.Ts,strtod(p,NULL));
							if (c_y2p1_ModelPara.d[j]<1)
							{
								c_y2p1_ModelPara.d[j]=1;
							//	c_ModelPara.flag[6]=0;
							}
							if(c_y2p1_ModelPara.d_max<c_y2p1_ModelPara.d[j])
								c_y2p1_ModelPara.d_max=c_y2p1_ModelPara.d[j];
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j<c_y2p1_ModelPara.totalmodel)
							c_y2p1_ModelPara.flag[6]=0;
				}
				else
				{
					c_y2p1_ModelPara.flag[6]=0;
				}
			}
			if(strncmp(str,"D2=",3)==0)
			{
				c_y2p1_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
					p=strtok(str,delims);
					j=0;
					while((p!=NULL)&&(j<=c_y2p1_ModelPara.totalmodel))
					{
						c_y2p1_ModelPara.d2[j]=max(m_TunePara.Ts,strtod(p,NULL));
						if (c_y2p1_ModelPara.d2[j]<1)
						{
							c_y2p1_ModelPara.d2[j]=1;
							//	c_ModelPara.flag[6]=0;
						}
						if(c_y2p1_ModelPara.d2_max<c_y2p1_ModelPara.d2[j])
							c_y2p1_ModelPara.d2_max=c_y2p1_ModelPara.d2[j];
						j++;
						p=strtok(NULL,delims); 
					} 
					if (j<c_y2p1_ModelPara.totalmodel)
						c_y2p1_ModelPara.flag[6]=0;
					break;
				}
				else
				{
					c_y2p1_ModelPara.flag[6]=0;
				}
			}
			//fputs(str,stdout);
		}
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Load OK.\n",FilePath);
		printf("%s",StrBuf);
		c_y2p1_ModelPara.T=m_TunePara.Ts;
		c_y2p1_ModelPara.flag[3]=1;//T���������ļ�����
		j=0;//AlgState=AlgState&0xFFFE;AlgState=AlgState|0x0001;		
		for (int i=0;i<7;i++){
			m_y2p1_ModelPara.flag[i]=c_y2p1_ModelPara.flag[i];
			if (1==c_y2p1_ModelPara.flag[i])
			{
				j++;
				AlgState=AlgState&(~(0x0001<<(i+1)));
			}
			else
				AlgState=AlgState|0x0001<<(i+1);
		}

		if (7==j){
			//����ϵͳ��ɢ��na>nb
			m_y2p1_ModelPara.totalmodel=c_y2p1_ModelPara.totalmodel;
			m_y2p1_ModelPara.na=c_y2p1_ModelPara.na;
			m_y2p1_ModelPara.nb=c_y2p1_ModelPara.na;
			m_y2p1_ModelPara.nb2=c_y2p1_ModelPara.na;
			m_y2p1_ModelPara.T=c_y2p1_ModelPara.T*1000;//����Ϊ��λ
			m_y2p1_ModelPara.d_max=1;
			m_y2p1_ModelPara.d2_max=1;

			for(int i=0;i<=m_y2p1_ModelPara.totalmodel;i++)
			{
				if(c_y2p1_ModelPara.nb2!=0)
				{
					myc2d2(c_y2p1_ModelPara.a[i],c_y2p1_ModelPara.b[i],c_y2p1_ModelPara.b2[i],c_y2p1_ModelPara.na,c_y2p1_ModelPara.nb,c_y2p1_ModelPara.nb2,m_y2p1_ModelPara.T/1000,m_y2p1_ModelPara.a[i],m_y2p1_ModelPara.b[i],m_y2p1_ModelPara.b2[i]);
				}
				else
				{
					myc2d(c_y2p1_ModelPara.a[i],c_y2p1_ModelPara.b[i],c_y2p1_ModelPara.na,c_y2p1_ModelPara.nb,m_y2p1_ModelPara.T/1000,m_y2p1_ModelPara.a[i],m_y2p1_ModelPara.b[i]);
				}
				m_y2p1_ModelPara.d[i]=c_y2p1_ModelPara.d[i]/(m_y2p1_ModelPara.T/1000);
				if(c_y2p1_ModelPara.nb2!=0)
				{
					m_y2p1_ModelPara.d2[i]=c_y2p1_ModelPara.d2[i]/(m_y2p1_ModelPara.T/1000);
				}
				if(m_y2p1_ModelPara.d_max<m_y2p1_ModelPara.d[i])
					m_y2p1_ModelPara.d_max=m_y2p1_ModelPara.d[i];
				if(c_y2p1_ModelPara.nb2!=0)
				{
					if(m_y2p1_ModelPara.d2_max<m_y2p1_ModelPara.d2[i])
						m_y2p1_ModelPara.d2_max=m_y2p1_ModelPara.d2[i];
					m_y2p1_ModelPara.d_max = max(m_y2p1_ModelPara.d_max,m_y2p1_ModelPara.d2_max);
				}
				double ta[2]={1,-1};
				conv(m_y2p1_ModelPara.a[i],ta,m_y2p1_ModelPara.na,2,m_y2p1_ModelPara.aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
			}
				return true;
		}
		else{
			//cout<< "input wrong!" << endl;
			return false;		
		}
	}
	else
	{
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
		printf("%s",StrBuf);
		return false;
	}
}

//�������ļ���ȡϵͳģ��
//�ļ���Ϊ��MgpcModel+���+.para;���Ϊ6λ���֣�ǰ��λΪҳ�ţ�������λ��0������λΪID�ţ�������λ��0;
//���磺��5ҳ���㷨IDΪ3�����Ϊ��005003�������ļ���Ϊ��MgpcModel005003.para;
//�ļ�Ŀ¼������Ŀ¼\domain013\data\ctrl\obj\drop001
bool AlgICSMIMOGPC::ReadParay2p2FromFile()//Ҫ�� MR ?
{	
	char filedir[]="\/edpf\/prj\/running\/ctrl\/obj";//���ļ�Ŀ¼
 //   char filedir2[]="E:\\project\\DPU59running\\ctrl\\obj";//������DPU�ļ�Ŀ¼
	char filedir2[]="C:\\Program Files\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
    char filedir3[]="C:\\Program Files (x86)\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
	//char filedir2[]="C:\\Projects\\edpf\\MGPCMatrix";//������DPU�ļ�Ŀ¼
	AlgIdType AlgId,SheetNo,m_RecId,FinalId;
	AlgId=GetId();//��ȡ�㷨Id
	SheetNo=AlgId>>16;
	m_RecId=AlgId&0x000000FF;
	SheetNo=AlgId>>16;

	short int	St;
	char	filename[MAX_PATH]="";
	char	FilePath[MAX_PATH]="";
	char	StrBuf[256];
	size_t	Len=MAX_PATH;
	short int	rc;
	short int	j;

	FILE *fp=NULL;
	char str[MAX_PATH];
	char *p;
	char *delims={" []=\n"};
	// �����ļ���
	sprintf(filename,"Mgpcc2p2Model%03d%03d.para",SheetNo,m_RecId);
	sprintf(FilePath,"%s/%s",filedir,filename);
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(FilePath,"%s\\%s",filedir2,filename);
		if((fp=fopen(FilePath,"r"))==NULL)
		{
			sprintf(FilePath,"%s\\%s",filedir3,filename);
		    if((fp=fopen(FilePath,"r"))==NULL)
		    {
			    sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
			    printf("%s",StrBuf);
			    AlgState=AlgState|0x0001;
			    return false;
		    }
		}
	}

	AlgState=AlgState&0xFFFE;

	for(int i=0;i<7;i++)
		c_y2p2_ModelPara.flag[i]=0;
	
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
		printf("%s",StrBuf);
		AlgState=AlgState|0x0001;
		return false;
	}

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"totalmodel=",11)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y2p2_ModelPara.totalmodel=strtod(p,NULL)-1;
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y2p2_ModelPara.totalmodel>(MIMOGPC_MODELNUM_MAX-1)||c_y2p2_ModelPara.totalmodel<0){
		c_y2p2_ModelPara.flag[0]=0;
		rewind(fp);
	}
	else{
		c_y2p2_ModelPara.flag[0]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"na=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y2p2_ModelPara.na=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y2p2_ModelPara.na>MIMOGPC_MODEL_NA_MAX||c_y2p2_ModelPara.na<2)
	{
		c_y2p2_ModelPara.flag[1]=0;
		rewind(fp);
	}
	else
	{
		c_y2p2_ModelPara.flag[1]=1;
	}	
	
	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y2p2_ModelPara.nb=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y2p2_ModelPara.nb>MIMOGPC_MODEL_NB_MAX||c_y2p2_ModelPara.nb<1||c_y2p2_ModelPara.nb>c_y2p2_ModelPara.na)
	{
		c_y2p2_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_y2p2_ModelPara.flag[2]=1;
	}

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb2=",4)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_y2p2_ModelPara.nb2=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_y2p2_ModelPara.nb2>MIMOGPC_MODEL_NB_MAX||c_y2p2_ModelPara.nb2>c_y2p2_ModelPara.na)
	{
		c_y2p2_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_y2p2_ModelPara.flag[2]=1;
	}


	AlgState=AlgState&0xFFFE;

	if (c_y2p2_ModelPara.flag[0]&&c_y2p2_ModelPara.flag[1]&&c_y2p2_ModelPara.flag[2])//must read "totalmodel" and "na" before
	{
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(str[0]=='A')
			{
				c_y2p2_ModelPara.flag[4]=1;
				for(int i=0;i<=c_y2p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_y2p2_ModelPara.na))
						  {
							c_y2p2_ModelPara.a[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_y2p2_ModelPara.na)//check a[i][0]=1;
							c_y2p2_ModelPara.flag[4]=0;
					}
					else
					{
						c_y2p2_ModelPara.flag[4]=0;
					}

				}
			}
			if(strncmp(str,"B=",2)==0)
			{
				c_y2p2_ModelPara.flag[5]=1;
				for(int i=0;i<=c_y2p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_y2p2_ModelPara.nb))
						  {
							c_y2p2_ModelPara.b[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_y2p2_ModelPara.nb)
							c_y2p2_ModelPara.flag[5]=0;
					}
					else
					{
						c_y2p2_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"B2=",3)==0)
			{
				c_y2p2_ModelPara.flag[5]=1;
				for(int i=0;i<=c_y2p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<c_y2p2_ModelPara.nb2))
						{
							c_y2p2_ModelPara.b2[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j!=c_y2p2_ModelPara.nb2)//check a[i][0]=1;
							c_y2p2_ModelPara.flag[5]=0;
					}
					else
					{
						c_y2p2_ModelPara.flag[5]=0;
					}
				}
			}

			if(strncmp(str,"D=",2)==0)
			{
				c_y2p2_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<=c_y2p2_ModelPara.totalmodel))
						{
							c_y2p2_ModelPara.d[j]=max(m_TunePara.Ts,strtod(p,NULL));
							if (c_y2p2_ModelPara.d[j]<1)
							{
								c_y2p2_ModelPara.d[j]=1;
							//	c_ModelPara.flag[6]=0;
							}
							if(c_y2p2_ModelPara.d_max<c_y2p2_ModelPara.d[j])
								c_y2p2_ModelPara.d_max=c_y2p2_ModelPara.d[j];
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j<c_y2p2_ModelPara.totalmodel)
							c_y2p2_ModelPara.flag[6]=0;
				}
				else
				{
					c_y2p2_ModelPara.flag[6]=0;
				}
			}
			if(strncmp(str,"D2=",3)==0)
			{
				c_y2p2_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
					p=strtok(str,delims);
					j=0;
					while((p!=NULL)&&(j<=c_y2p2_ModelPara.totalmodel))
					{
						c_y2p2_ModelPara.d2[j]=max(m_TunePara.Ts,strtod(p,NULL));
						if (c_y2p2_ModelPara.d2[j]<1)
						{
							c_y2p2_ModelPara.d2[j]=1;
							//	c_ModelPara.flag[6]=0;
						}
						if(c_y2p2_ModelPara.d2_max<c_y2p2_ModelPara.d2[j])
							c_y2p2_ModelPara.d2_max=c_y2p2_ModelPara.d2[j];
						j++;
						p=strtok(NULL,delims); 
					} 
					if (j<c_y2p2_ModelPara.totalmodel)
						c_y2p2_ModelPara.flag[6]=0;
					break;
				}
				else
				{
					c_y2p2_ModelPara.flag[6]=0;
				}
			}
			//fputs(str,stdout);
		}
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Load OK.\n",FilePath);
		printf("%s",StrBuf);
		c_y2p2_ModelPara.T=m_TunePara.Ts;
		c_y2p2_ModelPara.flag[3]=1;//T���������ļ�����
		j=0;//AlgState=AlgState&0xFFFE;AlgState=AlgState|0x0001;		
		for (int i=0;i<7;i++){
			m_y2p2_ModelPara.flag[i]=c_y2p2_ModelPara.flag[i];
			if (1==c_y2p2_ModelPara.flag[i])
			{
				j++;
				AlgState=AlgState&(~(0x0001<<(i+1)));
			}
			else
				AlgState=AlgState|0x0001<<(i+1);
		}

		if (7==j){
			//����ϵͳ��ɢ��na>nb
			m_y2p2_ModelPara.totalmodel=c_y2p2_ModelPara.totalmodel;
			m_y2p2_ModelPara.na=c_y2p2_ModelPara.na;
			m_y2p2_ModelPara.nb=c_y2p2_ModelPara.na;
			m_y2p2_ModelPara.nb2=c_y2p2_ModelPara.na;
			m_y2p2_ModelPara.T=c_y2p2_ModelPara.T*1000;//����Ϊ��λ
			m_y2p2_ModelPara.d_max=1;
			m_y2p2_ModelPara.d2_max=1;

			for(int i=0;i<=m_y2p2_ModelPara.totalmodel;i++)
			{
				if(c_y2p2_ModelPara.nb2!=0)
				{
					myc2d2(c_y2p2_ModelPara.a[i],c_y2p2_ModelPara.b[i],c_y2p2_ModelPara.b2[i],c_y2p2_ModelPara.na,c_y2p2_ModelPara.nb,c_y2p2_ModelPara.nb2,m_y2p2_ModelPara.T/1000,m_y2p2_ModelPara.a[i],m_y2p2_ModelPara.b[i],m_y2p2_ModelPara.b2[i]);
				}
				else
				{
					myc2d(c_y2p2_ModelPara.a[i],c_y2p2_ModelPara.b[i],c_y2p2_ModelPara.na,c_y2p2_ModelPara.nb,m_y2p2_ModelPara.T/1000,m_y2p2_ModelPara.a[i],m_y2p2_ModelPara.b[i]);
				}
				m_y2p2_ModelPara.d[i]=c_y2p2_ModelPara.d[i]/(m_y2p2_ModelPara.T/1000);
				if(c_y2p2_ModelPara.nb2!=0)
				{
					m_y2p2_ModelPara.d2[i]=c_y2p2_ModelPara.d2[i]/(m_y2p2_ModelPara.T/1000);
				}
				if(m_y2p2_ModelPara.d_max<m_y2p2_ModelPara.d[i])
					m_y2p2_ModelPara.d_max=m_y2p2_ModelPara.d[i];
				if(c_y2p2_ModelPara.nb2!=0)
				{
					if(m_y2p2_ModelPara.d2_max<m_y2p2_ModelPara.d2[i])
						m_y2p2_ModelPara.d2_max=m_y2p2_ModelPara.d2[i];
					m_y2p2_ModelPara.d_max = max(m_y2p2_ModelPara.d_max,m_y2p2_ModelPara.d2_max);
				}
				double ta[2]={1,-1};
				conv(m_y2p2_ModelPara.a[i],ta,m_y2p2_ModelPara.na,2,m_y2p2_ModelPara.aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
			}
				return true;
		}
		else{
			//cout<< "input wrong!" << endl;
			return false;		
		}
	}
	else
	{
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
		printf("%s",StrBuf);
		return false;
	}
}

/*==============================================����һ========================================================*/
//�������ļ���ȡϵͳģ��
//�ļ���Ϊ��MgpcModel+���+.para;���Ϊ6λ���֣�ǰ��λΪҳ�ţ�������λ��0������λΪID�ţ�������λ��0;
//���磺��5ҳ���㷨IDΪ3�����Ϊ��005003�������ļ���Ϊ��MgpcModel005003.para;
//�ļ�Ŀ¼������Ŀ¼\domain013\data\ctrl\obj\drop001
bool AlgICSMIMOGPC::ReadParav1p1FromFile()//Ҫ�� MR ?
{	
	char filedir[]="\/edpf\/prj\/running\/ctrl\/obj";//���ļ�Ŀ¼
    //char filedir2[]="E:\\project\\DPU59running\\ctrl\\obj";//������DPU�ļ�Ŀ¼
	char filedir2[]="C:\\Program Files\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
    char filedir3[]="C:\\Program Files (x86)\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
	//char filedir2[]="C:\\Projects\\edpf\\MGPCMatrix";//������DPU�ļ�Ŀ¼
	AlgIdType AlgId,SheetNo,m_RecId,FinalId;
	AlgId=GetId();//��ȡ�㷨Id
	SheetNo=AlgId>>16;
	m_RecId=AlgId&0x000000FF;
	SheetNo=AlgId>>16;

	short int	St;
	char	filename[MAX_PATH]="";
	char	FilePath[MAX_PATH]="";
	char	StrBuf[256];
	size_t	Len=MAX_PATH;
	short int	rc;
	short int	j;

	FILE *fp=NULL;
	char str[MAX_PATH];
	char *p;
	char *delims={" []=\n"};
	// �����ļ���
	sprintf(filename,"Mgpcv1p1Model%03d%03d.para",SheetNo,m_RecId);
	sprintf(FilePath,"%s/%s",filedir,filename);
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(FilePath,"%s\\%s",filedir2,filename);
		if((fp=fopen(FilePath,"r"))==NULL)
		{
			sprintf(FilePath,"%s\\%s",filedir3,filename);
		    if((fp=fopen(FilePath,"r"))==NULL)
		    {
			    sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
			    printf("%s",StrBuf);
			    AlgState=AlgState|0x0001;
			    return false;
		    }
		}
	}

	AlgState=AlgState&0xFFFE;

	for(int i=0;i<7;i++)
		c_v1p1_ModelPara.flag[i]=0;
	
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
		printf("%s",StrBuf);
		AlgState=AlgState|0x0001;
		return false;
	}

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"totalmodel=",11)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v1p1_ModelPara.totalmodel=strtod(p,NULL)-1;
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v1p1_ModelPara.totalmodel>(M1MIMOGPC_MODELNUM_MAX-1)||c_v1p1_ModelPara.totalmodel<0){
		c_v1p1_ModelPara.flag[0]=0;
		rewind(fp);
	}
	else{
		c_v1p1_ModelPara.flag[0]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"na=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v1p1_ModelPara.na=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v1p1_ModelPara.na>M1MIMOGPC_MODEL_NA_MAX||c_v1p1_ModelPara.na<2)
	{
		c_v1p1_ModelPara.flag[1]=0;
		rewind(fp);
	}
	else
	{
		c_v1p1_ModelPara.flag[1]=1;
	}	
	
	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v1p1_ModelPara.nb=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v1p1_ModelPara.nb>M1MIMOGPC_MODEL_NB_MAX||c_v1p1_ModelPara.nb<1||c_v1p1_ModelPara.nb>c_v1p1_ModelPara.na)
	{
		c_v1p1_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_v1p1_ModelPara.flag[2]=1;
	}

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb2=",4)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v1p1_ModelPara.nb2=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v1p1_ModelPara.nb2>M1MIMOGPC_MODEL_NB_MAX||c_v1p1_ModelPara.nb2>c_v1p1_ModelPara.na)
	{
		c_v1p1_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_v1p1_ModelPara.flag[2]=1;
	}

	AlgState=AlgState&0xFFFE;

	if (c_v1p1_ModelPara.flag[0]&&c_v1p1_ModelPara.flag[1]&&c_v1p1_ModelPara.flag[2])//must read "totalmodel" and "na" before
	{
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(str[0]=='A')
			{
				c_v1p1_ModelPara.flag[4]=1;
				for(int i=0;i<=c_v1p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_v1p1_ModelPara.na))
						  {
							c_v1p1_ModelPara.a[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_v1p1_ModelPara.na)//check a[i][0]=1;
							c_v1p1_ModelPara.flag[4]=0;
					}
					else
					{
						c_v1p1_ModelPara.flag[4]=0;
					}

				}
			}
			if(strncmp(str,"B=",2)==0)
			{
				c_v1p1_ModelPara.flag[5]=1;
				for(int i=0;i<=c_v1p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_v1p1_ModelPara.nb))
						  {
							c_v1p1_ModelPara.b[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_v1p1_ModelPara.nb)
							c_v1p1_ModelPara.flag[5]=0;
					}
					else
					{
						c_v1p1_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"B2=",3)==0)
			{
				c_v1p1_ModelPara.flag[5]=1;
				for(int i=0;i<=c_v1p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<c_v1p1_ModelPara.nb2))
						{
							c_v1p1_ModelPara.b2[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j!=c_v1p1_ModelPara.nb2)//check a[i][0]=1;
							c_v1p1_ModelPara.flag[5]=0;
					}
					else
					{
						c_v1p1_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"D=",2)==0)
			{
				c_v1p1_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<=c_v1p1_ModelPara.totalmodel))
						{
							c_v1p1_ModelPara.d[j]=max(m_TunePara.Ts,strtod(p,NULL));
							if (c_v1p1_ModelPara.d[j]<1)
							{
								c_v1p1_ModelPara.d[j]=1;
								//c1_ModelPara.flag[6]=0;
							}
							if(c_v1p1_ModelPara.d_max<c_v1p1_ModelPara.d[j])
								c_v1p1_ModelPara.d_max=c_v1p1_ModelPara.d[j];
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j<c_v1p1_ModelPara.totalmodel)
							c_v1p1_ModelPara.flag[6]=0;
				}
				else
				{
					c_v1p1_ModelPara.flag[6]=0;
				}
			}
			if(strncmp(str,"D2=",3)==0)
			{
				c_v1p1_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
					p=strtok(str,delims);
					j=0;
					while((p!=NULL)&&(j<=c_v1p1_ModelPara.totalmodel))
					{
						c_v1p1_ModelPara.d2[j]=max(m_TunePara.Ts,strtod(p,NULL));
						if (c_v1p1_ModelPara.d2[j]<1)
						{
							c_v1p1_ModelPara.d2[j]=1;
							//	c_ModelPara.flag[6]=0;
						}
						if(c_v1p1_ModelPara.d2_max<c_v1p1_ModelPara.d2[j])
							c_v1p1_ModelPara.d2_max=c_v1p1_ModelPara.d2[j];
						j++;
						p=strtok(NULL,delims); 
					} 
					if (j<c_v1p1_ModelPara.totalmodel)
						c_v1p1_ModelPara.flag[6]=0;
					break;
				}
				else
				{
					c_v1p1_ModelPara.flag[6]=0;
				}
			}
			//fputs(str,stdout);
		}
//��NLTD����
		rc=1;
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(strncmp(str,"NLTDFLAG=",9)==0)
			{
				p=strtok(str,delims);
				if(p!=NULL)
				{
					p=strtok(NULL,delims);
					c_v1p1_ModelPara.NLTDflag=strtod(p,NULL);
					rc=0;
					break;
				}
			}
		}
		if (rc||c_v1p1_ModelPara.NLTDflag>1||c_v1p1_ModelPara.NLTDflag<0)
		{
			c_v1p1_ModelPara.flag[0]=0;
			rewind(fp);
		}
		else
		{
			c_v1p1_ModelPara.flag[0]=1;
		}
		if(c_v1p1_ModelPara.NLTDflag==1)
		{
			rc=1;
			while(fgets(str,sizeof(str),fp)!=NULL)
			{
				if(strncmp(str,"NLTDna=",7)==0)
				{
					p=strtok(str,delims);
					if(p!=NULL)
					{
						p=strtok(NULL,delims);
						c_v1p1_ModelPara.NLTDna=strtod(p,NULL);
						rc=0;
						break;
					}
				}
			}
			if (rc||c_v1p1_ModelPara.NLTDna>M1MIMOGPC_MODEL_NA_MAX||c_v1p1_ModelPara.NLTDna<2)
			{
				c_v1p1_ModelPara.flag[1]=0;
				rewind(fp);
			}
			else
			{
				c_v1p1_ModelPara.flag[1]=1;
			}	
	
			rc=1;
			while(fgets(str,sizeof(str),fp)!=NULL)
			{
				if(strncmp(str,"NLTDnb=",7)==0)
				{
					p=strtok(str,delims);
					if(p!=NULL)
					{
						p=strtok(NULL,delims);
						c_v1p1_ModelPara.NLTDnb=strtod(p,NULL);
						rc=0;
						break;
					}
				}
			}
			if (rc||c_v1p1_ModelPara.NLTDnb>M1MIMOGPC_MODEL_NB_MAX||c_v1p1_ModelPara.NLTDnb<1||c_v1p1_ModelPara.NLTDnb>c_v1p1_ModelPara.NLTDna)
			{
				c_v1p1_ModelPara.flag[2]=0;
				rewind(fp);
			}
			else
			{
				c_v1p1_ModelPara.flag[2]=1;
			}	
			if (c_v1p1_ModelPara.flag[0]&&c_v1p1_ModelPara.flag[1]&&c_v1p1_ModelPara.flag[2])
			{
				while(fgets(str,sizeof(str),fp)!=NULL)
				{
					if(strncmp(str,"NLTDA=",6)==0)
					{
						c_v1p1_ModelPara.flag[4]=1;
						for(int i=0;i<=c_v1p1_ModelPara.totalmodel;i++)
						{
							if(fgets(str,sizeof(str),fp)!=NULL)
							{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<c_v1p1_ModelPara.NLTDna))
								{
									c_v1p1_ModelPara.NLTD_a[i][j]=strtod(p,NULL);
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j!=c_v1p1_ModelPara.NLTDna)//check a[i][0]=1;
									c_v1p1_ModelPara.flag[4]=0;
							}
							else
							{
								c_v1p1_ModelPara.flag[4]=0;
							}

						}
					}
					if(strncmp(str,"NLTDB=",6)==0)
					{
						c_v1p1_ModelPara.flag[5]=1;
						for(int i=0;i<=c_v1p1_ModelPara.totalmodel;i++)
						{
							if(fgets(str,sizeof(str),fp)!=NULL)
							{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<c_v1p1_ModelPara.NLTDnb))
								{
									c_v1p1_ModelPara.NLTD_b[i][j]=strtod(p,NULL);
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j!=c_v1p1_ModelPara.NLTDnb)
									c_v1p1_ModelPara.flag[5]=0;
							}
							else
							{
								c_v1p1_ModelPara.flag[5]=0;
							}
						}
					}
					if(strncmp(str,"NLTDD=",6)==0)
					{
						c_v1p1_ModelPara.flag[6]=1;
						if(fgets(str,sizeof(str),fp)!=NULL)
						{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<=c_v1p1_ModelPara.totalmodel))
								{
									c_v1p1_ModelPara.NLTD_d[j]=max(m_TunePara.Ts,strtod(p,NULL));
									if (c_v1p1_ModelPara.NLTD_d[j]<1)
									{
										c_v1p1_ModelPara.NLTD_d[j]=1;
										//c1_ModelPara.flag[6]=0;
									}
									if(c_v1p1_ModelPara.NLTD_d_max<c_v1p1_ModelPara.NLTD_d[j])
										c_v1p1_ModelPara.NLTD_d_max=c_v1p1_ModelPara.NLTD_d[j];
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j<c_v1p1_ModelPara.totalmodel)
									c_v1p1_ModelPara.flag[6]=0;
						}
						else
						{
							c_v1p1_ModelPara.flag[6]=0;
						}
					}
				}
			}
			else
			{
				fclose(fp);
				sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
				printf("%s",StrBuf);
				return false;
			}
		
		}
		
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Load OK.\n",FilePath);
		printf("%s",StrBuf);
		c_v1p1_ModelPara.T=m_TunePara.Ts;
		c_v1p1_ModelPara.flag[3]=1;//T���������ļ�����
		j=0;//AlgState=AlgState&0xFFFE;AlgState=AlgState|0x0001;		
		for (int i=0;i<7;i++){
			m_v1p1_ModelPara.flag[i]=c_v1p1_ModelPara.flag[i];
			if (1==c_v1p1_ModelPara.flag[i])
			{
				j++;
				AlgState=AlgState&(~(0x0001<<(i+1)));
			}
			else
				AlgState=AlgState|0x0001<<(i+1);
		}

		if (7==j)
		{
			//����ϵͳ��ɢ��na>nb
			m_v1p1_ModelPara.totalmodel=c_v1p1_ModelPara.totalmodel;
			m_v1p1_ModelPara.na=c_v1p1_ModelPara.na;
			m_v1p1_ModelPara.nb=c_v1p1_ModelPara.na;
			m_v1p1_ModelPara.nb2=c_v1p1_ModelPara.na;
			m_v1p1_ModelPara.T=c_v1p1_ModelPara.T*1000;//����Ϊ��λ
			m_v1p1_ModelPara.d_max=1;
			m_v1p1_ModelPara.d2_max=1;

			m_v1p1_ModelPara.NLTDna=c_v1p1_ModelPara.NLTDna;
			m_v1p1_ModelPara.NLTDnb=c_v1p1_ModelPara.NLTDna;
			m_v1p1_ModelPara.NLTD_d_max=1;

			for(int i=0;i<=m_v1p1_ModelPara.totalmodel;i++)
			{
				if(c_v1p1_ModelPara.nb2!=0)
				{
					myc2d2(c_v1p1_ModelPara.a[i],c_v1p1_ModelPara.b[i],c_v1p1_ModelPara.b2[i],c_v1p1_ModelPara.na,c_v1p1_ModelPara.nb,c_v1p1_ModelPara.nb2,m_v1p1_ModelPara.T/1000,m_v1p1_ModelPara.a[i],m_v1p1_ModelPara.b[i],m_v1p1_ModelPara.b2[i]);
				}
				else
				{
					myc2d(c_v1p1_ModelPara.a[i],c_v1p1_ModelPara.b[i],c_v1p1_ModelPara.na,c_v1p1_ModelPara.nb,m_v1p1_ModelPara.T/1000,m_v1p1_ModelPara.a[i],m_v1p1_ModelPara.b[i]);
				}
				m_v1p1_ModelPara.d[i]=c_v1p1_ModelPara.d[i]/(m_v1p1_ModelPara.T/1000);
				if(c_v1p1_ModelPara.nb2!=0)
				{
					m_v1p1_ModelPara.d2[i]=c_v1p1_ModelPara.d2[i]/(m_v1p1_ModelPara.T/1000);
				}
				if(m_v1p1_ModelPara.d_max<m_v1p1_ModelPara.d[i])
					m_v1p1_ModelPara.d_max=m_v1p1_ModelPara.d[i];
				if(c_v1p1_ModelPara.nb2!=0)
				{
					if(m_v1p1_ModelPara.d2_max<m_v1p1_ModelPara.d2[i])
						m_v1p1_ModelPara.d2_max=m_v1p1_ModelPara.d2[i];
					m_v1p1_ModelPara.d_max = max(m_v1p1_ModelPara.d_max,m_v1p1_ModelPara.d2_max);
				}
				double ta[2]={1,-1};
				conv(m_v1p1_ModelPara.a[i],ta,m_v1p1_ModelPara.na,2,m_v1p1_ModelPara.aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
				if(c_v1p1_ModelPara.NLTDflag==1)
				{
					myc2d(c_v1p1_ModelPara.NLTD_a[i],c_v1p1_ModelPara.NLTD_b[i],c_v1p1_ModelPara.NLTDna,c_v1p1_ModelPara.NLTDnb,m_v1p1_ModelPara.T/1000,m_v1p1_ModelPara.NLTD_a[i],m_v1p1_ModelPara.NLTD_b[i]);
					m_v1p1_ModelPara.NLTD_d[i]=c_v1p1_ModelPara.NLTD_d[i]/(m_v1p1_ModelPara.T/1000);
					if(m_v1p1_ModelPara.NLTD_d_max<m_v1p1_ModelPara.NLTD_d[i])
						m_v1p1_ModelPara.NLTD_d_max=m_v1p1_ModelPara.NLTD_d[i];
					conv(m_v1p1_ModelPara.NLTD_a[i],ta,m_v1p1_ModelPara.NLTDna,2,m_v1p1_ModelPara.NLTD_aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
				}
				
			}
				return true;
		}
		else{
			//cout<< "input wrong!" << endl;
			return false;		
		}
	}
	else
	{
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
		printf("%s",StrBuf);
		return false;
	}
}

/*==============================================����һ========================================================*/
//�������ļ���ȡϵͳģ��
//�ļ���Ϊ��MgpcModel+���+.para;���Ϊ6λ���֣�ǰ��λΪҳ�ţ�������λ��0������λΪID�ţ�������λ��0;
//���磺��5ҳ���㷨IDΪ3�����Ϊ��005003�������ļ���Ϊ��MgpcModel005003.para;
//�ļ�Ŀ¼������Ŀ¼\domain013\data\ctrl\obj\drop001
bool AlgICSMIMOGPC::ReadParav1p2FromFile()//Ҫ�� MR ?
{	
	char filedir[]="\/edpf\/prj\/running\/ctrl\/obj";//���ļ�Ŀ¼
    //char filedir2[]="E:\\project\\DPU59running\\ctrl\\obj";//������DPU�ļ�Ŀ¼
	char filedir2[]="C:\\Program Files\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
    char filedir3[]="C:\\Program Files (x86)\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
	//char filedir2[]="C:\\Projects\\edpf\\MGPCMatrix";//������DPU�ļ�Ŀ¼
	AlgIdType AlgId,SheetNo,m_RecId,FinalId;
	AlgId=GetId();//��ȡ�㷨Id
	SheetNo=AlgId>>16;
	m_RecId=AlgId&0x000000FF;
	SheetNo=AlgId>>16;

	short int	St;
	char	filename[MAX_PATH]="";
	char	FilePath[MAX_PATH]="";
	char	StrBuf[256];
	size_t	Len=MAX_PATH;
	short int	rc;
	short int	j;

	FILE *fp=NULL;
	char str[MAX_PATH];
	char *p;
	char *delims={" []=\n"};
	// �����ļ���
	sprintf(filename,"Mgpcv1p2Model%03d%03d.para",SheetNo,m_RecId);
	sprintf(FilePath,"%s/%s",filedir,filename);
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(FilePath,"%s\\%s",filedir2,filename);
		if((fp=fopen(FilePath,"r"))==NULL)
		{
			sprintf(FilePath,"%s\\%s",filedir3,filename);
		    if((fp=fopen(FilePath,"r"))==NULL)
		    {
			    sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
			    printf("%s",StrBuf);
			    AlgState=AlgState|0x0001;
			    return false;
		    }
		}
	}

	AlgState=AlgState&0xFFFE;

	for(int i=0;i<7;i++)
		c_v1p2_ModelPara.flag[i]=0;
	
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
		printf("%s",StrBuf);
		AlgState=AlgState|0x0001;
		return false;
	}

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"totalmodel=",11)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v1p2_ModelPara.totalmodel=strtod(p,NULL)-1;
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v1p2_ModelPara.totalmodel>(M1MIMOGPC_MODELNUM_MAX-1)||c_v1p2_ModelPara.totalmodel<0){
		c_v1p2_ModelPara.flag[0]=0;
		rewind(fp);
	}
	else{
		c_v1p2_ModelPara.flag[0]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"na=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v1p2_ModelPara.na=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v1p2_ModelPara.na>M1MIMOGPC_MODEL_NA_MAX||c_v1p2_ModelPara.na<2)
	{
		c_v1p2_ModelPara.flag[1]=0;
		rewind(fp);
	}
	else
	{
		c_v1p2_ModelPara.flag[1]=1;
	}	
	
	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v1p2_ModelPara.nb=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v1p2_ModelPara.nb>M1MIMOGPC_MODEL_NB_MAX||c_v1p2_ModelPara.nb<1||c_v1p2_ModelPara.nb>c_v1p2_ModelPara.na)
	{
		c_v1p2_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_v1p2_ModelPara.flag[2]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb2=",4)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v1p2_ModelPara.nb2=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v1p2_ModelPara.nb2>M1MIMOGPC_MODEL_NB_MAX||c_v1p2_ModelPara.nb2>c_v1p2_ModelPara.na)
	{
		c_v1p2_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_v1p2_ModelPara.flag[2]=1;
	}

	AlgState=AlgState&0xFFFE;

	if (c_v1p2_ModelPara.flag[0]&&c_v1p2_ModelPara.flag[1]&&c_v1p2_ModelPara.flag[2])//must read "totalmodel" and "na" before
	{
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(str[0]=='A')
			{
				c_v1p2_ModelPara.flag[4]=1;
				for(int i=0;i<=c_v1p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_v1p2_ModelPara.na))
						  {
							c_v1p2_ModelPara.a[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_v1p2_ModelPara.na)//check a[i][0]=1;
							c_v1p2_ModelPara.flag[4]=0;
					}
					else
					{
						c_v1p2_ModelPara.flag[4]=0;
					}

				}
			}
			if(strncmp(str,"B=",2)==0)
			{
				c_v1p2_ModelPara.flag[5]=1;
				for(int i=0;i<=c_v1p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_v1p2_ModelPara.nb))
						  {
							c_v1p2_ModelPara.b[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_v1p2_ModelPara.nb)
							c_v1p2_ModelPara.flag[5]=0;
					}
					else
					{
						c_v1p2_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"B2=",3)==0)
			{
				c_v1p2_ModelPara.flag[5]=1;
				for(int i=0;i<=c_v1p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<c_v1p2_ModelPara.nb2))
						{
							c_v1p2_ModelPara.b2[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j!=c_v1p2_ModelPara.nb2)//check a[i][0]=1;
							c_v1p2_ModelPara.flag[5]=0;
					}
					else
					{
						c_v1p2_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"D=",2)==0)
			{
				c_v1p2_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<=c_v1p2_ModelPara.totalmodel))
						{
							c_v1p2_ModelPara.d[j]=max(m_TunePara.Ts,strtod(p,NULL));
							if (c_v1p2_ModelPara.d[j]<1)
							{
								c_v1p2_ModelPara.d[j]=1;
								//c1_ModelPara.flag[6]=0;
							}
							if(c_v1p2_ModelPara.d_max<c_v1p2_ModelPara.d[j])
								c_v1p2_ModelPara.d_max=c_v1p2_ModelPara.d[j];
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j<c_v1p2_ModelPara.totalmodel)
							c_v1p2_ModelPara.flag[6]=0;
						//break;
				}
				else
				{
					c_v1p2_ModelPara.flag[6]=0;
				}
			}
			if(strncmp(str,"D2=",3)==0)
			{
				c_v1p2_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
					p=strtok(str,delims);
					j=0;
					while((p!=NULL)&&(j<=c_v1p2_ModelPara.totalmodel))
					{
						c_v1p2_ModelPara.d2[j]=max(m_TunePara.Ts,strtod(p,NULL));
						if (c_v1p2_ModelPara.d2[j]<1)
						{
							c_v1p2_ModelPara.d2[j]=1;
							//	c_ModelPara.flag[6]=0;
						}
						if(c_v1p2_ModelPara.d2_max<c_v1p2_ModelPara.d2[j])
							c_v1p2_ModelPara.d2_max=c_v1p2_ModelPara.d2[j];
						j++;
						p=strtok(NULL,delims); 
					} 
					if (j<c_v1p2_ModelPara.totalmodel)
						c_v1p2_ModelPara.flag[6]=0;
					break;
				}
				else
				{
					c_v1p2_ModelPara.flag[6]=0;
				}
			}
			//fputs(str,stdout);
		}
//��NLTD����
		rc=1;
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(strncmp(str,"NLTDFLAG=",9)==0)
			{
				p=strtok(str,delims);
				if(p!=NULL)
				{
					p=strtok(NULL,delims);
					c_v1p2_ModelPara.NLTDflag=strtod(p,NULL);
					rc=0;
					break;
				}
			}
		}
		if (rc||c_v1p2_ModelPara.NLTDflag>1||c_v1p2_ModelPara.NLTDflag<0)
		{
			c_v1p2_ModelPara.flag[0]=0;
			rewind(fp);
		}
		else
		{
			c_v1p2_ModelPara.flag[0]=1;
		}
		if(c_v1p2_ModelPara.NLTDflag==1)
		{
			rc=1;
			while(fgets(str,sizeof(str),fp)!=NULL)
			{
				if(strncmp(str,"NLTDna=",7)==0)
				{
					p=strtok(str,delims);
					if(p!=NULL)
					{
						p=strtok(NULL,delims);
						c_v1p2_ModelPara.NLTDna=strtod(p,NULL);
						rc=0;
						break;
					}
				}
			}
			if (rc||c_v1p2_ModelPara.NLTDna>M1MIMOGPC_MODEL_NA_MAX||c_v1p2_ModelPara.NLTDna<2)
			{
				c_v1p2_ModelPara.flag[1]=0;
				rewind(fp);
			}
			else
			{
				c_v1p2_ModelPara.flag[1]=1;
			}	
	
			rc=1;
			while(fgets(str,sizeof(str),fp)!=NULL)
			{
				if(strncmp(str,"NLTDnb=",7)==0)
				{
					p=strtok(str,delims);
					if(p!=NULL)
					{
						p=strtok(NULL,delims);
						c_v1p2_ModelPara.NLTDnb=strtod(p,NULL);
						rc=0;
						break;
					}
				}
			}
			if (rc||c_v1p2_ModelPara.NLTDnb>M1MIMOGPC_MODEL_NB_MAX||c_v1p2_ModelPara.NLTDnb<1||c_v1p2_ModelPara.NLTDnb>c_v1p2_ModelPara.NLTDna)
			{
				c_v1p2_ModelPara.flag[2]=0;
				rewind(fp);
			}
			else
			{
				c_v1p2_ModelPara.flag[2]=1;
			}	
			if (c_v1p2_ModelPara.flag[0]&&c_v1p2_ModelPara.flag[1]&&c_v1p2_ModelPara.flag[2])
			{
				while(fgets(str,sizeof(str),fp)!=NULL)
				{
					if(strncmp(str,"NLTDA=",6)==0)
					{
						c_v1p2_ModelPara.flag[4]=1;
						for(int i=0;i<=c_v1p2_ModelPara.totalmodel;i++)
						{
							if(fgets(str,sizeof(str),fp)!=NULL)
							{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<c_v1p2_ModelPara.NLTDna))
								{
									c_v1p2_ModelPara.NLTD_a[i][j]=strtod(p,NULL);
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j!=c_v1p2_ModelPara.NLTDna)//check a[i][0]=1;
									c_v1p2_ModelPara.flag[4]=0;
							}
							else
							{
								c_v1p2_ModelPara.flag[4]=0;
							}

						}
					}
					if(strncmp(str,"NLTDB=",6)==0)
					{
						c_v1p2_ModelPara.flag[5]=1;
						for(int i=0;i<=c_v1p2_ModelPara.totalmodel;i++)
						{
							if(fgets(str,sizeof(str),fp)!=NULL)
							{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<c_v1p2_ModelPara.NLTDnb))
								{
									c_v1p2_ModelPara.NLTD_b[i][j]=strtod(p,NULL);
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j!=c_v1p2_ModelPara.NLTDnb)
									c_v1p2_ModelPara.flag[5]=0;
							}
							else
							{
								c_v1p2_ModelPara.flag[5]=0;
							}
						}
					}
					if(strncmp(str,"NLTDD=",6)==0)
					{
						c_v1p2_ModelPara.flag[6]=1;
						if(fgets(str,sizeof(str),fp)!=NULL)
						{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<=c_v1p2_ModelPara.totalmodel))
								{
									c_v1p2_ModelPara.NLTD_d[j]=max(m_TunePara.Ts,strtod(p,NULL));
									if (c_v1p2_ModelPara.NLTD_d[j]<1)
									{
										c_v1p2_ModelPara.NLTD_d[j]=1;
										//c1_ModelPara.flag[6]=0;
									}
									if(c_v1p2_ModelPara.NLTD_d_max<c_v1p2_ModelPara.NLTD_d[j])
										c_v1p2_ModelPara.NLTD_d_max=c_v1p2_ModelPara.NLTD_d[j];
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j<c_v1p2_ModelPara.totalmodel)
									c_v1p2_ModelPara.flag[6]=0;
						}
						else
						{
							c_v1p2_ModelPara.flag[6]=0;
						}
					}
				}
			}
			else
			{
				fclose(fp);
				sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
				printf("%s",StrBuf);
				return false;
			}
		
		}
		
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Load OK.\n",FilePath);
		printf("%s",StrBuf);
		c_v1p2_ModelPara.T=m_TunePara.Ts;
		c_v1p2_ModelPara.flag[3]=1;//T���������ļ�����
		j=0;//AlgState=AlgState&0xFFFE;AlgState=AlgState|0x0001;		
		for (int i=0;i<7;i++){
			m_v1p2_ModelPara.flag[i]=c_v1p2_ModelPara.flag[i];
			if (1==c_v1p2_ModelPara.flag[i])
			{
				j++;
				AlgState=AlgState&(~(0x0001<<(i+1)));
			}
			else
				AlgState=AlgState|0x0001<<(i+1);
		}

		if (7==j)
		{
			//����ϵͳ��ɢ��na>nb
			m_v1p2_ModelPara.totalmodel=c_v1p2_ModelPara.totalmodel;
			m_v1p2_ModelPara.na=c_v1p2_ModelPara.na;
			m_v1p2_ModelPara.nb=c_v1p2_ModelPara.na;
			m_v1p2_ModelPara.nb2=c_v1p2_ModelPara.na;
			m_v1p2_ModelPara.T=c_v1p2_ModelPara.T*1000;//����Ϊ��λ
			m_v1p2_ModelPara.d_max=1;
			m_v1p2_ModelPara.d2_max=1;

			m_v1p2_ModelPara.NLTDna=c_v1p2_ModelPara.NLTDna;
			m_v1p2_ModelPara.NLTDnb=c_v1p2_ModelPara.NLTDna;
			m_v1p2_ModelPara.NLTD_d_max=1;

			for(int i=0;i<=m_v1p2_ModelPara.totalmodel;i++)
			{
				if(c_v1p2_ModelPara.nb2!=0)
				{
					myc2d2(c_v1p2_ModelPara.a[i],c_v1p2_ModelPara.b[i],c_v1p2_ModelPara.b2[i],c_v1p2_ModelPara.na,c_v1p2_ModelPara.nb,c_v1p2_ModelPara.nb2,m_v1p2_ModelPara.T/1000,m_v1p2_ModelPara.a[i],m_v1p2_ModelPara.b[i],m_v1p2_ModelPara.b2[i]);
				}
				else
				{
					myc2d(c_v1p2_ModelPara.a[i],c_v1p2_ModelPara.b[i],c_v1p2_ModelPara.na,c_v1p2_ModelPara.nb,m_v1p2_ModelPara.T/1000,m_v1p2_ModelPara.a[i],m_v1p2_ModelPara.b[i]);
				}
				m_v1p2_ModelPara.d[i]=c_v1p2_ModelPara.d[i]/(m_v1p2_ModelPara.T/1000);
				if(c_v1p2_ModelPara.nb2!=0)
				{
					m_v1p2_ModelPara.d2[i]=c_v1p2_ModelPara.d2[i]/(m_v1p2_ModelPara.T/1000);
				}
				if(m_v1p2_ModelPara.d_max<m_v1p2_ModelPara.d[i])
					m_v1p2_ModelPara.d_max=m_v1p2_ModelPara.d[i];
				if(c_v1p2_ModelPara.nb2!=0)
				{
					if(m_v1p2_ModelPara.d2_max<m_v1p2_ModelPara.d2[i])
						m_v1p2_ModelPara.d2_max=m_v1p2_ModelPara.d2[i];
					m_v1p2_ModelPara.d_max = max(m_v1p2_ModelPara.d_max,m_v1p2_ModelPara.d2_max);
				}
				double ta[2]={1,-1};
				conv(m_v1p2_ModelPara.a[i],ta,m_v1p2_ModelPara.na,2,m_v1p2_ModelPara.aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
				if(c_v1p2_ModelPara.NLTDflag==1)
				{
					myc2d(c_v1p2_ModelPara.NLTD_a[i],c_v1p2_ModelPara.NLTD_b[i],c_v1p2_ModelPara.NLTDna,c_v1p2_ModelPara.NLTDnb,m_v1p2_ModelPara.T/1000,m_v1p2_ModelPara.NLTD_a[i],m_v1p2_ModelPara.NLTD_b[i]);
					m_v1p2_ModelPara.NLTD_d[i]=c_v1p2_ModelPara.NLTD_d[i]/(m_v1p2_ModelPara.T/1000);
					if(m_v1p2_ModelPara.NLTD_d_max<m_v1p2_ModelPara.NLTD_d[i])
						m_v1p2_ModelPara.NLTD_d_max=m_v1p2_ModelPara.NLTD_d[i];
					conv(m_v1p2_ModelPara.NLTD_a[i],ta,m_v1p2_ModelPara.NLTDna,2,m_v1p2_ModelPara.NLTD_aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
				}
				
			}
				return true;
		}
		else{
			//cout<< "input wrong!" << endl;
			return false;		
		}
	}
	else
	{
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
		printf("%s",StrBuf);
		return false;
	}
}



/*==============================================���Ŷ�========================================================*/
//�������ļ���ȡϵͳģ��
//�ļ���Ϊ��MgpcModel+���+.para;���Ϊ6λ���֣�ǰ��λΪҳ�ţ�������λ��0������λΪID�ţ�������λ��0;
//���磺��5ҳ���㷨IDΪ3�����Ϊ��005003�������ļ���Ϊ��MgpcModel005003.para;
//�ļ�Ŀ¼������Ŀ¼\domain013\data\ctrl\obj\drop001
bool AlgICSMIMOGPC::ReadParav2p1FromFile()//Ҫ�� MR ?
{	
	char filedir[]="\/edpf\/prj\/running\/ctrl\/obj";//���ļ�Ŀ¼
    //char filedir2[]="E:\\project\\DPU59running\\ctrl\\obj";//������DPU�ļ�Ŀ¼
	char filedir2[]="C:\\Program Files\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
    char filedir3[]="C:\\Program Files (x86)\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
	//char filedir2[]="C:\\Projects\\edpf\\MGPCMatrix";//������DPU�ļ�Ŀ¼
	AlgIdType AlgId,SheetNo,m_RecId,FinalId;
	AlgId=GetId();//��ȡ�㷨Id
	SheetNo=AlgId>>16;
	m_RecId=AlgId&0x000000FF;
	SheetNo=AlgId>>16;

	short int	St;
	char	filename[MAX_PATH]="";
	char	FilePath[MAX_PATH]="";
	char	StrBuf[256];
	size_t	Len=MAX_PATH;
	short int	rc;
	short int	j;

	FILE *fp=NULL;
	char str[MAX_PATH];
	char *p;
	char *delims={" []=\n"};
	// �����ļ���
	sprintf(filename,"Mgpcv2p1Model%03d%03d.para",SheetNo,m_RecId);
	sprintf(FilePath,"%s/%s",filedir,filename);
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(FilePath,"%s\\%s",filedir2,filename);
		if((fp=fopen(FilePath,"r"))==NULL)
		{
			sprintf(FilePath,"%s\\%s",filedir3,filename);
		    if((fp=fopen(FilePath,"r"))==NULL)
		    {
			    sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
			    printf("%s",StrBuf);
			    AlgState=AlgState|0x0001;
			    return false;
		    }
		}
	}

	AlgState=AlgState&0xFFFE;

	for(int i=0;i<7;i++)
		c_v2p1_ModelPara.flag[i]=0;
	
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
		printf("%s",StrBuf);
		AlgState=AlgState|0x0001;
		return false;
	}

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"totalmodel=",11)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v2p1_ModelPara.totalmodel=strtod(p,NULL)-1;
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v2p1_ModelPara.totalmodel>(M2MIMOGPC_MODELNUM_MAX-1)||c_v2p1_ModelPara.totalmodel<0){
		c_v2p1_ModelPara.flag[0]=0;
		rewind(fp);
	}
	else{
		c_v2p1_ModelPara.flag[0]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"na=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v2p1_ModelPara.na=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v2p1_ModelPara.na>M2MIMOGPC_MODEL_NA_MAX||c_v2p1_ModelPara.na<2)
	{
		c_v2p1_ModelPara.flag[1]=0;
		rewind(fp);
	}
	else
	{
		c_v2p1_ModelPara.flag[1]=1;
	}	
	
	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v2p1_ModelPara.nb=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v2p1_ModelPara.nb>M2MIMOGPC_MODEL_NB_MAX||c_v2p1_ModelPara.nb<1||c_v2p1_ModelPara.nb>c_v2p1_ModelPara.na)
	{
		c_v2p1_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_v2p1_ModelPara.flag[2]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb2=",4)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v2p1_ModelPara.nb2=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v2p1_ModelPara.nb2>M2MIMOGPC_MODEL_NB_MAX||c_v2p1_ModelPara.nb2>c_v2p1_ModelPara.na)
	{
		c_v2p1_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_v2p1_ModelPara.flag[2]=1;
	}

	AlgState=AlgState&0xFFFE;

	if (c_v2p1_ModelPara.flag[0]&&c_v2p1_ModelPara.flag[1]&&c_v2p1_ModelPara.flag[2])//must read "totalmodel" and "na" before
	{
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(str[0]=='A')
			{
				c_v2p1_ModelPara.flag[4]=1;
				for(int i=0;i<=c_v2p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_v2p1_ModelPara.na))
						  {
							c_v2p1_ModelPara.a[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_v2p1_ModelPara.na)//check a[i][0]=1;
							c_v2p1_ModelPara.flag[4]=0;
					}
					else
					{
						c_v2p1_ModelPara.flag[4]=0;
					}

				}
			}
			if(strncmp(str,"B=",2)==0)
			{
				c_v2p1_ModelPara.flag[5]=1;
				for(int i=0;i<=c_v2p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_v2p1_ModelPara.nb))
						  {
							c_v2p1_ModelPara.b[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_v2p1_ModelPara.nb)
							c_v2p1_ModelPara.flag[5]=0;
					}
					else
					{
						c_v2p1_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"B2=",3)==0)
			{
				c_v2p1_ModelPara.flag[5]=1;
				for(int i=0;i<=c_v2p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<c_v2p1_ModelPara.nb2))
						{
							c_v2p1_ModelPara.b2[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j!=c_v2p1_ModelPara.nb2)//check a[i][0]=1;
							c_v2p1_ModelPara.flag[5]=0;
					}
					else
					{
						c_v2p1_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"D=",2)==0)
			{
				c_v2p1_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<=c_v2p1_ModelPara.totalmodel))
						{
							c_v2p1_ModelPara.d[j]=max(m_TunePara.Ts,strtod(p,NULL));
							if (c_v2p1_ModelPara.d[j]<1)
							{
								c_v2p1_ModelPara.d[j]=1;
							//	c2_ModelPara.flag[6]=0;
							}
							if(c_v2p1_ModelPara.d_max<c_v2p1_ModelPara.d[j])
								c_v2p1_ModelPara.d_max=c_v2p1_ModelPara.d[j];
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j<c_v2p1_ModelPara.totalmodel)
							c_v2p1_ModelPara.flag[6]=0;
						//break;
				}
				else
				{
					c_v2p1_ModelPara.flag[6]=0;
				}
			}
			if(strncmp(str,"D2=",3)==0)
			{
				c_v2p1_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
					p=strtok(str,delims);
					j=0;
					while((p!=NULL)&&(j<=c_v2p1_ModelPara.totalmodel))
					{
						c_v2p1_ModelPara.d2[j]=max(m_TunePara.Ts,strtod(p,NULL));
						if (c_v2p1_ModelPara.d2[j]<1)
						{
							c_v2p1_ModelPara.d2[j]=1;
							//	c_ModelPara.flag[6]=0;
						}
						if(c_v2p1_ModelPara.d2_max<c_v2p1_ModelPara.d2[j])
							c_v2p1_ModelPara.d2_max=c_v2p1_ModelPara.d2[j];
						j++;
						p=strtok(NULL,delims); 
					} 
					if (j<c_v2p1_ModelPara.totalmodel)
						c_v2p1_ModelPara.flag[6]=0;
					break;
				}
				else
				{
					c_v2p1_ModelPara.flag[6]=0;
				}
			}
			//fputs(str,stdout);
		}
		//��NLTD����
		rc=1;
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(strncmp(str,"NLTDFLAG=",9)==0)
			{
				p=strtok(str,delims);
				if(p!=NULL)
				{
					p=strtok(NULL,delims);
					c_v2p1_ModelPara.NLTDflag=strtod(p,NULL);
					rc=0;
					break;
				}
			}
		}
		if (rc||c_v2p1_ModelPara.NLTDflag>1||c_v2p1_ModelPara.NLTDflag<0)
		{
			c_v2p1_ModelPara.flag[0]=0;
			rewind(fp);
		}
		else
		{
			c_v2p1_ModelPara.flag[0]=1;
		}
		if(c_v2p1_ModelPara.NLTDflag==1)
		{
			rc=1;
			while(fgets(str,sizeof(str),fp)!=NULL)
			{
				if(strncmp(str,"NLTDna=",7)==0)
				{
					p=strtok(str,delims);
					if(p!=NULL)
					{
						p=strtok(NULL,delims);
						c_v2p1_ModelPara.NLTDna=strtod(p,NULL);
						rc=0;
						break;
					}
				}
			}
			if (rc||c_v2p1_ModelPara.NLTDna>M2MIMOGPC_MODEL_NA_MAX||c_v2p1_ModelPara.NLTDna<2)
			{
				c_v2p1_ModelPara.flag[1]=0;
				rewind(fp);
			}
			else
			{
				c_v2p1_ModelPara.flag[1]=1;
			}	
	
			rc=1;
			while(fgets(str,sizeof(str),fp)!=NULL)
			{
				if(strncmp(str,"NLTDnb=",7)==0)
				{
					p=strtok(str,delims);
					if(p!=NULL)
					{
						p=strtok(NULL,delims);
						c_v2p1_ModelPara.NLTDnb=strtod(p,NULL);
						rc=0;
						break;
					}
				}
			}
			if (rc||c_v2p1_ModelPara.NLTDnb>M2MIMOGPC_MODEL_NB_MAX||c_v2p1_ModelPara.NLTDnb<1||c_v2p1_ModelPara.NLTDnb>c_v2p1_ModelPara.NLTDna)
			{
				c_v2p1_ModelPara.flag[2]=0;
				rewind(fp);
			}
			else
			{
				c_v2p1_ModelPara.flag[2]=1;
			}	
			if (c_v2p1_ModelPara.flag[0]&&c_v2p1_ModelPara.flag[1]&&c_v2p1_ModelPara.flag[2])
			{
				while(fgets(str,sizeof(str),fp)!=NULL)
				{
					if(strncmp(str,"NLTDA=",6)==0)
					{
						c_v2p1_ModelPara.flag[4]=1;
						for(int i=0;i<=c_v2p1_ModelPara.totalmodel;i++)
						{
							if(fgets(str,sizeof(str),fp)!=NULL)
							{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<c_v2p1_ModelPara.NLTDna))
								{
									c_v2p1_ModelPara.NLTD_a[i][j]=strtod(p,NULL);
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j!=c_v2p1_ModelPara.NLTDna)//check a[i][0]=1;
									c_v2p1_ModelPara.flag[4]=0;
							}
							else
							{
								c_v2p1_ModelPara.flag[4]=0;
							}

						}
					}
					if(strncmp(str,"NLTDB=",6)==0)
					{
						c_v2p1_ModelPara.flag[5]=1;
						for(int i=0;i<=c_v2p1_ModelPara.totalmodel;i++)
						{
							if(fgets(str,sizeof(str),fp)!=NULL)
							{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<c_v2p1_ModelPara.NLTDnb))
								{
									c_v2p1_ModelPara.NLTD_b[i][j]=strtod(p,NULL);
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j!=c_v2p1_ModelPara.NLTDnb)
									c_v2p1_ModelPara.flag[5]=0;
							}
							else
							{
								c_v2p1_ModelPara.flag[5]=0;
							}
						}
					}
					if(strncmp(str,"NLTDD=",6)==0)
					{
						c_v2p1_ModelPara.flag[6]=1;
						if(fgets(str,sizeof(str),fp)!=NULL)
						{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<=c_v2p1_ModelPara.totalmodel))
								{
									c_v2p1_ModelPara.NLTD_d[j]=max(m_TunePara.Ts,strtod(p,NULL));
									if (c_v2p1_ModelPara.NLTD_d[j]<1)
									{
										c_v2p1_ModelPara.NLTD_d[j]=1;
										//c2_ModelPara.flag[6]=0;
									}
									if(c_v2p1_ModelPara.NLTD_d_max<c_v2p1_ModelPara.NLTD_d[j])
										c_v2p1_ModelPara.NLTD_d_max=c_v2p1_ModelPara.NLTD_d[j];
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j<c_v2p1_ModelPara.totalmodel)
									c_v2p1_ModelPara.flag[6]=0;
						}
						else
						{
							c_v2p1_ModelPara.flag[6]=0;
						}
					}
				}
			}
			else
			{
				fclose(fp);
				sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
				printf("%s",StrBuf);
				return false;
			}
		
		}



		fclose(fp);
		sprintf(StrBuf,"Para File: %s Load OK.\n",FilePath);
		printf("%s",StrBuf);
		c_v2p1_ModelPara.T=m_TunePara.Ts;
		c_v2p1_ModelPara.flag[3]=1;//T���������ļ�����
		j=0;//AlgState=AlgState&0xFFFE;AlgState=AlgState|0x0001;		
		for (int i=0;i<7;i++){
			m_v2p1_ModelPara.flag[i]=c_v2p1_ModelPara.flag[i];
			if (1==c_v2p1_ModelPara.flag[i])
			{
				j++;
				AlgState=AlgState&(~(0x0001<<(i+1)));
			}
			else
				AlgState=AlgState|0x0001<<(i+1);
		}

		if (7==j){
			//����ϵͳ��ɢ��na>nb
			m_v2p1_ModelPara.totalmodel=c_v2p1_ModelPara.totalmodel;
			m_v2p1_ModelPara.na=c_v2p1_ModelPara.na;
			m_v2p1_ModelPara.nb=c_v2p1_ModelPara.na;
			m_v2p1_ModelPara.nb2=c_v2p1_ModelPara.na;
			m_v2p1_ModelPara.T=c_v2p1_ModelPara.T*1000;//����Ϊ��λ
			m_v2p1_ModelPara.d_max=1;
			m_v2p1_ModelPara.d2_max=1;

			m_v2p1_ModelPara.NLTDna=c_v2p1_ModelPara.NLTDna;
			m_v2p1_ModelPara.NLTDnb=c_v2p1_ModelPara.NLTDna;
			m_v2p1_ModelPara.NLTD_d_max=1;

			for(int i=0;i<=m_v2p1_ModelPara.totalmodel;i++)
			{
				if(c_v2p1_ModelPara.nb2!=0)
				{
					myc2d2(c_v2p1_ModelPara.a[i],c_v2p1_ModelPara.b[i],c_v2p1_ModelPara.b2[i],c_v2p1_ModelPara.na,c_v2p1_ModelPara.nb,c_v2p1_ModelPara.nb2,m_v2p1_ModelPara.T/1000,m_v2p1_ModelPara.a[i],m_v2p1_ModelPara.b[i],m_v2p1_ModelPara.b2[i]);
				}
				else
				{
					myc2d(c_v2p1_ModelPara.a[i],c_v2p1_ModelPara.b[i],c_v2p1_ModelPara.na,c_v2p1_ModelPara.nb,m_v2p1_ModelPara.T/1000,m_v2p1_ModelPara.a[i],m_v2p1_ModelPara.b[i]);
				}
				m_v2p1_ModelPara.d[i]=c_v2p1_ModelPara.d[i]/(m_v2p1_ModelPara.T/1000);
				if(c_v2p1_ModelPara.nb2!=0)
				{
					m_v2p1_ModelPara.d2[i]=c_v2p1_ModelPara.d2[i]/(m_v2p1_ModelPara.T/1000);
				}
				if(m_v2p1_ModelPara.d_max<m_v2p1_ModelPara.d[i])
					m_v2p1_ModelPara.d_max=m_v2p1_ModelPara.d[i];
				if(c_v2p1_ModelPara.nb2!=0)
				{
					if(m_v2p1_ModelPara.d2_max<m_v2p1_ModelPara.d2[i])
						m_v2p1_ModelPara.d2_max=m_v2p1_ModelPara.d2[i];
					m_v2p1_ModelPara.d_max = max(m_v2p1_ModelPara.d_max,m_v2p1_ModelPara.d2_max);
				}
				double ta[2]={1,-1};
				conv(m_v2p1_ModelPara.a[i],ta,m_v2p1_ModelPara.na,2,m_v2p1_ModelPara.aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
				if(c_v2p1_ModelPara.NLTDflag==1)
				{
					myc2d(c_v2p1_ModelPara.NLTD_a[i],c_v2p1_ModelPara.NLTD_b[i],c_v2p1_ModelPara.NLTDna,c_v2p1_ModelPara.NLTDnb,m_v2p1_ModelPara.T/1000,m_v2p1_ModelPara.NLTD_a[i],m_v2p1_ModelPara.NLTD_b[i]);
					m_v2p1_ModelPara.NLTD_d[i]=c_v2p1_ModelPara.NLTD_d[i]/(m_v2p1_ModelPara.T/1000);
					if(m_v2p1_ModelPara.NLTD_d_max<m_v2p1_ModelPara.NLTD_d[i])
						m_v2p1_ModelPara.NLTD_d_max=m_v2p1_ModelPara.NLTD_d[i];
					conv(m_v2p1_ModelPara.NLTD_a[i],ta,m_v2p1_ModelPara.NLTDna,2,m_v2p1_ModelPara.NLTD_aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
				}
				

			}
				return true;
		}
		else{
			//cout<< "input wrong!" << endl;
			return false;		
		}
	}
	else
	{
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
		printf("%s",StrBuf);
		return false;
	}
}


/*==============================================���Ŷ�========================================================*/
//�������ļ���ȡϵͳģ��
//�ļ���Ϊ��MgpcModel+���+.para;���Ϊ6λ���֣�ǰ��λΪҳ�ţ�������λ��0������λΪID�ţ�������λ��0;
//���磺��5ҳ���㷨IDΪ3�����Ϊ��005003�������ļ���Ϊ��MgpcModel005003.para;
//�ļ�Ŀ¼������Ŀ¼\domain013\data\ctrl\obj\drop001
bool AlgICSMIMOGPC::ReadParav2p2FromFile()//Ҫ�� MR ?
{	
	char filedir[]="\/edpf\/prj\/running\/ctrl\/obj";//���ļ�Ŀ¼
    //char filedir2[]="E:\\project\\DPU59running\\ctrl\\obj";//������DPU�ļ�Ŀ¼
	char filedir2[]="C:\\Program Files\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
    char filedir3[]="C:\\Program Files (x86)\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
	//char filedir2[]="C:\\Projects\\edpf\\MGPCMatrix";//������DPU�ļ�Ŀ¼
	AlgIdType AlgId,SheetNo,m_RecId,FinalId;
	AlgId=GetId();//��ȡ�㷨Id
	SheetNo=AlgId>>16;
	m_RecId=AlgId&0x000000FF;
	SheetNo=AlgId>>16;

	short int	St;
	char	filename[MAX_PATH]="";
	char	FilePath[MAX_PATH]="";
	char	StrBuf[256];
	size_t	Len=MAX_PATH;
	short int	rc;
	short int	j;

	FILE *fp=NULL;
	char str[MAX_PATH];
	char *p;
	char *delims={" []=\n"};
	// �����ļ���
	sprintf(filename,"Mgpcv2p2Model%03d%03d.para",SheetNo,m_RecId);
	sprintf(FilePath,"%s/%s",filedir,filename);
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(FilePath,"%s\\%s",filedir2,filename);
		if((fp=fopen(FilePath,"r"))==NULL)
		{
			sprintf(FilePath,"%s\\%s",filedir3,filename);
		    if((fp=fopen(FilePath,"r"))==NULL)
		    {
			    sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
			    printf("%s",StrBuf);
			    AlgState=AlgState|0x0001;
			    return false;
		    }
		}
	}

	AlgState=AlgState&0xFFFE;

	for(int i=0;i<7;i++)
		c_v2p2_ModelPara.flag[i]=0;
	
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
		printf("%s",StrBuf);
		AlgState=AlgState|0x0001;
		return false;
	}

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"totalmodel=",11)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v2p2_ModelPara.totalmodel=strtod(p,NULL)-1;
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v2p2_ModelPara.totalmodel>(M2MIMOGPC_MODELNUM_MAX-1)||c_v2p2_ModelPara.totalmodel<0){
		c_v2p2_ModelPara.flag[0]=0;
		rewind(fp);
	}
	else{
		c_v2p2_ModelPara.flag[0]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"na=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v2p2_ModelPara.na=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v2p2_ModelPara.na>M2MIMOGPC_MODEL_NA_MAX||c_v2p2_ModelPara.na<2)
	{
		c_v2p2_ModelPara.flag[1]=0;
		rewind(fp);
	}
	else
	{
		c_v2p2_ModelPara.flag[1]=1;
	}	
	
	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v2p2_ModelPara.nb=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v2p2_ModelPara.nb>M2MIMOGPC_MODEL_NB_MAX||c_v2p2_ModelPara.nb<1||c_v2p2_ModelPara.nb>c_v2p2_ModelPara.na)
	{
		c_v2p2_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_v2p2_ModelPara.flag[2]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb2=",4)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v2p2_ModelPara.nb2=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v2p2_ModelPara.nb2>M2MIMOGPC_MODEL_NB_MAX||c_v2p2_ModelPara.nb2>c_v2p2_ModelPara.na)
	{
		c_v2p2_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_v2p2_ModelPara.flag[2]=1;
	}

	AlgState=AlgState&0xFFFE;

	if (c_v2p2_ModelPara.flag[0]&&c_v2p2_ModelPara.flag[1]&&c_v2p2_ModelPara.flag[2])//must read "totalmodel" and "na" before
	{
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(str[0]=='A')
			{
				c_v2p2_ModelPara.flag[4]=1;
				for(int i=0;i<=c_v2p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_v2p2_ModelPara.na))
						  {
							c_v2p2_ModelPara.a[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_v2p2_ModelPara.na)//check a[i][0]=1;
							c_v2p2_ModelPara.flag[4]=0;
					}
					else
					{
						c_v2p2_ModelPara.flag[4]=0;
					}

				}
			}
			if(strncmp(str,"B=",2)==0)
			{
				c_v2p2_ModelPara.flag[5]=1;
				for(int i=0;i<=c_v2p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_v2p2_ModelPara.nb))
						  {
							c_v2p2_ModelPara.b[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_v2p2_ModelPara.nb)
							c_v2p2_ModelPara.flag[5]=0;
					}
					else
					{
						c_v2p2_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"B2=",3)==0)
			{
				c_v2p2_ModelPara.flag[5]=1;
				for(int i=0;i<=c_v2p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<c_v2p2_ModelPara.nb2))
						{
							c_v2p2_ModelPara.b2[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j!=c_v2p2_ModelPara.nb2)//check a[i][0]=1;
							c_v2p2_ModelPara.flag[5]=0;
					}
					else
					{
						c_v2p2_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"D=",2)==0)
			{
				c_v2p2_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<=c_v2p2_ModelPara.totalmodel))
						{
							c_v2p2_ModelPara.d[j]=max(m_TunePara.Ts,strtod(p,NULL));
							if (c_v2p2_ModelPara.d[j]<1)
							{
								c_v2p2_ModelPara.d[j]=1;
							//	c2_ModelPara.flag[6]=0;
							}
							if(c_v2p2_ModelPara.d_max<c_v2p2_ModelPara.d[j])
								c_v2p2_ModelPara.d_max=c_v2p2_ModelPara.d[j];
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j<c_v2p2_ModelPara.totalmodel)
							c_v2p2_ModelPara.flag[6]=0;
						//break;
				}
				else
				{
					c_v2p2_ModelPara.flag[6]=0;
				}
			}
			if(strncmp(str,"D2=",3)==0)
			{
				c_v2p2_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
					p=strtok(str,delims);
					j=0;
					while((p!=NULL)&&(j<=c_v2p2_ModelPara.totalmodel))
					{
						c_v2p2_ModelPara.d2[j]=max(m_TunePara.Ts,strtod(p,NULL));
						if (c_v2p2_ModelPara.d2[j]<1)
						{
							c_v2p2_ModelPara.d2[j]=1;
							//	c_ModelPara.flag[6]=0;
						}
						if(c_v2p2_ModelPara.d2_max<c_v2p2_ModelPara.d2[j])
							c_v2p2_ModelPara.d2_max=c_v2p2_ModelPara.d2[j];
						j++;
						p=strtok(NULL,delims); 
					} 
					if (j<c_v2p2_ModelPara.totalmodel)
						c_v2p2_ModelPara.flag[6]=0;
					break;
				}
				else
				{
					c_v2p2_ModelPara.flag[6]=0;
				}
			}
			//fputs(str,stdout);
		}
		//��NLTD����
		rc=1;
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(strncmp(str,"NLTDFLAG=",9)==0)
			{
				p=strtok(str,delims);
				if(p!=NULL)
				{
					p=strtok(NULL,delims);
					c_v2p2_ModelPara.NLTDflag=strtod(p,NULL);
					rc=0;
					break;
				}
			}
		}
		if (rc||c_v2p2_ModelPara.NLTDflag>1||c_v2p2_ModelPara.NLTDflag<0)
		{
			c_v2p2_ModelPara.flag[0]=0;
			rewind(fp);
		}
		else
		{
			c_v2p2_ModelPara.flag[0]=1;
		}
		if(c_v2p2_ModelPara.NLTDflag==1)
		{
			rc=1;
			while(fgets(str,sizeof(str),fp)!=NULL)
			{
				if(strncmp(str,"NLTDna=",7)==0)
				{
					p=strtok(str,delims);
					if(p!=NULL)
					{
						p=strtok(NULL,delims);
						c_v2p2_ModelPara.NLTDna=strtod(p,NULL);
						rc=0;
						break;
					}
				}
			}
			if (rc||c_v2p2_ModelPara.NLTDna>M2MIMOGPC_MODEL_NA_MAX||c_v2p2_ModelPara.NLTDna<2)
			{
				c_v2p2_ModelPara.flag[1]=0;
				rewind(fp);
			}
			else
			{
				c_v2p2_ModelPara.flag[1]=1;
			}	
	
			rc=1;
			while(fgets(str,sizeof(str),fp)!=NULL)
			{
				if(strncmp(str,"NLTDnb=",7)==0)
				{
					p=strtok(str,delims);
					if(p!=NULL)
					{
						p=strtok(NULL,delims);
						c_v2p2_ModelPara.NLTDnb=strtod(p,NULL);
						rc=0;
						break;
					}
				}
			}
			if (rc||c_v2p2_ModelPara.NLTDnb>M2MIMOGPC_MODEL_NB_MAX||c_v2p2_ModelPara.NLTDnb<1||c_v2p2_ModelPara.NLTDnb>c_v2p2_ModelPara.NLTDna)
			{
				c_v2p2_ModelPara.flag[2]=0;
				rewind(fp);
			}
			else
			{
				c_v2p2_ModelPara.flag[2]=1;
			}	
			if (c_v2p2_ModelPara.flag[0]&&c_v2p2_ModelPara.flag[1]&&c_v2p2_ModelPara.flag[2])
			{
				while(fgets(str,sizeof(str),fp)!=NULL)
				{
					if(strncmp(str,"NLTDA=",6)==0)
					{
						c_v2p2_ModelPara.flag[4]=1;
						for(int i=0;i<=c_v2p2_ModelPara.totalmodel;i++)
						{
							if(fgets(str,sizeof(str),fp)!=NULL)
							{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<c_v2p2_ModelPara.NLTDna))
								{
									c_v2p2_ModelPara.NLTD_a[i][j]=strtod(p,NULL);
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j!=c_v2p2_ModelPara.NLTDna)//check a[i][0]=1;
									c_v2p2_ModelPara.flag[4]=0;
							}
							else
							{
								c_v2p2_ModelPara.flag[4]=0;
							}

						}
					}
					if(strncmp(str,"NLTDB=",6)==0)
					{
						c_v2p2_ModelPara.flag[5]=1;
						for(int i=0;i<=c_v2p2_ModelPara.totalmodel;i++)
						{
							if(fgets(str,sizeof(str),fp)!=NULL)
							{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<c_v2p2_ModelPara.NLTDnb))
								{
									c_v2p2_ModelPara.NLTD_b[i][j]=strtod(p,NULL);
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j!=c_v2p2_ModelPara.NLTDnb)
									c_v2p2_ModelPara.flag[5]=0;
							}
							else
							{
								c_v2p2_ModelPara.flag[5]=0;
							}
						}
					}
					if(strncmp(str,"NLTDD=",6)==0)
					{
						c_v2p2_ModelPara.flag[6]=1;
						if(fgets(str,sizeof(str),fp)!=NULL)
						{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<=c_v2p2_ModelPara.totalmodel))
								{
									c_v2p2_ModelPara.NLTD_d[j]=max(m_TunePara.Ts,strtod(p,NULL));
									if (c_v2p2_ModelPara.NLTD_d[j]<1)
									{
										c_v2p2_ModelPara.NLTD_d[j]=1;
										//c2_ModelPara.flag[6]=0;
									}
									if(c_v2p2_ModelPara.NLTD_d_max<c_v2p2_ModelPara.NLTD_d[j])
										c_v2p2_ModelPara.NLTD_d_max=c_v2p2_ModelPara.NLTD_d[j];
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j<c_v2p2_ModelPara.totalmodel)
									c_v2p2_ModelPara.flag[6]=0;
						}
						else
						{
							c_v2p2_ModelPara.flag[6]=0;
						}
					}
				}
			}
			else
			{
				fclose(fp);
				sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
				printf("%s",StrBuf);
				return false;
			}
		
		}



		fclose(fp);
		sprintf(StrBuf,"Para File: %s Load OK.\n",FilePath);
		printf("%s",StrBuf);
		c_v2p2_ModelPara.T=m_TunePara.Ts;
		c_v2p2_ModelPara.flag[3]=1;//T���������ļ�����
		j=0;//AlgState=AlgState&0xFFFE;AlgState=AlgState|0x0001;		
		for (int i=0;i<7;i++){
			m_v2p2_ModelPara.flag[i]=c_v2p2_ModelPara.flag[i];
			if (1==c_v2p2_ModelPara.flag[i])
			{
				j++;
				AlgState=AlgState&(~(0x0001<<(i+1)));
			}
			else
				AlgState=AlgState|0x0001<<(i+1);
		}

		if (7==j){
			//����ϵͳ��ɢ��na>nb
			m_v2p2_ModelPara.totalmodel=c_v2p2_ModelPara.totalmodel;
			m_v2p2_ModelPara.na=c_v2p2_ModelPara.na;
			m_v2p2_ModelPara.nb=c_v2p2_ModelPara.na;
			m_v2p2_ModelPara.nb2=c_v2p2_ModelPara.na;
			m_v2p2_ModelPara.T=c_v2p2_ModelPara.T*1000;//����Ϊ��λ
			m_v2p2_ModelPara.d_max=1;
			m_v2p2_ModelPara.d2_max=1;

			m_v2p2_ModelPara.NLTDna=c_v2p2_ModelPara.NLTDna;
			m_v2p2_ModelPara.NLTDnb=c_v2p2_ModelPara.NLTDna;
			m_v2p2_ModelPara.NLTD_d_max=1;

			for(int i=0;i<=m_v2p2_ModelPara.totalmodel;i++)
			{
				if(c_v2p2_ModelPara.nb2!=0)
				{
					myc2d2(c_v2p2_ModelPara.a[i],c_v2p2_ModelPara.b[i],c_v2p2_ModelPara.b2[i],c_v2p2_ModelPara.na,c_v2p2_ModelPara.nb,c_v2p2_ModelPara.nb2,m_v2p2_ModelPara.T/1000,m_v2p2_ModelPara.a[i],m_v2p2_ModelPara.b[i],m_v2p2_ModelPara.b2[i]);
				}
				else
				{
					myc2d(c_v2p2_ModelPara.a[i],c_v2p2_ModelPara.b[i],c_v2p2_ModelPara.na,c_v2p2_ModelPara.nb,m_v2p2_ModelPara.T/1000,m_v2p2_ModelPara.a[i],m_v2p2_ModelPara.b[i]);
				}
				m_v2p2_ModelPara.d[i]=c_v2p2_ModelPara.d[i]/(m_v2p2_ModelPara.T/1000);
				if(c_v2p2_ModelPara.nb2!=0)
				{
					m_v2p2_ModelPara.d2[i]=c_v2p2_ModelPara.d2[i]/(m_v2p2_ModelPara.T/1000);
				}
				if(m_v2p2_ModelPara.d_max<m_v2p2_ModelPara.d[i])
					m_v2p2_ModelPara.d_max=m_v2p2_ModelPara.d[i];
				if(c_v2p2_ModelPara.nb2!=0)
				{
					if(m_v2p2_ModelPara.d2_max<m_v2p2_ModelPara.d2[i])
						m_v2p2_ModelPara.d2_max=m_v2p2_ModelPara.d2[i];
					m_v2p2_ModelPara.d_max = max(m_v2p2_ModelPara.d_max,m_v2p2_ModelPara.d2_max);
				}
				double ta[2]={1,-1};
				conv(m_v2p2_ModelPara.a[i],ta,m_v2p2_ModelPara.na,2,m_v2p2_ModelPara.aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
				if(c_v2p2_ModelPara.NLTDflag==1)
				{
					myc2d(c_v2p2_ModelPara.NLTD_a[i],c_v2p2_ModelPara.NLTD_b[i],c_v2p2_ModelPara.NLTDna,c_v2p2_ModelPara.NLTDnb,m_v2p2_ModelPara.T/1000,m_v2p2_ModelPara.NLTD_a[i],m_v2p2_ModelPara.NLTD_b[i]);
					m_v2p2_ModelPara.NLTD_d[i]=c_v2p2_ModelPara.NLTD_d[i]/(m_v2p2_ModelPara.T/1000);
					if(m_v2p2_ModelPara.NLTD_d_max<m_v2p2_ModelPara.NLTD_d[i])
						m_v2p2_ModelPara.NLTD_d_max=m_v2p2_ModelPara.NLTD_d[i];
					conv(m_v2p2_ModelPara.NLTD_a[i],ta,m_v2p2_ModelPara.NLTDna,2,m_v2p2_ModelPara.NLTD_aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
				}
				

			}
				return true;
		}
		else{
			//cout<< "input wrong!" << endl;
			return false;		
		}
	}
	else
	{
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
		printf("%s",StrBuf);
		return false;
	}
}



/*==============================================������========================================================*/
//�������ļ���ȡϵͳģ��
//�ļ���Ϊ��MgpcModel+���+.para;���Ϊ6λ���֣�ǰ��λΪҳ�ţ�������λ��0������λΪID�ţ�������λ��0;
//���磺��5ҳ���㷨IDΪ3�����Ϊ��005003�������ļ���Ϊ��MgpcModel005003.para;
//�ļ�Ŀ¼������Ŀ¼\domain013\data\ctrl\obj\drop001
bool AlgICSMIMOGPC::ReadParav3p1FromFile()//Ҫ�� MR ?
{	
	char filedir[]="\/edpf\/prj\/running\/ctrl\/obj";//���ļ�Ŀ¼
    //char filedir2[]="E:\\project\\DPU59running\\ctrl\\obj";//������DPU�ļ�Ŀ¼
	char filedir2[]="C:\\Program Files\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
    char filedir3[]="C:\\Program Files (x86)\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
	//char filedir2[]="C:\\Projects\\edpf\\MGPCMatrix";//������DPU�ļ�Ŀ¼
	AlgIdType AlgId,SheetNo,m_RecId,FinalId;
	AlgId=GetId();//��ȡ�㷨Id
	SheetNo=AlgId>>16;
	m_RecId=AlgId&0x000000FF;
	SheetNo=AlgId>>16;

	short int	St;
	char	filename[MAX_PATH]="";
	char	FilePath[MAX_PATH]="";
	char	StrBuf[256];
	size_t	Len=MAX_PATH;
	short int	rc;
	short int	j;

	FILE *fp=NULL;
	char str[MAX_PATH];
	char *p;
	char *delims={" []=\n"};
	// �����ļ���
	sprintf(filename,"Mgpcv3p1Model%03d%03d.para",SheetNo,m_RecId);
	sprintf(FilePath,"%s/%s",filedir,filename);
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(FilePath,"%s\\%s",filedir2,filename);
		if((fp=fopen(FilePath,"r"))==NULL)
		{
			sprintf(FilePath,"%s\\%s",filedir3,filename);
		    if((fp=fopen(FilePath,"r"))==NULL)
		    {
			    sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
			    printf("%s",StrBuf);
			    AlgState=AlgState|0x0001;
			    return false;
		    }
		}
	}

	AlgState=AlgState&0xFFFE;

	for(int i=0;i<7;i++)
		c_v3p1_ModelPara.flag[i]=0;
	
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
		printf("%s",StrBuf);
		AlgState=AlgState|0x0001;
		return false;
	}

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"totalmodel=",11)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v3p1_ModelPara.totalmodel=strtod(p,NULL)-1;
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v3p1_ModelPara.totalmodel>(M3MIMOGPC_MODELNUM_MAX-1)||c_v3p1_ModelPara.totalmodel<0){
		c_v3p1_ModelPara.flag[0]=0;
		rewind(fp);
	}
	else{
		c_v3p1_ModelPara.flag[0]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"na=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v3p1_ModelPara.na=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v3p1_ModelPara.na>M3MIMOGPC_MODEL_NA_MAX||c_v3p1_ModelPara.na<2)
	{
		c_v3p1_ModelPara.flag[1]=0;
		rewind(fp);
	}
	else
	{
		c_v3p1_ModelPara.flag[1]=1;
	}	
	
	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v3p1_ModelPara.nb=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v3p1_ModelPara.nb>M3MIMOGPC_MODEL_NB_MAX||c_v3p1_ModelPara.nb<1||c_v3p1_ModelPara.nb>c_v3p1_ModelPara.na)
	{
		c_v3p1_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_v3p1_ModelPara.flag[2]=1;
	}

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb2=",4)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v3p1_ModelPara.nb2=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v3p1_ModelPara.nb2>M3MIMOGPC_MODEL_NB_MAX||c_v3p1_ModelPara.nb2>c_v3p1_ModelPara.na)
	{
		c_v3p1_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_v3p1_ModelPara.flag[2]=1;
	}

	AlgState=AlgState&0xFFFE;

	if (c_v3p1_ModelPara.flag[0]&&c_v3p1_ModelPara.flag[1]&&c_v3p1_ModelPara.flag[2])//must read "totalmodel" and "na" before
	{
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(str[0]=='A')
			{
				c_v3p1_ModelPara.flag[4]=1;
				for(int i=0;i<=c_v3p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_v3p1_ModelPara.na))
						  {
							c_v3p1_ModelPara.a[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_v3p1_ModelPara.na)//check a[i][0]=1;
							c_v3p1_ModelPara.flag[4]=0;
					}
					else
					{
						c_v3p1_ModelPara.flag[4]=0;
					}

				}
			}
			if(strncmp(str,"B=",2)==0)
			{
				c_v3p1_ModelPara.flag[5]=1;
				for(int i=0;i<=c_v3p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_v3p1_ModelPara.nb))
						  {
							c_v3p1_ModelPara.b[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_v3p1_ModelPara.nb)
							c_v3p1_ModelPara.flag[5]=0;
					}
					else
					{
						c_v3p1_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"B2=",3)==0)
			{
				c_v3p1_ModelPara.flag[5]=1;
				for(int i=0;i<=c_v3p1_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<c_v3p1_ModelPara.nb2))
						{
							c_v3p1_ModelPara.b2[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j!=c_v3p1_ModelPara.nb2)//check a[i][0]=1;
							c_v3p1_ModelPara.flag[5]=0;
					}
					else
					{
						c_v3p1_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"D=",2)==0)
			{
				c_v3p1_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<=c_v3p1_ModelPara.totalmodel))
						{
							c_v3p1_ModelPara.d[j]=max(m_TunePara.Ts,strtod(p,NULL));
							if (c_v3p1_ModelPara.d[j]<1)
							{
								c_v3p1_ModelPara.d[j]=1;
							//	c3_ModelPara.flag[6]=0;
							}
							if(c_v3p1_ModelPara.d_max<c_v3p1_ModelPara.d[j])
								c_v3p1_ModelPara.d_max=c_v3p1_ModelPara.d[j];
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j<c_v3p1_ModelPara.totalmodel)
							c_v3p1_ModelPara.flag[6]=0;
						//break;
				}
				else
				{
					c_v3p1_ModelPara.flag[6]=0;
				}
			}
			if(strncmp(str,"D2=",3)==0)
			{
				c_v3p1_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
					p=strtok(str,delims);
					j=0;
					while((p!=NULL)&&(j<=c_v3p1_ModelPara.totalmodel))
					{
						c_v3p1_ModelPara.d2[j]=max(m_TunePara.Ts,strtod(p,NULL));
						if (c_v3p1_ModelPara.d2[j]<1)
						{
							c_v3p1_ModelPara.d2[j]=1;
							//	c_ModelPara.flag[6]=0;
						}
						if(c_v3p1_ModelPara.d2_max<c_v3p1_ModelPara.d2[j])
							c_v3p1_ModelPara.d2_max=c_v3p1_ModelPara.d2[j];
						j++;
						p=strtok(NULL,delims); 
					} 
					if (j<c_v3p1_ModelPara.totalmodel)
						c_v3p1_ModelPara.flag[6]=0;
					break;
				}
				else
				{
					c_v3p1_ModelPara.flag[6]=0;
				}
			}
			//fputs(str,stdout);
		}

		//��NLTD����
		rc=1;
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(strncmp(str,"NLTDFLAG=",9)==0)
			{
				p=strtok(str,delims);
				if(p!=NULL)
				{
					p=strtok(NULL,delims);
					c_v3p1_ModelPara.NLTDflag=strtod(p,NULL);
					rc=0;
					break;
				}
			}
		}
		if (rc||c_v3p1_ModelPara.NLTDflag>1||c_v3p1_ModelPara.NLTDflag<0)
		{
			c_v3p1_ModelPara.flag[0]=0;
			rewind(fp);
		}
		else
		{
			c_v3p1_ModelPara.flag[0]=1;
		}
		if(c_v3p1_ModelPara.NLTDflag==1)
		{
			rc=1;
			while(fgets(str,sizeof(str),fp)!=NULL)
			{
				if(strncmp(str,"NLTDna=",7)==0)
				{
					p=strtok(str,delims);
					if(p!=NULL)
					{
						p=strtok(NULL,delims);
						c_v3p1_ModelPara.NLTDna=strtod(p,NULL);
						rc=0;
						break;
					}
				}
			}
			if (rc||c_v3p1_ModelPara.NLTDna>M3MIMOGPC_MODEL_NA_MAX||c_v3p1_ModelPara.NLTDna<2)
			{
				c_v3p1_ModelPara.flag[1]=0;
				rewind(fp);
			}
			else
			{
				c_v3p1_ModelPara.flag[1]=1;
			}	
	
			rc=1;
			while(fgets(str,sizeof(str),fp)!=NULL)
			{
				if(strncmp(str,"NLTDnb=",7)==0)
				{
					p=strtok(str,delims);
					if(p!=NULL)
					{
						p=strtok(NULL,delims);
						c_v3p1_ModelPara.NLTDnb=strtod(p,NULL);
						rc=0;
						break;
					}
				}
			}
			if (rc||c_v3p1_ModelPara.NLTDnb>M3MIMOGPC_MODEL_NB_MAX||c_v3p1_ModelPara.NLTDnb<1||c_v3p1_ModelPara.NLTDnb>c_v3p1_ModelPara.NLTDna)
			{
				c_v3p1_ModelPara.flag[2]=0;
				rewind(fp);
			}
			else
			{
				c_v3p1_ModelPara.flag[2]=1;
			}	
			if (c_v3p1_ModelPara.flag[0]&&c_v3p1_ModelPara.flag[1]&&c_v3p1_ModelPara.flag[2])
			{
				while(fgets(str,sizeof(str),fp)!=NULL)
				{
					if(strncmp(str,"NLTDA=",6)==0)
					{
						c_v3p1_ModelPara.flag[4]=1;
						for(int i=0;i<=c_v3p1_ModelPara.totalmodel;i++)
						{
							if(fgets(str,sizeof(str),fp)!=NULL)
							{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<c_v3p1_ModelPara.NLTDna))
								{
									c_v3p1_ModelPara.NLTD_a[i][j]=strtod(p,NULL);
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j!=c_v3p1_ModelPara.NLTDna)//check a[i][0]=1;
									c_v3p1_ModelPara.flag[4]=0;
							}
							else
							{
								c_v3p1_ModelPara.flag[4]=0;
							}

						}
					}
					if(strncmp(str,"NLTDB=",6)==0)
					{
						c_v3p1_ModelPara.flag[5]=1;
						for(int i=0;i<=c_v3p1_ModelPara.totalmodel;i++)
						{
							if(fgets(str,sizeof(str),fp)!=NULL)
							{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<c_v3p1_ModelPara.NLTDnb))
								{
									c_v3p1_ModelPara.NLTD_b[i][j]=strtod(p,NULL);
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j!=c_v3p1_ModelPara.NLTDnb)
									c_v3p1_ModelPara.flag[5]=0;
							}
							else
							{
								c_v3p1_ModelPara.flag[5]=0;
							}
						}
					}
					if(strncmp(str,"NLTDD=",6)==0)
					{
						c_v3p1_ModelPara.flag[6]=1;
						if(fgets(str,sizeof(str),fp)!=NULL)
						{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<=c_v3p1_ModelPara.totalmodel))
								{
									c_v3p1_ModelPara.NLTD_d[j]=max(m_TunePara.Ts,strtod(p,NULL));
									if (c_v3p1_ModelPara.NLTD_d[j]<1)
									{
										c_v3p1_ModelPara.NLTD_d[j]=1;
										//c3_ModelPara.flag[6]=0;
									}
									if(c_v3p1_ModelPara.NLTD_d_max<c_v3p1_ModelPara.NLTD_d[j])
										c_v3p1_ModelPara.NLTD_d_max=c_v3p1_ModelPara.NLTD_d[j];
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j<c_v3p1_ModelPara.totalmodel)
									c_v3p1_ModelPara.flag[6]=0;
						}
						else
						{
							c_v3p1_ModelPara.flag[6]=0;
						}
					}
				}
			}
			else
			{
				fclose(fp);
				sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
				printf("%s",StrBuf);
				return false;
			}
		
		}
		

		fclose(fp);
		sprintf(StrBuf,"Para File: %s Load OK.\n",FilePath);
		printf("%s",StrBuf);
		c_v3p1_ModelPara.T=m_TunePara.Ts;
		c_v3p1_ModelPara.flag[3]=1;//T���������ļ�����
		j=0;//AlgState=AlgState&0xFFFE;AlgState=AlgState|0x0001;		
		for (int i=0;i<7;i++){
			m_v3p1_ModelPara.flag[i]=c_v3p1_ModelPara.flag[i];
			if (1==c_v3p1_ModelPara.flag[i])
			{
				j++;
				AlgState=AlgState&(~(0x0001<<(i+1)));
			}
			else
				AlgState=AlgState|0x0001<<(i+1);
		}

		if (7==j)
		{
			//����ϵͳ��ɢ��na>nb
			m_v3p1_ModelPara.totalmodel=c_v3p1_ModelPara.totalmodel;
			m_v3p1_ModelPara.na=c_v3p1_ModelPara.na;
			m_v3p1_ModelPara.nb=c_v3p1_ModelPara.na;
			m_v3p1_ModelPara.nb2=c_v3p1_ModelPara.na;
			m_v3p1_ModelPara.T=c_v3p1_ModelPara.T*1000;//����Ϊ��λ
			m_v3p1_ModelPara.d_max=1;
			m_v3p1_ModelPara.d2_max=1;

			m_v3p1_ModelPara.NLTDna=c_v3p1_ModelPara.NLTDna;
			m_v3p1_ModelPara.NLTDnb=c_v3p1_ModelPara.NLTDna;
			m_v3p1_ModelPara.NLTD_d_max=1;

			for(int i=0;i<=m_v3p1_ModelPara.totalmodel;i++)
			{
				if(c_v3p1_ModelPara.nb2!=0)
				{
					myc2d2(c_v3p1_ModelPara.a[i],c_v3p1_ModelPara.b[i],c_v3p1_ModelPara.b2[i],c_v3p1_ModelPara.na,c_v3p1_ModelPara.nb,c_v3p1_ModelPara.nb2,m_v3p1_ModelPara.T/1000,m_v3p1_ModelPara.a[i],m_v3p1_ModelPara.b[i],m_v3p1_ModelPara.b2[i]);
				}
				else
				{
					myc2d(c_v3p1_ModelPara.a[i],c_v3p1_ModelPara.b[i],c_v3p1_ModelPara.na,c_v3p1_ModelPara.nb,m_v3p1_ModelPara.T/1000,m_v3p1_ModelPara.a[i],m_v3p1_ModelPara.b[i]);
				}
				m_v3p1_ModelPara.d[i]=c_v3p1_ModelPara.d[i]/(m_v3p1_ModelPara.T/1000);
				if(c_v3p1_ModelPara.nb2!=0)
				{
					m_v3p1_ModelPara.d2[i]=c_v3p1_ModelPara.d2[i]/(m_v3p1_ModelPara.T/1000);
				}
				if(m_v3p1_ModelPara.d_max<m_v3p1_ModelPara.d[i])
					m_v3p1_ModelPara.d_max=m_v3p1_ModelPara.d[i];
				if(c_v3p1_ModelPara.nb2!=0)
				{
					if(m_v3p1_ModelPara.d2_max<m_v3p1_ModelPara.d2[i])
						m_v3p1_ModelPara.d2_max=m_v3p1_ModelPara.d2[i];
					m_v3p1_ModelPara.d_max = max(m_v3p1_ModelPara.d_max,m_v3p1_ModelPara.d2_max);
				}
				double ta[2]={1,-1};
				conv(m_v3p1_ModelPara.a[i],ta,m_v3p1_ModelPara.na,2,m_v3p1_ModelPara.aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
				if(c_v3p1_ModelPara.NLTDflag==1)
				{
					myc2d(c_v3p1_ModelPara.NLTD_a[i],c_v3p1_ModelPara.NLTD_b[i],c_v3p1_ModelPara.NLTDna,c_v3p1_ModelPara.NLTDnb,m_v3p1_ModelPara.T/1000,m_v3p1_ModelPara.NLTD_a[i],m_v3p1_ModelPara.NLTD_b[i]);
					m_v3p1_ModelPara.NLTD_d[i]=c_v3p1_ModelPara.NLTD_d[i]/(m_v3p1_ModelPara.T/1000);
					if(m_v3p1_ModelPara.NLTD_d_max<m_v3p1_ModelPara.NLTD_d[i])
						m_v3p1_ModelPara.NLTD_d_max=m_v3p1_ModelPara.NLTD_d[i];
					conv(m_v3p1_ModelPara.NLTD_a[i],ta,m_v3p1_ModelPara.NLTDna,2,m_v3p1_ModelPara.NLTD_aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};	
				}
				
			}
				return true;
		}
		else{
			//cout<< "input wrong!" << endl;
			return false;		
		}
	}
	else
	{
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
		printf("%s",StrBuf);
		return false;
	}
}

/*==============================================������========================================================*/
//�������ļ���ȡϵͳģ��
//�ļ���Ϊ��MgpcModel+���+.para;���Ϊ6λ���֣�ǰ��λΪҳ�ţ�������λ��0������λΪID�ţ�������λ��0;
//���磺��5ҳ���㷨IDΪ3�����Ϊ��005003�������ļ���Ϊ��MgpcModel005003.para;
//�ļ�Ŀ¼������Ŀ¼\domain013\data\ctrl\obj\drop001
bool AlgICSMIMOGPC::ReadParav3p2FromFile()//Ҫ�� MR ?
{	
	char filedir[]="\/edpf\/prj\/running\/ctrl\/obj";//���ļ�Ŀ¼
    //char filedir2[]="E:\\project\\DPU59running\\ctrl\\obj";//������DPU�ļ�Ŀ¼
	char filedir2[]="C:\\Program Files\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
    char filedir3[]="C:\\Program Files (x86)\\EDPF-NT plus\\MGPCMatrix\\";//������DPU�ļ�Ŀ¼
	//char filedir2[]="C:\\Projects\\edpf\\MGPCMatrix";//������DPU�ļ�Ŀ¼
	AlgIdType AlgId,SheetNo,m_RecId,FinalId;
	AlgId=GetId();//��ȡ�㷨Id
	SheetNo=AlgId>>16;
	m_RecId=AlgId&0x000000FF;
	SheetNo=AlgId>>16;

	short int	St;
	char	filename[MAX_PATH]="";
	char	FilePath[MAX_PATH]="";
	char	StrBuf[256];
	size_t	Len=MAX_PATH;
	short int	rc;
	short int	j;

	FILE *fp=NULL;
	char str[MAX_PATH];
	char *p;
	char *delims={" []=\n"};
	// �����ļ���
	sprintf(filename,"Mgpcv3p2Model%03d%03d.para",SheetNo,m_RecId);
	sprintf(FilePath,"%s/%s",filedir,filename);
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(FilePath,"%s\\%s",filedir2,filename);
		if((fp=fopen(FilePath,"r"))==NULL)
		{
			sprintf(FilePath,"%s\\%s",filedir3,filename);
		    if((fp=fopen(FilePath,"r"))==NULL)
		    {
			    sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
			    printf("%s",StrBuf);
			    AlgState=AlgState|0x0001;
			    return false;
		    }
		}
	}

	AlgState=AlgState&0xFFFE;

	for(int i=0;i<7;i++)
		c_v3p2_ModelPara.flag[i]=0;
	
	if((fp=fopen(FilePath,"r"))==NULL)
	{
		sprintf(StrBuf,"Can not open Para File: %s.\n",FilePath);
		printf("%s",StrBuf);
		AlgState=AlgState|0x0001;
		return false;
	}

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"totalmodel=",11)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v3p2_ModelPara.totalmodel=strtod(p,NULL)-1;
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v3p2_ModelPara.totalmodel>(M3MIMOGPC_MODELNUM_MAX-1)||c_v3p2_ModelPara.totalmodel<0){
		c_v3p2_ModelPara.flag[0]=0;
		rewind(fp);
	}
	else{
		c_v3p2_ModelPara.flag[0]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"na=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v3p2_ModelPara.na=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v3p2_ModelPara.na>M3MIMOGPC_MODEL_NA_MAX||c_v3p2_ModelPara.na<2)
	{
		c_v3p2_ModelPara.flag[1]=0;
		rewind(fp);
	}
	else
	{
		c_v3p2_ModelPara.flag[1]=1;
	}	
	
	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb=",2)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v3p2_ModelPara.nb=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v3p2_ModelPara.nb>M3MIMOGPC_MODEL_NB_MAX||c_v3p2_ModelPara.nb<1||c_v3p2_ModelPara.nb>c_v3p2_ModelPara.na)
	{
		c_v3p2_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_v3p2_ModelPara.flag[2]=1;
	}	

	rc=1;
	while(fgets(str,sizeof(str),fp)!=NULL)
	{
		if(strncmp(str,"nb2=",4)==0)
		{
			 p=strtok(str,delims);
			 if(p!=NULL)
			 {
				p=strtok(NULL,delims);
				c_v3p2_ModelPara.nb2=strtod(p,NULL);
				rc=0;
				break;
			 }
		}
	}
	if (rc||c_v3p2_ModelPara.nb2>M3MIMOGPC_MODEL_NB_MAX||c_v3p2_ModelPara.nb2>c_v3p2_ModelPara.na)
	{
		c_v3p2_ModelPara.flag[2]=0;
		rewind(fp);
	}
	else
	{
		c_v3p2_ModelPara.flag[2]=1;
	}

	AlgState=AlgState&0xFFFE;

	if (c_v3p2_ModelPara.flag[0]&&c_v3p2_ModelPara.flag[1]&&c_v3p2_ModelPara.flag[2])//must read "totalmodel" and "na" before
	{
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(str[0]=='A')
			{
				c_v3p2_ModelPara.flag[4]=1;
				for(int i=0;i<=c_v3p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_v3p2_ModelPara.na))
						  {
							c_v3p2_ModelPara.a[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_v3p2_ModelPara.na)//check a[i][0]=1;
							c_v3p2_ModelPara.flag[4]=0;
					}
					else
					{
						c_v3p2_ModelPara.flag[4]=0;
					}

				}
			}
			if(strncmp(str,"B=",2)==0)
			{
				c_v3p2_ModelPara.flag[5]=1;
				for(int i=0;i<=c_v3p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						 p=strtok(str,delims);
						 j=0;
						  while((p!=NULL)&&(j<c_v3p2_ModelPara.nb))
						  {
							c_v3p2_ModelPara.b[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						  } 
						  if (j!=c_v3p2_ModelPara.nb)
							c_v3p2_ModelPara.flag[5]=0;
					}
					else
					{
						c_v3p2_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"B2=",3)==0)
			{
				c_v3p2_ModelPara.flag[5]=1;
				for(int i=0;i<=c_v3p2_ModelPara.totalmodel;i++)
				{
					if(fgets(str,sizeof(str),fp)!=NULL)
					{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<c_v3p2_ModelPara.nb2))
						{
							c_v3p2_ModelPara.b2[i][j]=strtod(p,NULL);
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j!=c_v3p2_ModelPara.nb2)//check a[i][0]=1;
							c_v3p2_ModelPara.flag[5]=0;
					}
					else
					{
						c_v3p2_ModelPara.flag[5]=0;
					}
				}
			}
			if(strncmp(str,"D=",2)==0)
			{
				c_v3p2_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
						p=strtok(str,delims);
						j=0;
						while((p!=NULL)&&(j<=c_v3p2_ModelPara.totalmodel))
						{
							c_v3p2_ModelPara.d[j]=max(m_TunePara.Ts,strtod(p,NULL));
							if (c_v3p2_ModelPara.d[j]<1)
							{
								c_v3p2_ModelPara.d[j]=1;
							//	c3_ModelPara.flag[6]=0;
							}
							if(c_v3p2_ModelPara.d_max<c_v3p2_ModelPara.d[j])
								c_v3p2_ModelPara.d_max=c_v3p2_ModelPara.d[j];
							j++;
							p=strtok(NULL,delims); 
						} 
						if (j<c_v3p2_ModelPara.totalmodel)
							c_v3p2_ModelPara.flag[6]=0;
						//break;
				}
				else
				{
					c_v3p2_ModelPara.flag[6]=0;
				}
			}
			if(strncmp(str,"D2=",3)==0)
			{
				c_v3p2_ModelPara.flag[6]=1;
				if(fgets(str,sizeof(str),fp)!=NULL)
				{
					p=strtok(str,delims);
					j=0;
					while((p!=NULL)&&(j<=c_v3p2_ModelPara.totalmodel))
					{
						c_v3p2_ModelPara.d2[j]=max(m_TunePara.Ts,strtod(p,NULL));
						if (c_v3p2_ModelPara.d2[j]<1)
						{
							c_v3p2_ModelPara.d2[j]=1;
							//	c_ModelPara.flag[6]=0;
						}
						if(c_v3p2_ModelPara.d2_max<c_v3p2_ModelPara.d2[j])
							c_v3p2_ModelPara.d2_max=c_v3p2_ModelPara.d2[j];
						j++;
						p=strtok(NULL,delims); 
					} 
					if (j<c_v3p2_ModelPara.totalmodel)
						c_v3p2_ModelPara.flag[6]=0;
					break;
				}
				else
				{
					c_v3p2_ModelPara.flag[6]=0;
				}
			}

			//fputs(str,stdout);
		}

		//��NLTD����
		rc=1;
		while(fgets(str,sizeof(str),fp)!=NULL)
		{
			if(strncmp(str,"NLTDFLAG=",9)==0)
			{
				p=strtok(str,delims);
				if(p!=NULL)
				{
					p=strtok(NULL,delims);
					c_v3p2_ModelPara.NLTDflag=strtod(p,NULL);
					rc=0;
					break;
				}
			}
		}
		if (rc||c_v3p2_ModelPara.NLTDflag>1||c_v3p2_ModelPara.NLTDflag<0)
		{
			c_v3p2_ModelPara.flag[0]=0;
			rewind(fp);
		}
		else
		{
			c_v3p2_ModelPara.flag[0]=1;
		}
		if(c_v3p2_ModelPara.NLTDflag==1)
		{
			rc=1;
			while(fgets(str,sizeof(str),fp)!=NULL)
			{
				if(strncmp(str,"NLTDna=",7)==0)
				{
					p=strtok(str,delims);
					if(p!=NULL)
					{
						p=strtok(NULL,delims);
						c_v3p2_ModelPara.NLTDna=strtod(p,NULL);
						rc=0;
						break;
					}
				}
			}
			if (rc||c_v3p2_ModelPara.NLTDna>M3MIMOGPC_MODEL_NA_MAX||c_v3p2_ModelPara.NLTDna<2)
			{
				c_v3p2_ModelPara.flag[1]=0;
				rewind(fp);
			}
			else
			{
				c_v3p2_ModelPara.flag[1]=1;
			}	
	
			rc=1;
			while(fgets(str,sizeof(str),fp)!=NULL)
			{
				if(strncmp(str,"NLTDnb=",7)==0)
				{
					p=strtok(str,delims);
					if(p!=NULL)
					{
						p=strtok(NULL,delims);
						c_v3p2_ModelPara.NLTDnb=strtod(p,NULL);
						rc=0;
						break;
					}
				}
			}
			if (rc||c_v3p2_ModelPara.NLTDnb>M3MIMOGPC_MODEL_NB_MAX||c_v3p2_ModelPara.NLTDnb<1||c_v3p2_ModelPara.NLTDnb>c_v3p2_ModelPara.NLTDna)
			{
				c_v3p2_ModelPara.flag[2]=0;
				rewind(fp);
			}
			else
			{
				c_v3p2_ModelPara.flag[2]=1;
			}	
			if (c_v3p2_ModelPara.flag[0]&&c_v3p2_ModelPara.flag[1]&&c_v3p2_ModelPara.flag[2])
			{
				while(fgets(str,sizeof(str),fp)!=NULL)
				{
					if(strncmp(str,"NLTDA=",6)==0)
					{
						c_v3p2_ModelPara.flag[4]=1;
						for(int i=0;i<=c_v3p2_ModelPara.totalmodel;i++)
						{
							if(fgets(str,sizeof(str),fp)!=NULL)
							{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<c_v3p2_ModelPara.NLTDna))
								{
									c_v3p2_ModelPara.NLTD_a[i][j]=strtod(p,NULL);
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j!=c_v3p2_ModelPara.NLTDna)//check a[i][0]=1;
									c_v3p2_ModelPara.flag[4]=0;
							}
							else
							{
								c_v3p2_ModelPara.flag[4]=0;
							}

						}
					}
					if(strncmp(str,"NLTDB=",6)==0)
					{
						c_v3p2_ModelPara.flag[5]=1;
						for(int i=0;i<=c_v3p2_ModelPara.totalmodel;i++)
						{
							if(fgets(str,sizeof(str),fp)!=NULL)
							{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<c_v3p2_ModelPara.NLTDnb))
								{
									c_v3p2_ModelPara.NLTD_b[i][j]=strtod(p,NULL);
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j!=c_v3p2_ModelPara.NLTDnb)
									c_v3p2_ModelPara.flag[5]=0;
							}
							else
							{
								c_v3p2_ModelPara.flag[5]=0;
							}
						}
					}
					if(strncmp(str,"NLTDD=",6)==0)
					{
						c_v3p2_ModelPara.flag[6]=1;
						if(fgets(str,sizeof(str),fp)!=NULL)
						{
								p=strtok(str,delims);
								j=0;
								while((p!=NULL)&&(j<=c_v3p2_ModelPara.totalmodel))
								{
									c_v3p2_ModelPara.NLTD_d[j]=max(m_TunePara.Ts,strtod(p,NULL));
									if (c_v3p2_ModelPara.NLTD_d[j]<1)
									{
										c_v3p2_ModelPara.NLTD_d[j]=1;
										//c3_ModelPara.flag[6]=0;
									}
									if(c_v3p2_ModelPara.NLTD_d_max<c_v3p2_ModelPara.NLTD_d[j])
										c_v3p2_ModelPara.NLTD_d_max=c_v3p2_ModelPara.NLTD_d[j];
									j++;
									p=strtok(NULL,delims); 
								} 
								if (j<c_v3p2_ModelPara.totalmodel)
									c_v3p2_ModelPara.flag[6]=0;
						}
						else
						{
							c_v3p2_ModelPara.flag[6]=0;
						}
					}
				}
			}
			else
			{
				fclose(fp);
				sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
				printf("%s",StrBuf);
				return false;
			}
		
		}
		

		fclose(fp);
		sprintf(StrBuf,"Para File: %s Load OK.\n",FilePath);
		printf("%s",StrBuf);
		c_v3p2_ModelPara.T=m_TunePara.Ts;
		c_v3p2_ModelPara.flag[3]=1;//T���������ļ�����
		j=0;//AlgState=AlgState&0xFFFE;AlgState=AlgState|0x0001;		
		for (int i=0;i<7;i++){
			m_v3p2_ModelPara.flag[i]=c_v3p2_ModelPara.flag[i];
			if (1==c_v3p2_ModelPara.flag[i])
			{
				j++;
				AlgState=AlgState&(~(0x0001<<(i+1)));
			}
			else
				AlgState=AlgState|0x0001<<(i+1);
		}

		if (7==j)
		{
			//����ϵͳ��ɢ��na>nb
			m_v3p2_ModelPara.totalmodel=c_v3p2_ModelPara.totalmodel;
			m_v3p2_ModelPara.na=c_v3p2_ModelPara.na;
			m_v3p2_ModelPara.nb=c_v3p2_ModelPara.na;
			m_v3p2_ModelPara.nb2=c_v3p2_ModelPara.na;
			m_v3p2_ModelPara.T=c_v3p2_ModelPara.T*1000;//����Ϊ��λ
			m_v3p2_ModelPara.d_max=1;
			m_v3p2_ModelPara.d2_max=1;

			m_v3p2_ModelPara.NLTDna=c_v3p2_ModelPara.NLTDna;
			m_v3p2_ModelPara.NLTDnb=c_v3p2_ModelPara.NLTDna;
			m_v3p2_ModelPara.NLTD_d_max=1;

			for(int i=0;i<=m_v3p2_ModelPara.totalmodel;i++)
			{
				if(c_v3p2_ModelPara.nb2!=0)
				{
					myc2d2(c_v3p2_ModelPara.a[i],c_v3p2_ModelPara.b[i],c_v3p2_ModelPara.b2[i],c_v3p2_ModelPara.na,c_v3p2_ModelPara.nb,c_v3p2_ModelPara.nb2,m_v3p2_ModelPara.T/1000,m_v3p2_ModelPara.a[i],m_v3p2_ModelPara.b[i],m_v3p2_ModelPara.b2[i]);
				}
				else
				{
					myc2d(c_v3p2_ModelPara.a[i],c_v3p2_ModelPara.b[i],c_v3p2_ModelPara.na,c_v3p2_ModelPara.nb,m_v3p2_ModelPara.T/1000,m_v3p2_ModelPara.a[i],m_v3p2_ModelPara.b[i]);
				}
				m_v3p2_ModelPara.d[i]=c_v3p2_ModelPara.d[i]/(m_v3p2_ModelPara.T/1000);
				if(c_v3p2_ModelPara.nb2!=0)
				{
					m_v3p2_ModelPara.d2[i]=c_v3p2_ModelPara.d2[i]/(m_v3p2_ModelPara.T/1000);
				}
				if(m_v3p2_ModelPara.d_max<m_v3p2_ModelPara.d[i])
					m_v3p2_ModelPara.d_max=m_v3p2_ModelPara.d[i];
				if(c_v3p2_ModelPara.nb2!=0)
				{
					if(m_v3p2_ModelPara.d2_max<m_v3p2_ModelPara.d2[i])
						m_v3p2_ModelPara.d2_max=m_v3p2_ModelPara.d2[i];
					m_v3p2_ModelPara.d_max = max(m_v3p2_ModelPara.d_max,m_v3p2_ModelPara.d2_max);
				}
				double ta[2]={1,-1};
				conv(m_v3p2_ModelPara.a[i],ta,m_v3p2_ModelPara.na,2,m_v3p2_ModelPara.aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};
				if(c_v3p2_ModelPara.NLTDflag==1)
				{
					myc2d(c_v3p2_ModelPara.NLTD_a[i],c_v3p2_ModelPara.NLTD_b[i],c_v3p2_ModelPara.NLTDna,c_v3p2_ModelPara.NLTDnb,m_v3p2_ModelPara.T/1000,m_v3p2_ModelPara.NLTD_a[i],m_v3p2_ModelPara.NLTD_b[i]);
					m_v3p2_ModelPara.NLTD_d[i]=c_v3p2_ModelPara.NLTD_d[i]/(m_v3p2_ModelPara.T/1000);
					if(m_v3p2_ModelPara.NLTD_d_max<m_v3p2_ModelPara.NLTD_d[i])
						m_v3p2_ModelPara.NLTD_d_max=m_v3p2_ModelPara.NLTD_d[i];
					conv(m_v3p2_ModelPara.NLTD_a[i],ta,m_v3p2_ModelPara.NLTDna,2,m_v3p2_ModelPara.NLTD_aa[i]);//aa=conv(a,[1 -1]);double ta[2]={1,-1};	
				}
				
			}
				return true;
		}
		else{
			//cout<< "input wrong!" << endl;
			return false;		
		}
	}
	else
	{
		fclose(fp);
		sprintf(StrBuf,"Para File: %s Format Error.\n",FilePath);
		printf("%s",StrBuf);
		return false;
	}
}



//����ϵͳ��ɢ�����򣬲���˫���Ա任����tustin����
//s=(2/Ts)*(1-z^-1)/(1+z^-1),
//aΪ����ϵͳ��ĸ����ʽ,bΪ����ϵͳ���Ӷ���ʽ
//d_aΪ��ɢϵͳ��ĸ����ʽ,d_bΪ��ɢϵͳ���Ӷ���ʽ,���鳤��Ϊna
void AlgICSMIMOGPC::myc2d(double* a,double* b,int na,int nb,double Ts,double* d_a,double* d_b)
{
	int i,j;
	double p1[2]={1,-1};
	double p2[2]={1,1};
	for (i=1;i<na;i++)
	{
		a[i-1]=pow(2/Ts,na-i)*a[i-1];
	}
	if (nb>1)
	{
		for (i=1;i<nb;i++)
		{
			b[i-1]=pow(2/Ts,nb-i)*b[i-1];
		}
	}
	double** A,**B;
	malloc_nm(na,na,sizeof(double),sizeof(double*),(double***)&A);
	malloc_nm(na,na,sizeof(double),sizeof(double*),(double***)&B);
	for(i=1;i<=na;i++)
	{
		A[i-1][0]=a[i-1];
	}
	for(i=1;i<=nb;i++)
	{
		B[na-nb+i-1][0]=b[i-1];
	}

	for(i=1;i<=na;i++)
	{
		for (j=1;j<=na-i;j++)
		{
			double* tmp1;
			tmp1=(double*)malloc(j*sizeof(double));
			for(int j1=0;j1<j;j1++)
			{
				tmp1[j1]=A[i-1][j1];
			}
			conv(tmp1,p1,j,2,A[i-1]);
			for(int j1=0;j1<j;j1++)
			{
				tmp1[j1]=B[i-1][j1];
			}
			conv(tmp1,p1,j,2,B[i-1]);
			free(tmp1);
		}
		if(i>1)
		{
			for(j=1;j<=i-1;j++)
			{
				double *tmp1;
				tmp1=(double*)malloc((na-i+j)*sizeof(double));
				for(int j1=0;j1<na-i+j;j1++)
				{
					tmp1[j1]=A[i-1][j1];
				}
				conv(tmp1,p2,na-i+j,2,A[i-1]);
				for(int j1=0;j1<na-i+j;j1++)
				{
					tmp1[j1]=B[i-1][j1];
				}
				conv(tmp1,p2,na-i+j,2,B[i-1]);
				free(tmp1);
			}
		}
	}
	
	for (i=2;i<=na;i++)
	{
		for(j=1;j<=na;j++)
		{
			A[0][j-1]=A[0][j-1]+A[i-1][j-1];
			B[0][j-1]=B[0][j-1]+B[i-1][j-1];
		}
	}
	double k=A[0][0];
	if (k!=0)
	{
		for (i=1;i<=na;i++)
		{
			A[0][i-1]=A[0][i-1]/k;
			d_a[i-1]=A[0][i-1];
			B[0][i-1]=B[0][i-1]/k;
			d_b[i-1]=B[0][i-1];
		}
	}

	
/* 	cout << "Discretization completed !" << endl;
	cout << "Den=" << endl;
	for(i=0;i<na;i++)
	{
		cout << A[0][i] << " ";
	}
	cout << endl << "Num=" << endl;
	for(i=0;i<na;i++)
	{
		cout << B[0][i] << " ";
	}
	cout << endl; */
	
	free_nm(na,(double**)A);
	free_nm(na,(double**)B);
	
}

//LUP�ֽ�
void AlgICSMIMOGPC::LUP_Descomposition(double *A,double *L,double *U,int *P,int N)
{
    int row=0;
    for(int i=0;i<N;i++)
    {
        P[i]=i;
    }
    for(int i=0;i<N-1;i++)
    {
        double p=0;
        for(int j=i;j<N;j++)
        {
            if(fabs(A[j*N+i])>p)
            {
                p=fabs(A[j*N+i]);
                row=j;
            }
        }
        if(0==p)
        {
            //cout<< "�������죬�޷�������" <<endl;
            return ;
        }

        //����P[i]��P[row]
        int tmp=P[i];
        P[i]=P[row];
        P[row]=tmp;

        double tmp2=0;
        for(int j=0;j<N;j++)
        {
            //����A[i][j]�� A[row][j]
            tmp2=A[i*N+j];
            A[i*N+j]=A[row*N+j];
            A[row*N+j]=tmp2;
        }

        //����ͬLU�ֽ�
        double u=A[i*N+i],l=0;
        for(int j=i+1;j<N;j++)
        {
            l=A[j*N+i]/u;
            A[j*N+i]=l;
            for(int k=i+1;k<N;k++)
            {
                A[j*N+k]=A[j*N+k]-A[i*N+k]*l;
            }
        }

    }

    //����L��U
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<=i;j++)
        {
            if(i!=j)
            {
                L[i*N+j]=A[i*N+j];
            }
            else
            {
                L[i*N+j]=1;
            }
        }
        for(int k=i;k<N;k++)
        {
            U[i*N+k]=A[i*N+k];
        }
    }

}

//LUP��ⷽ��
void AlgICSMIMOGPC::LUP_Solve(double *L,double *U,int *P,double *b,int N,double *x)
{
    //double *x=new double[N]();
    double *y=new double[N]();

    //�����滻
    for(int i = 0;i < N;i++)
    {
        y[i] = b[P[i]];
        for(int j = 0;j < i;j++)
        {
            y[i] = y[i] - L[i*N+j]*y[j];
        }
    }
    //�����滻
    for(int i = N-1;i >= 0; i--)
    {
        x[i]=y[i];
        for(int j = N-1;j > i;j--)
        {
            x[i] = x[i] - U[i*N+j]*x[j];
        }
        x[i] /= U[i*N+i];
	}
	//MGPC2_Matrix.x = x;
	//delete []x;
	delete []y;
    return;
}

/*****************����ԭ��ת��BEGIN********************/

/* ��� */
int AlgICSMIMOGPC::getNext(int i, int m, int n)
{
  return (i%n)*m + i/n;
}

/* ǰ�� */
int AlgICSMIMOGPC::getPre(int i, int m, int n)
{
  return (i%m)*n + i/m;
}

/* �������±�iΪ���Ļ� */
void AlgICSMIMOGPC::movedata(double *mtx, int i, int m, int n)
{
  double temp = mtx[i]; // �ݴ�
  int cur = i;    // ��ǰ�±�
  int pre = getPre(cur, m, n);
  while(pre != i)
  {
    mtx[cur] = mtx[pre];
    cur = pre;
    pre = getPre(cur, m, n);
  }
  mtx[cur] = temp;
}

/* ת�ã���ѭ���������л� */
void AlgICSMIMOGPC::transpose(double *mtx, int m, int n)
{
  for(int i=0; i<m*n; ++i)
  {
    int next = getNext(i, m, n);
    while(next > i) // �����ں��С��i˵���ظ�,�Ͳ�������ȥ�ˣ�ֻ�в��ظ�ʱ����whileѭ����
      next = getNext(next, m, n);
    if(next == i)  // ����ǰ��
      movedata(mtx, i, m, n);
  }
}
/*****************����ԭ��ת��END********************/

//LUP����(��ÿ��b����ĸ���x������װ)
void AlgICSMIMOGPC::LUP_solve_inverse(double *A , int N , double *inv_A)
{
    //��������A�ĸ�����ע�ⲻ��ֱ����A���㣬��ΪLUP�ֽ��㷨�ѽ���ı�
    double *A_mirror = new double[N*N]();
    //double *inv_A=new double[N*N]();//���յ�����󣨻���Ҫת�ã�
	double *inv_A_each=new double[N]();//������ĸ���
    double *b    =new double[N]();//b��ΪB����о������

    for(int i=0;i<N;i++)
    {
        double *L=new double[N*N]();
        double *U=new double[N*N]();
        int *P=new int[N]();

        //���쵥λ���ÿһ��
        for(int j=0;j<N;j++)
        {
            b[j]=0;
        }
        b[i]=1;

        //ÿ�ζ���Ҫ���½�A����һ��
        for(int j=0;j<N*N;j++)
        {
            A_mirror[j]=A[j];
        }

        LUP_Descomposition(A_mirror,L,U,P,N);

        LUP_Solve (L,U,P,b,N,inv_A_each);
        memcpy(inv_A+i*N,inv_A_each,N*sizeof(double));//������ƴ������
		delete []L;
		delete []U;
		delete []P;
		
	}
    transpose(inv_A,N,N);//�������ڸ���ÿ��b�����x���д洢�������ת��
	//MGPC2_Matrix.inv_A = inv_A;
	delete []A_mirror;
	//delete []inv_A;
	delete []b;
	delete []inv_A_each;
    return;
}

void AlgICSMIMOGPC::CRateCalculate(double **A ,double **B , vector<double>& result1 ,vector<double>& result2 ,int GMA ,int GMB ,int GN,float Gama1,float Gama2,float eta1,float eta2)
{
	double *GA_TEMP = new double[GN*GN]();
	double *GA_TEMP2 = new double[GN*GN]();
	double *GB_TEMP = new double[GN*GN]();
	double *GB_TEMP2 = new double[GN*GN]();
	double *GZ_TEMP = new double[GN*GN]();
	double *GZ_TEMP2 = new double[GN*GN]();
	vector<double> GA_transpose;
	vector<double> GB_transpose;
	GA_transpose.resize(GN*GMA,0);
	GB_transpose.resize(GN*GMB,0);
	
	for(int i=0;i<GN;i++)
    {
        for(int j=0;j<GMA;j++)
		{
			result1[i*GMA+j] = 0 ;
		}
		for(int j=0;j<GMB;j++)
		{
			result2[i*GMB+j] = 0 ;
		}
	}


//GA��ת��
	for(int i=0;i<GN;i++)
    {
        for(int j=0;j<GMA;j++)
        {
            GA_transpose[i*GMA+j] = A[j][i];
        }
    }
//GB��ת��
	for(int i=0;i<GN;i++)
    {
        for(int j=0;j<GMB;j++)
        {
            GB_transpose[i*GMB+j] = B[j][i];
        }
    }
//�Լ���ת����GA���(GA_transpose*GA)
	for(int i=0;i<GN;i++)
    {
        for(int j=0;j<GN;j++)
        {
            for(int k=0;k<GMA;k++)
            {
                GA_TEMP[i*GN+j] += GA_transpose[i*GMA+k]*A[k][j];
            }
        }
    }
//�Լ���ת����GB���(GB_transpose*GB)
	for(int i=0;i<GN;i++)
    {
        for(int j=0;j<GN;j++)
        {
            for(int k=0;k<GMB;k++)
            {
                GB_TEMP[i*GN+j] += GB_transpose[i*GMB+k]*B[k][j];
            }
        }
    }

//���Թ�����Ȩ�ز����(QA*GA_transpose*GA+QB*GB_transpose*GB)
	for(int i=0;i<GN*GN;i++)
	{
		GA_TEMP[i]=GA_TEMP[i]*eta1;
		GB_TEMP[i]=GB_TEMP[i]*eta2;
	}
	for(int i=0;i<GN*GN;i++)
	{
		GZ_TEMP[i]=GA_TEMP[i]+GB_TEMP[i];
	}
//QA*GA_transpose*GA+QB*GB_transpose*GB+V
	for(int i=0;i<m_TunePara.Nu1;i++)
	{
		 GZ_TEMP[i*GN+i]=GZ_TEMP[i*GN+i]+Gama1;
	}
	for(int i=m_TunePara.Nu1;i<GN;i++)
	{
		 GZ_TEMP[i*GN+i]=GZ_TEMP[i*GN+i]+Gama2;
	}
	
//�������
	LUP_solve_inverse( GZ_TEMP,GN,GZ_TEMP2);
/*	
	for(int ia=0;ia<GN;ia++)
	{
		for(int ja=0;ja<GN;ja++)
		{	
			printf("%f ",GZ_TEMP2[ia*GN+ja]);
		}
		printf("\n");
	}
*/

	for(int i=0;i<GN*GN;i++)
	{
		GA_TEMP2[i]=GZ_TEMP2[i]*eta1;
		GB_TEMP2[i]=GZ_TEMP2[i]*eta2;
	}
	//LUP_solve_inverse( GZ_TEMP,GN,GZ_TEMP2);	
//���Լ���ת�����
	for(int i=0;i<GN;i++)
    {
        for(int j=0;j<GMA;j++)
        {
            for(int k=0;k<GN;k++)
            {
                result1[i*GMA+j] += GA_TEMP2[i*GN+k]*GA_transpose[k*GMA+j];
            }
        }
    }
	for(int i=0;i<GN;i++)
    {
        for(int j=0;j<GMB;j++)
        {
            for(int k=0;k<GN;k++)
            {
                result2[i*GMB+j] += GB_TEMP2[i*GN+k]*GB_transpose[k*GMB+j];
            }
        }
    }
	delete []GA_TEMP;
	delete []GA_TEMP2;
	delete []GB_TEMP;
	delete []GB_TEMP2;
	delete []GZ_TEMP;
	delete []GZ_TEMP2;
	//DeletePnt();
	printf("CRateCalculate OK\n");
	return;
}

int AlgICSMIMOGPC::Calculate_max(int a,int b)
{
	if (a>b)
		return a;
	else
		return b;
}
int AlgICSMIMOGPC::Calculate_min(int a,int b)
{
	if (a>b)
		return b;
	else
		return a;
}
/*============================================================================================================
����������G1��G2����ƴ�ӳɾ���GF���ײ�����
���磺
G1={1,0,       G2={1,0,0,       GF={0,0,1,0,0,
    1,1}           1,1,0,           1,0,1,1,0,
	               1,1,1}           1,1,1,1,1}
===============================================================================================================*/
void AlgICSMIMOGPC::MatrixCombine(double** G1,double** G2,int G1N,int G1M ,int G2N,int G2M ,double** GF)
{
	int dev=0;//����������ƫ��

	if(G1N>G2N)
	{
		dev=G1N-G2N;
		for(int i=0;i<G1N;i++)
		{	
			for(int j=0;j<G1M;j++)
			{
				GF[i][j]=G1[i][j];
			}
			if(i>=dev)
			{
				for(int x=G1M;x<G1M+G2M;x++)
				{
					GF[i][x]=G2[i-dev][x-G1M];
				}
			}
		}

	}
	else
	{
		dev=G2N-G1N;
		for(int i=0;i<G2N;i++)
		{	
			
			if(i>=dev)
			{
				for(int j=0;j<G1M;j++)
				{
					GF[i][j]=G1[i-dev][j];
				}
			}
			for(int x=G1M;x<G1M+G2M;x++)
			{
				GF[i][x]=G2[i][x-G1M];
			}
			
		}
	}
	return;
}

void AlgICSMIMOGPC::myc2d2(double* a,double* b1,double* b2,int na,int nb1,int nb2,double Ts,double* d_a,double* d_b1,double* d_b2)
{
	int i,j;
	double p1[2]={1,-1};
	double p2[2]={1,1};
	for (i=1;i<na;i++)
	{
		a[i-1]=pow(2/Ts,na-i)*a[i-1];
	}
	if (nb1>1)
	{
		for (i=1;i<nb1;i++)
		{
			b1[i-1]=pow(2/Ts,nb1-i)*b1[i-1];
		}
	}
	if (nb2>1)
	{
		for (i=1;i<nb2;i++)
		{
			b2[i-1]=pow(2/Ts,nb2-i)*b2[i-1];
		}
	}

	double** A,**B1,**B2;
	malloc_nm(na,na,sizeof(double),sizeof(double*),(double***)&A);
	malloc_nm(na,na,sizeof(double),sizeof(double*),(double***)&B1);
	malloc_nm(na,na,sizeof(double),sizeof(double*),(double***)&B2);
	for(i=1;i<=na;i++)
	{
		A[i-1][0]=a[i-1];
	}
	for(i=1;i<=nb1;i++)
	{
		B1[na-nb1+i-1][0]=b1[i-1];
	}
	for(i=1;i<=nb2;i++)
	{
		B2[na-nb2+i-1][0]=b2[i-1];
	}

	for(i=1;i<=na;i++)
	{
		for (j=1;j<=na-i;j++)
		{
			double* tmp1;
			tmp1=(double*)malloc(j*sizeof(double));
			for(int j1=0;j1<j;j1++)
			{
				tmp1[j1]=A[i-1][j1];
			}
			conv(tmp1,p1,j,2,A[i-1]);
			for(int j1=0;j1<j;j1++)
			{
				tmp1[j1]=B1[i-1][j1];
			}
			conv(tmp1,p1,j,2,B1[i-1]);
			for(int j1=0;j1<j;j1++)
			{
				tmp1[j1]=B2[i-1][j1];
			}
			conv(tmp1,p1,j,2,B2[i-1]);
			free(tmp1);
		}
		if(i>1)
		{
			for(j=1;j<=i-1;j++)
			{
				double *tmp1;
				tmp1=(double*)malloc((na-i+j)*sizeof(double));
				for(int j1=0;j1<na-i+j;j1++)
				{
					tmp1[j1]=A[i-1][j1];
				}
				conv(tmp1,p2,na-i+j,2,A[i-1]);
				for(int j1=0;j1<na-i+j;j1++)
				{
					tmp1[j1]=B1[i-1][j1];
				}
				conv(tmp1,p2,na-i+j,2,B1[i-1]);
				for(int j1=0;j1<na-i+j;j1++)
				{
					tmp1[j1]=B2[i-1][j1];
				}
				conv(tmp1,p2,na-i+j,2,B2[i-1]);
				free(tmp1);
			
			
			}
		}
	}
	
	for (i=2;i<=na;i++)
	{
		for(j=1;j<=na;j++)
		{
			A[0][j-1]=A[0][j-1]+A[i-1][j-1];
			B1[0][j-1]=B1[0][j-1]+B1[i-1][j-1];
			B2[0][j-1]=B2[0][j-1]+B2[i-1][j-1];
		}
	}
	double k=A[0][0];
	if (k!=0)
	{
		for (i=1;i<=na;i++)
		{
			A[0][i-1]=A[0][i-1]/k;
			d_a[i-1]=A[0][i-1];
			B1[0][i-1]=B1[0][i-1]/k;
			d_b1[i-1]=B1[0][i-1];
			B2[0][i-1]=B2[0][i-1]/k;
			d_b2[i-1]=B2[0][i-1];
		}
	}
	
	free_nm(na,(double**)A);
	free_nm(na,(double**)B1);
	free_nm(na,(double**)B2);
}


