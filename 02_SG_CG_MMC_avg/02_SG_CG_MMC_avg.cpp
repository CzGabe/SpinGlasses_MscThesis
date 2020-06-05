//#include "stdafx.h"
#include <fstream>
#include <list>
#include <vector>
#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>
//***Mersenne Twister Random Number Generator
#include "MersenneTwister.h"
//***Sajat header-ek
#include "saveFiles.cpp"
#include "readFromFiles.cpp"
#include "variousFunctions.cpp"
#include "equationFunctions.cpp"
#include "spec_MMC_functions.cpp"

using namespace std;

//PDDF fuggvenyek kellenek az atlagolashoz
vector<double> getPDDF_forC(vector<vector<double> > C, int distSize, double bin_C)
{
	vector<double> cUpperTri = getUpperTri_fromSqMatrix_doub(C);
	int dataNum = cUpperTri.size();
	double tmp;
	vector<double> cPDDF(distSize,0.0);

	for (int i=0; i<dataNum; i++)
	{
		tmp = -1-(bin_C/2);
		for (int j=0; j<distSize; j++)
		{
			if (cUpperTri[i] >= tmp && cUpperTri[i] <= tmp+bin_C)
			{
				cPDDF[j]++;
				break;
			}
			else
			{
				tmp += bin_C;
			}
		}
	}

	//Normalas
	for (int i=0; i<distSize; i++)
		cPDDF[i] /= dataNum;			//igy az integral 1 lesz...

	//A szeleken fel-binek vannak
	cPDDF[0] *= 2;
	cPDDF[distSize-1] *= 2;

	return cPDDF;
}

vector<double> getPDDF_forC_alter(vector<vector<double> > C, int distSize, double bin_C)
{
	vector<double> cUpperTri = getUpperTri_fromSqMatrix_doub(C);
	int dataNum = cUpperTri.size();
	double tmp;
	vector<double> cPDDF(distSize,0.0);

	for (int i=0; i<dataNum; i++)
	{
		tmp = -1;
		for (int j=0; j<distSize; j++)
		{
			if (cUpperTri[i] >= tmp && cUpperTri[i] <= tmp+bin_C)
			{
				cPDDF[j]++;
				break;
			}
			else
			{
				tmp += bin_C;
			}
		}
	}

	//Normalas
	for (int i=0; i<distSize; i++)
		cPDDF[i] /= dataNum;			//igy az integral 1 lesz...

	return cPDDF;
}

int main()
{
	///--------------Parameterek beolvasasa-----------------///
	cout << "Parameterek beolvasasa" << endl << endl;
	string fileP = "02_SG_Param_MMC_avg.dat";
	int N = loadParameter_int(fileP,3);	
	string fileSpin = loadParameter_string(fileP,6);
	double Tmin = loadParameter_double(fileP,9);
	double Tmax = loadParameter_double(fileP,12);
	double dT = loadParameter_double(fileP,15);
	double Hmin = loadParameter_double(fileP,18);
	double Hmax = loadParameter_double(fileP,21);
	double dH = loadParameter_double(fileP,24);
	double kB = loadParameter_double(fileP,27);
	int timeMax = loadParameter_int(fileP,30);
	double bin_C = loadParameter_double(fileP,33);
	int mSampling = loadParameter_int(fileP,36);
	int eSampling = loadParameter_int(fileP,39);
	int qEASampling = loadParameter_int(fileP,42); 
	int avgNum_overJ = loadParameter_int(fileP,45);

	int pairingNum = int((N*(N-1))/2);
	int eSize =  int(floor((double(timeMax)/double(eSampling))))+1;		//energia vektor merete; +1 a 0. idopillanat miatt
	int mSize = int(floor((double(timeMax)/double(mSampling))))+1;		//magnesezettseg vektor merete; +1 a 0. idopillanat miatt
	int qEAsize = int(floor((double(timeMax)/double(qEASampling))))+1;	//qEA vektor merete; +1 a 0. idopillanat miatt

	Tmax += dT/2.0; Hmax += dH/2.0;

	//*** Vektorok deklaralasa ***//
	vector<double> vecZeros_D(N,0.0);
	vector<int> vecZeros_I(N,0);
	vector<bool> vecFalses(N,false);

	vecZeros_D.resize(eSize,0.0); vector<vector<double> > E(2,vecZeros_D);
	vecZeros_I.resize(mSize,0); vector<vector<int> > m(2,vecZeros_I);
	vecZeros_D.resize(mSize,0.0); vector<vector<double> > m_avg(2,vecZeros_D);		//atlagolva es a vegen normalva is
	vecZeros_D.resize(qEAsize,0.0); vector<vector<double> > q_EA(2,vecZeros_D);

	vecZeros_D.resize(N,0.0); vecZeros_I.resize(N,0);
	vector<double> hField(N,0.0);
	vector<vector<double> >  K(N,vecZeros_D);		//K=<sisj>
	vector<vector<double> >  C(N,vecZeros_D);		//C=<sisj>-<si><sj>
	vector<double> C_CDF(pairingNum,0.0);

	vecZeros_D.resize(pairingNum,0.0);
	vector<vector<double> > K_ValuesAll(avgNum_overJ, vecZeros_D);
	vector<vector<double> > C_ValuesAll(avgNum_overJ, vecZeros_D);
	vecZeros_D.resize(N,0.0);

	vector<bool> spin(N,false);
	vector<bool> spin_1st(N,false);
	vector<vector<int> >  J(N,vecZeros_I);
	vector<vector<bool> >  J_pos(N,vecFalses);
	vector<vector<bool> >  J_val(N,vecFalses);
	vector<int> s_sum(N,0);
	vector<vector<int> >  s2_sum(N,vecZeros_I);
	vector<double> s_avg(N,0.0);
	vector<vector<double> >  s2_avg(N,vecZeros_D);
	
	int H_num = int(((Hmax-Hmin)/dH)+1);
	int T_num = int(((Tmax-Tmin)/dT)+1);

	vecZeros_D.resize(H_num,0.0);
	vector<vector<double> >  C_Var(T_num,vecZeros_D);		//szoras
	vector<vector<double> >  K_Var(T_num,vecZeros_D);		//szoras
	vector<vector<double> >  d_CkJ(T_num,vecZeros_D);
	vector<vector<double> >  d_KkJ(T_num,vecZeros_D);
	vecZeros_D.resize(N,0.0);
	int Tindex = -1;
	int Hindex = -1;

	int cBinNum = int(floor((2/bin_C)+1));
	int cBinNum_alter = int(floor(2/bin_C));
	int mBinNum = N+1;

	vecZeros_D.resize(cBinNum,0.0);			//PDDF matrixok alabb
	vector<vector<double> >  C_AVG(2,vecZeros_D);
	vector<vector<double> >  K_AVG(2,vecZeros_D);
	vector<double>  Ctmp_AVG(2,0.0);
	vecZeros_D.resize(cBinNum_alter,0.0);
	vector<vector<double> >  K_AVG_alter(2,vecZeros_D);
	vector<vector<double> >  C_AVG_alter(2,vecZeros_D);
	vecZeros_I.resize(mBinNum,0); vector<vector<int> >  m_PDDF_sum(2,vecZeros_I);
	vecZeros_D.resize(mBinNum,0.0);  vector<vector<double> >  m_PDDF_avg(2,vecZeros_D);	//A vegen ebben kerul m_PDDF_sum normalva
	vecZeros_D.resize(N,0.0); vecZeros_I.resize(N,0);
	
	double betha;
	int sweep10percent = int(floor(double(timeMax/10.0)));

	MTRand MTrnd;		//MT rnd szam generator

	string strN = "_N"+convert_ItoS(N);
	string strT;
	string strAVG = "_avg" + convert_ItoS(avgNum_overJ);
	string strH;

	//kiindulasi Spin-vektor beolvasasa / generalasa
	if (fileSpin != "nem")
		spin_1st = load_SpinBool_CG(N,fileSpin);
	else
		spin_1st = create_rndSpinBool_CG(N,MTrnd);

	setXaxis_forVecInMat_D(E,0.0,double(eSampling));
	setXaxis_forVecInMat_I(m,0,mSampling);
	setXaxis_forVecInMat_D(m_avg,0.0,double(mSampling));
	setXaxis_forVecInMat_D(q_EA,0.0,double(qEASampling));
	setXaxis_forVecInMat_D(C_AVG,-1.0,bin_C);
	setXaxis_forVecInMat_D(K_AVG,-1.0,bin_C);
	setXaxis_forVecInMat_D(C_AVG_alter,-1.0+(bin_C/2.0),bin_C);
	setXaxis_forVecInMat_D(K_AVG_alter,-1.0+(bin_C/2.0),bin_C);
	setXaxis_forVecInMat_I(m_PDDF_sum,(-1)*N,2);
	
	for (double h=Hmin; h<=Hmax; h+=dH)
	{
		cout << "Ter: " << h << endl;
		strH = "_H"+convert_DtoFix2S(h);
		Hindex++; Tindex = -1;
		for (int i=0; i<N; i++)
			hField[i] = h;

		for (double t=Tmin; t<=Tmax; t+=dT)
		{
			cout << "homerseklet: " << t << endl;
			strT = "_T"+convert_DtoFix2S(t);
			Tindex++;
			betha = 1/(kB*t);

			for (int i=0; i<eSize; i++)
				E[1][i] = 0.0;

			for (int i=0; i<mSize; i++)
				m_avg[1][i] = 0.0;
			
			for (int i=0; i<qEAsize; i++)
				q_EA[1][i] = 0.0;

			for (int i=0; i<cBinNum; i++)
				C_AVG[1][i] = K_AVG[1][i] = 0.0;

			for (int i=0; i<cBinNum_alter; i++)
				C_AVG_alter[1][i] = K_AVG_alter[1][i] = 0.0;

			for (int i=0; i<mBinNum; i++)
			{
				m_PDDF_avg[1][i] = 0.0;
				m_PDDF_sum[1][i] = 0;
			}

			for (int avgIndex=0; avgIndex<avgNum_overJ; avgIndex++)
				cout << ".";
			cout << endl;

/*AVG*/		for (int avgIndex=0; avgIndex<avgNum_overJ; avgIndex++)
			{
				cout << "|";
				J = create_rndJint_CG(N,MTrnd);
				giveVal_and_Pos_matrixes(J,J_pos,J_val);
				spin = spin_1st;
				
				//INICIALIZALAS
				s_sum.clear(); s_sum.resize(N,0);	s2_sum.clear(); s2_sum.resize(N,vecZeros_I);

				for (int st=0; st<sweep10percent; st++)			//Elso 10%-nak megfelelo lepes kihagyasa a statisztikabol
					spinFlop_MMC(J_pos,J_val,spin,betha,hField,MTrnd);

				///-------------ido fejlesztes----------------///
				for (int time=0; time<=timeMax; time++)	
				{
					spinFlop_MMC(J_pos,J_val,spin,betha,hField,MTrnd);	//[12] Metropolis Monte Carlo
					giveSpinSum_MMC(spin,s_sum,s2_sum);					//[2][3] atlagok szamitasahoz osszegzes
					
					if (time % eSampling == 0) E[1][int(time/eSampling)] += getEnergy(J_pos,J_val,spin,hField);	//[1] energia szamitasa
					if (time % mSampling == 0) m[1][int(time/mSampling)] = getMagnetism(spin);					//[4] magnesseg szamitasa 
					if (time % qEASampling == 0) q_EA[1][int(time/qEASampling)] += getQEA_sSUM(time,s_sum);		//[16]
				}
			
				giveSpinAVG_MMC(timeMax,s_sum,s2_sum,s_avg,s2_avg);			//[2][3] atlagok szamitasa korrelaciokhoz es magnesseghez
			
				//Connected-Correlation; C=<sisj>-<si><sj>
				C = getConCorrMatrix(s_avg,s2_avg);		//[10]
				C_CDF = getUpperTri_fromSqMatrix_doub(C);
				
				for (int i=0; i<pairingNum; i++)
					C_ValuesAll[avgIndex][i] = C_CDF[i];
				
				//Atlagolas
				Ctmp_AVG = getPDDF_forC(C,cBinNum,bin_C);
				for (int i=0; i<cBinNum; i++)
					C_AVG[1][i] += Ctmp_AVG[i];

				Ctmp_AVG = getPDDF_forC_alter(C,cBinNum_alter,bin_C);
				for (int i=0; i<cBinNum_alter; i++)
					C_AVG_alter[1][i] += Ctmp_AVG[i];

				d_CkJ[Tindex][Hindex] += get_dCkJ(C,J);			//[15]

				//Correlation; K=<sisj>
				K = getCorrMatrix(s2_avg);		//[10]
				C_CDF = getUpperTri_fromSqMatrix_doub(K);
				
				for (int i=0; i<pairingNum; i++)
					K_ValuesAll[avgIndex][i] = C_CDF[i];
				
				//PDDF-ek kerese:
				Ctmp_AVG = getPDDF_forC(K,cBinNum,bin_C);
				for (int i=0; i<cBinNum; i++)
					K_AVG[1][i] += Ctmp_AVG[i];

				Ctmp_AVG = getPDDF_forC_alter(K,cBinNum_alter,bin_C);
				for (int i=0; i<cBinNum_alter; i++)
					K_AVG_alter[1][i] += Ctmp_AVG[i];

				d_KkJ[Tindex][Hindex] += get_dCkJ(K,J);			//[15]

				for (int i=0; i<mSize; i++)
					m_avg[1][i] += double(m[1][i]);

				getPDDF_forM(m,m_PDDF_sum);

			} //End of AVG loop
			cout << endl;

			for (int i=0; i<cBinNum; i++)
			{
				C_AVG[1][i] /= avgNum_overJ;
				K_AVG[1][i] /= avgNum_overJ;
			}

			for (int i=0; i<cBinNum_alter; i++)
			{
				C_AVG_alter[1][i] /= avgNum_overJ;
				K_AVG_alter[1][i] /= avgNum_overJ;
			}

			for (int i=0; i<eSize; i++)
				E[1][i] /= double(avgNum_overJ);

			for (int i=0; i<mSize; i++)
				m_avg[1][i] /= double(avgNum_overJ*N);		//Normalunk N-nel is

			for (int i=0; i<qEAsize; i++)
				q_EA[1][i] /= double(avgNum_overJ);

			for (int i=0; i<mBinNum; i++)
			{
				m_PDDF_avg[0][i] = double(m_PDDF_sum[0][i])/double(N);
				m_PDDF_avg[1][i] = double(m_PDDF_sum[1][i])/double(avgNum_overJ*mSize);
			}

			C_Var[Tindex][Hindex] = get_SqrtVar_ofMat(C_ValuesAll);
			K_Var[Tindex][Hindex] = get_SqrtVar_ofMat(K_ValuesAll);

			cout << "Adatok mentese" << endl;
			saveMatDoub(C_ValuesAll,"C_Values" + strN + strT + strH + strAVG,"Egy oszlopban: egy rnd valasztott J-hez tartozo kapcsolt-korrelcio ertekek (C matrixok felso haromszogei)");
			saveMatDoubToVec(C_ValuesAll,"C_Values" + strN + strT + strH + strAVG + "_CDF","");
			saveMatDoub(C_AVG,"C_PDDF" + strN + strT + strH + strAVG,"Kapcsolt-Korrelacio matrix (C=<sisj>-<si><sj>) PDDF");
			saveMatDoub(C_AVG_alter,"C_PDDF_alter" + strN + strT + strH + strAVG,"Kapcsolt-Korrelacio matrix (C=<sisj>-<si><sj>) PDDF");
			saveMatDoub(K_ValuesAll,"K_Values" + strN + strT + strH + strAVG,"Egy oszlopban: egy rnd valasztott J-hez tartozo korrelacio ertekek (K matrixok felso haromszogei)");
			saveMatDoubToVec(K_ValuesAll,"K_Values" + strN + strT + strH + strAVG + "_CDF","");
			saveMatDoub(K_AVG,"K_PDDF" + strN + strT + strH + strAVG,"Korrelacio matrix (C=<sisj>) PDDF");
			saveMatDoub(K_AVG_alter,"K_PDDF_alter" + strN + strT + strH + strAVG,"Korrelacio matrix (C=<sisj>) PDDF");

			saveMatDoub(E,"E" + strN + strT + strH + strAVG,"");
			saveMatDoub(m_avg,"m" + strN + strT + strH + strAVG,"");
			saveMatDoub(m_PDDF_avg,"m_PDDF" + strN + strT + strH + strAVG,"");
			saveMatDoub(q_EA,"q_EA" + strN + strT + strH + strAVG,"");

			cout << endl;
		} //End of T loop
	} //End of H loop

	for (int i=0; i<T_num; i++)
	{
		for (int j=0; j<H_num; j++)
		{
			C_Var[i][j] /= avgNum_overJ;
			K_Var[i][j] /= avgNum_overJ;
			d_CkJ[i][j] /= avgNum_overJ;
			d_KkJ[i][j] /= avgNum_overJ;
		}
	}

	saveMatDoub(C_Var,"C_Var" + strN + strAVG, "Kapcsolt-Korrelacios eloszlasok szorasa; sor ~ T tengely; oszlop ~ H tengely");
	saveMatDoub3d_wAxes(C_Var,Tmin,Tmax,dT,Hmin,Hmax,dH,"C_Var" + strN + strAVG + "_3D","A szorasok szorasa= "+ convert_DtoS(get_SqrtVar_ofMat(C_Var)) +"; x tengely = T; y tengely = H");

	saveMatDoub(K_Var,"K_Var" + strN + strAVG, "Korrelacios eloszlasok szorasa; sor ~ T tengely; oszlop ~ H tengely");
	saveMatDoub3d_wAxes(K_Var,Tmin,Tmax,dT,Hmin,Hmax,dH,"K_Var" + strN + strAVG + "_3D","A szorasok szorasa= "+ convert_DtoS(get_SqrtVar_ofMat(K_Var)) +"; x tengely = T; y tengely = H");

	saveMatDoub(d_CkJ,"dCkJ" + strN + strAVG, "d_CkJ; sor ~ T tengely; oszlop ~ H tengely");
	saveMatDoub3d_wAxes(d_CkJ,Tmin,Tmax,dT,Hmin,Hmax,dH,"dCkJ" + strN + strAVG + "_3D","x tengely = T; y tengely = H");

	saveMatDoub(d_KkJ,"dKkJ" + strN + strAVG, "d_KkJ; sor ~ T tengely; oszlop ~ H tengely");
	saveMatDoub3d_wAxes(d_KkJ,Tmin,Tmax,dT,Hmin,Hmax,dH,"dKkJ" + strN + strAVG + "_3D","x tengely = T; y tengely = H");

	system("pause");
}

