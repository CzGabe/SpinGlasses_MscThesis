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

//[12] keplet: MMC eseten energia-különbseget szamol az elozo lepeshez kepest es csinal egy spin-flop-ot, ha szukseges
void spinFlop_MMC_inHomField(vector<vector<bool> >  J_pos, vector<vector<bool> >  J_val, vector<bool>& s1, vector<bool>& s2, double betha, vector<double> h, MTRand& MTrnd)
{
	int N = s1.size();
	double dE1 = 0; double dE2 = 0;
	double expBdE1; double expBdE2;
	
	int i = int(floor(MTrnd() * (N)));	//rnd kivalasztjuk az i-edik spint

	//Energiaváltozás			Az MMC reszig ekvivalens ezzel: for (int j=0; j<N; j++) {dE += 1/N(J[i][j] * s[i] * s[j]);}  dE += h * s[i];
	for (int j=0; j<N; j++)
	{
		if(J_pos[i][j] == true)
		{
			if((J_val[i][j] ^ (s1[i] ^ s1[j])) == false)	//'==false', mert igy vesszuk figyelembe, hogy az i-edig spint megforditottuk
				dE1-=2;
			else
				dE1+=2;

			if((J_val[i][j] ^ (s2[i] ^ s2[j])) == false)	//'==false', mert igy vesszuk figyelembe, hogy az i-edig spint megforditottuk
				dE2-=2;
			else
				dE2+=2;
		}
	}

	dE1 /= sqrt(double(N));
									
	//Metropolis-Monte-carlo
	expBdE1 = exp((-1.0)*betha*dE1);
	
	if((dE1 < 0) || (MTrnd() < expBdE1)) 
		s1[i] = s1[i] ^ true;		//A spin atforditasa (bool negalasa ekvivalens ezzel: s_tmp[i] = (s_tmp[i] == true) ? false : true;)

	dE2 /= sqrt(double(N));

	if (s2[i] == false)
		dE2 -= 2*h[i];
	else
		dE2 += 2*h[i];
									
	//Metropolis-Monte-carlo
	expBdE2 = exp((-1.0)*betha*dE2);
	
	if((dE2 < 0) || (MTrnd() < expBdE2)) 
		s2[i] = s2[i] ^ true;		//A spin atforditasa (bool negalasa ekvivalens ezzel: s_tmp[i] = (s_tmp[i] == true) ? false : true;)

	
}

void giveSpinDiff(vector<bool> spin1, vector<bool> spin2, int& spinDiffTime, vector<bool>& spinDiff, int actualTime)
{
	int N=spin1.size();
	bool cont;
	for (int i=0; i<N; i++)
		if ((spin1[i] ^ spin2[i]) == true) spinDiff[i] = false;

	if (spinDiffTime == 0)
	{
		cont = true;
		for (int i=0; i<N; i++)
			if (spinDiff[i] == true) cont = false;

		if (cont == true) spinDiffTime = actualTime;
	}
}

vector<vector<int> > getCompareMatrix_inhomField(vector<bool> spin1, vector<bool> spin2, vector<bool> spinDiff, int sideX, int sideY)
{
	vector<int> vecGray(2*sideX+2,127);
	vector<vector<int> > Comp(2*sideY+2,vecGray);
	int N = spin1.size();
	int sizeX_inc = sideX+2;
	int sizeY_inc = sideY+2;
	
	//bal felul: Alap rendszer
	int count = 0;
	int i = 0;
	int j = 0;
	for (int s=0; s<N; s++)
	{
		count++;

		Comp[i][j] = (spin1[s] == true) ? 255 : 0;

		j++;

		if(count==sideX)
		{
			i++;
			count = 0;
			j = 0;
		}
	}

	//jobb felul: 1 spinPinned rendszer
	count = 0;
	i = 0;
	j = 0;
	for (int s=0; s<N; s++)
	{
		count++;

		Comp[i][j+sizeX_inc] = (spin2[s] == true) ? 255 : 0;
		//if (s == chosenSpin) Comp[i][j+sizeX_inc] = chSpinVal;

		j++;

		if(count==sideX)
		{
			i++;
			count = 0;
			j = 0;
		}
	}



	//bal alul: A 2 rendszer (kulonbsege) elemeinek szorzata
	count = 0;
	i = 0;
	j = 0;

	for (int s=0; s<N; s++)
	{
		count++;

		Comp[i+sizeY_inc][j] = ((spin1[s] ^ spin2[s]) == false) ? 255 : 0;		//Comp[i][j] = spin1[s] * spin2[s] --> feher, ha +1

		j++;

		if(count==sideX)
		{
			i++;
			count = 0;
			j = 0;
		}
	}

	//jobb alul: Teljes atfedes
	count = 0;
	i = 0;
	j = 0;

	for (int s=0; s<N; s++)
	{
		count++;

		Comp[i+sizeY_inc][j+sizeX_inc] = (spinDiff[s] == true) ? 255 : 0;

		j++;

		if(count==sideX)
		{
			i++;
			count = 0;
			j = 0;
		}
	}


	return Comp;
}

int main()
{
	///--------------Parameterek beolvasasa-----------------///
	cout << "Parameterek beolvasasa" << endl << endl;
	int N;
	string fileP = "05_SG_Param_MMC_vsInhomField.dat";
	string fileJ = loadParameter_string(fileP,3);
	string fileSpin = loadParameter_string(fileP,6);
	string fileH = "";
	
	if (fileJ == "nem") 
		N = loadParameter_int(fileP,9);
	else
		N = load_N_fromJfile(fileJ);

	double Tmin = loadParameter_double(fileP,12);
	double Tmax = loadParameter_double(fileP,15);
	double dT = loadParameter_double(fileP,18);
	fileH = loadParameter_string(fileP,21);
	double kB = loadParameter_double(fileP,24);
	int timeMax = loadParameter_int(fileP,27);
	double bin_C = loadParameter_double(fileP,30);
	int eSampling = loadParameter_int(fileP,33);
	int mSampling = loadParameter_int(fileP,36);
	int qEASampling = loadParameter_int(fileP,39); 
	int qABSampling = loadParameter_int(fileP,42); 
	double bin_qAB = loadParameter_double(fileP,45);
	int spinX = loadParameter_int(fileP,48);
	int spinY = loadParameter_int(fileP,51);
	int spinPixelSize = loadParameter_int(fileP,54);
	
	int eSize =  int(floor((double(timeMax)/double(eSampling))))+1;		//energia vektor merete; +1 a 0. idopillanat miatt
	int mSize = int(floor((double(timeMax)/double(mSampling))))+1;		//magnesezettseg vektor merete; +1 a 0. idopillanat miatt
	int qEAsize = int(floor((double(timeMax)/double(qEASampling))))+1;	//qEA vektor merete; +1 a 0. idopillanat miatt
	int qABsize = int(floor((double(timeMax)/double(qABSampling))))+1;	//qEA vektor merete; +1 a 0. idopillanat miatt

	Tmax += dT/2.0;
	 
	list<string> infos;

	//*** Vektorok deklaralasa ***//
	vector<double> vecZeros_D(N,0.0);
	vector<int> vecZeros_I(N,0);
	vector<bool> vecFalses(N,false);

	vecZeros_D.resize(eSize,0.0); vector<vector<double> > E1(2,vecZeros_D); vector<vector<double> > E2(2,vecZeros_D);
	vecZeros_I.resize(mSize,0); vector<vector<int> > m1(2,vecZeros_I); vector<vector<int> > m2(2,vecZeros_I);
	vecZeros_D.resize(mSize,0.0); vector<vector<double> > m_norm1(2,vecZeros_D); vector<vector<double> > m_norm2(2,vecZeros_D);
	vecZeros_D.resize(qEAsize,0.0); vector<vector<double> > q_EA1(2,vecZeros_D); vector<vector<double> > q_EA2(2,vecZeros_D);
	vecZeros_D.resize(qABsize,0.0); vector<vector<double> > q_AB(2,vecZeros_D);
	
	vecZeros_D.resize(N,0.0); vecZeros_I.resize(N,0);
	vector<double> hField(N,0.0);
	vector<vector<double> >  C(N,vecZeros_D);
	vector<double> C_CDF(int((N*(N-1))/2),0.0);
	vector<bool> spin1(N,false);
	vector<bool> spin2(N,false);
	vector<bool> spin_1st(N,false);
	vector<bool> spinDiff(N,true);
	vector<vector<int> >  J(N,vecZeros_I);
	vector<vector<bool> >  J_pos(N,vecFalses);
	vector<vector<bool> >  J_val(N,vecFalses);
	vector<int> s_sum1(N,0);
	vector<int> s_sum2(N,0);
	vector<vector<int> >  s2_sum1(N,vecZeros_I);
	vector<vector<int> >  s2_sum2(N,vecZeros_I);
	vector<double> s_avg1(N,0.0);
	vector<double> s_avg2(N,0.0);
	vector<vector<double> >  s2_avg1(N,vecZeros_D);
	vector<vector<double> >  s2_avg2(N,vecZeros_D);	
	vector<double> effField1(N,0.0);
	vector<double> effField2(N,0.0);

	vecZeros_I.clear(); vecZeros_I.resize(2*spinX+2,0);
	vector<vector<int> > spinCompareM(2*spinY+2,vecZeros_I);

	int lineMin; int lineMax;
	double corrMin; double corrMax;
	
	int mBinNum = N+1;
	vecZeros_I.resize(mBinNum,0); vector<vector<int> >  m_PDDF1(2,vecZeros_I); vector<vector<int> >  m_PDDF2(2,vecZeros_I);
	vecZeros_D.resize(mBinNum,0.0); vector<vector<double> >  m_PDDF_norm1(2,vecZeros_D); vector<vector<double> >  m_PDDF_norm2(2,vecZeros_D);
	
	vecZeros_D.resize(N,0.0); vecZeros_I.resize(N,0); vecFalses.resize(N,false);

	setXaxis_forVecInMat_D(E1,0.0,double(eSampling)); setXaxis_forVecInMat_D(E2,0.0,double(eSampling));
	setXaxis_forVecInMat_D(m_norm1,0.0,double(mSampling)); setXaxis_forVecInMat_D(m_norm2,0.0,double(mSampling));
	setXaxis_forVecInMat_I(m1,0,mSampling); setXaxis_forVecInMat_I(m2,0,mSampling);
	setXaxis_forVecInMat_D(q_EA1,0.0,double(qEASampling)); setXaxis_forVecInMat_D(q_EA2,0.0,double(qEASampling));
	setXaxis_forVecInMat_I(m_PDDF1,(-1)*N,2); setXaxis_forVecInMat_I(m_PDDF2,(-1)*N,2);
	setXaxis_forVecInMat_D(m_PDDF_norm1,-1.0,(2.0/double(N))); setXaxis_forVecInMat_D(m_PDDF_norm2,-1.0,(2.0/double(N)));
	setXaxis_forVecInMat_D(q_AB,0.0,double(qABSampling));
	
	double betha;
	int spinDiffTime = 0;
	bool snapshotCheck;
	
	string strT;
	string strN = "_N"+convert_ItoS(N);
	string frTriang;
	MTRand MTrnd;		//MT rnd szam generator
	
	//J beolvasasa / generalasa
	if (fileJ != "nem")
		J = load_Jint(N,fileJ);
	else
		J = create_rndJint_CG(N,MTrnd);

	//J felbontasa poziciora (J_pos[i][j]: 1, ha van kapcsolat; 0, ha nincs) es ertekre (J_val[i][j]: 1, ha +1-es a kapcsolat; 0, ha -1-es)
	giveVal_and_Pos_matrixes(J,J_pos,J_val);

	//kiindulasi Spin-vektor beolvasasa / generalasa
	if (fileSpin != "nem")
		spin_1st = load_SpinBool_CG(N,fileSpin);
	else
		spin_1st = create_rndSpinBool_CG(N,MTrnd);

	if(fileSpin == "nem") saveVecBool(spin_1st,"spinIn_N" + strN,"");		//kiindulasi vektor kimentese

	hField = load_inhomogenH(N, fileH);
	
	cout << "MMC modszer inditasa N=" + convert_ItoS(N) + " eseten" << endl;
	cout << endl;
		
	for (double t=Tmin; t<=Tmax; t+=dT)
	{
		strT = "_T"+convert_DtoFix2S(t);		//T homerseklet atirasa a file nevekhez
		betha = 1/(kB * t);

		///----------inicializalas-----------///
		cout << "(Inicializalas, elokeszites...)" << endl;
		infos.clear(); snapshotCheck = true;
		spinDiff.clear(); spinDiff.resize(N,true); spinDiffTime = 0;
		s_sum1.clear(); s_sum1.resize(N,0); s_sum2.clear(); s_sum2.resize(N,0);
		s2_sum1.clear(); s2_sum1.resize(N,vecZeros_I); s2_sum2.clear(); s2_sum2.resize(N,vecZeros_I);
		effField1.clear(); effField1.resize(N,0.0); effField2.clear(); effField2.resize(N,0.0);
		for (int i=0; i<mBinNum; i++)
				m_PDDF1[1][i] = m_PDDF2[1][i] = 0;

		spin1 = spin_1st; spin2 = spin_1st;				//Minten adott N-re es T-re ugyanabbol a spinallapotbol indulunk
		spinCompareM = getCompareMatrix_inhomField(spin1,spin2,spinDiff,spinX,spinY);
		saveMatInt_meshPGM(spinCompareM,"spinCompare" + strN + strT + "_" + convert_Ito5str(0),spinPixelSize);

	///-------------ido fejlesztes----------------///
		cout << "Rendszer idofejlesztese T = " << t  << " eseten" << endl << "\"|\" = 100 000 SpinFlop ~ " + convert_DtoS(100000.0/double(N)) + " sweep: ";
		for (int time=0; time<timeMax; time++)	
		{
			if (time % 100000 == 0) cout << "|";
			spinFlop_MMC_inHomField(J_pos,J_val,spin1,spin2,betha,hField,MTrnd);
			giveSpinSum_MMC(spin1,s_sum1,s2_sum1); giveSpinSum_MMC(spin2,s_sum2,s2_sum2);	//[2][3] atlagok szamitasahoz osszegzes
			giveEffFieldSum_MMC(effField1,J_pos,J_val,spin1); giveEffFieldSum_MMC(effField2,J_pos,J_val,spin2);	//[13] belso eff.terhez szum

			if (time % eSampling == 0)
			{
				E1[1][int(time/eSampling)] = getEnergy_justPairing(J_pos,J_val,spin1);	//[1] energia szamitasa
				E2[1][int(time/eSampling)] = getEnergy(J_pos,J_val,spin2,hField);
			}
			if (time % mSampling == 0)
			{
				m1[1][int(time/mSampling)] = getMagnetism(spin1);					//[4] magnesseg szamitasa 
				m2[1][int(time/mSampling)] = getMagnetism(spin2);
			}
			if (time % qEASampling == 0)
			{
				q_EA1[1][int(time/qEASampling)] = getQEA_sSUM(time,s_sum1);		//[16]
				q_EA2[1][int(time/qEASampling)] = getQEA_sSUM(time,s_sum2);
			}
			if (time % qABSampling == 0)
			{
				q_AB[1][int(time/qABSampling)] = getQAB(spin1,spin2);
			}
			
			if (snapshotCheck == true)
			{
				giveSpinDiff(spin1,spin2,spinDiffTime,spinDiff,time);
				spinCompareM = getCompareMatrix_inhomField(spin1,spin2,spinDiff,spinX,spinY);
				//saveVecBool(spin1,"spinSnapshot" + strN + strT + "_" + convert_Ito5str(time),"");
				//saveVecBool(spin2,"spinSnapshot" + strN + strT + "_" + convert_Ito5str(time) + chosenSpinIndex,"");
				//if (time <= 1000)
				saveMatInt_meshPGM(spinCompareM,"spinCompare" + strN + strT + "_" + convert_Ito5str(time+1),spinPixelSize);
				if (spinDiffTime != 0) snapshotCheck = false;
			}
		}
			
		cout << endl;
			
		///-----Atlagok, belso ter es megnesseg-----///
		giveSpinAVG_MMC(timeMax,s_sum1,s2_sum1,s_avg1,s2_avg1); giveSpinAVG_MMC(timeMax,s_sum2,s2_sum2,s_avg2,s2_avg2);
		giveEffFieldAvg_MMC(effField1,timeMax); giveEffFieldAvg_MMC(effField2,timeMax);
		
		for (int i=0; i<mSize; i++)
		{
			m_norm1[1][i] = double(m1[1][i])/double(N);
			m_norm2[1][i] = double(m2[1][i])/double(N);
		}

		getPDDF_forM(m1,m_PDDF1);
		getPDDF_forM(m2,m_PDDF2);

		for (int i=0; i<mBinNum; i++)
		{
			m_PDDF_norm1[1][i] = double(m_PDDF1[1][i])/double(mSize);
			m_PDDF_norm2[1][i] = double(m_PDDF2[1][i])/double(mSize);
		}

		///-----Adatok kimentese------///
		cout << "Alap rendszer adatainak kimentese" << endl;

		saveMatDoub(E1,"E" + strN + strT,"Energy values, sampling freq="+convert_ItoS(eSampling));
		saveMatDoub(m_norm1,"m" + strN + strT, "Magnetism values, sampling freq="+convert_ItoS(mSampling));
		saveMatDoub(m_PDDF_norm1,"m" + strN + strT + "_PDDF","Magnesezettseg Vsz.-suruseg eo.-fv.");
		saveMatDoub(q_EA1,"q_EA" + strN + strT,"q_EA, Edward-Anderson order-parameter, sampling freq="+convert_ItoS(qEASampling));
		saveVecDoub(s_avg1,"s_avg" + strN + strT,"<si> vector");
		saveVecDoub(effField1,"effField" + strN + strT,"");
		
		infos.push_back("Atfedes teljes elvesztesenek ideje:"); infos.push_back(convert_ItoS(spinDiffTime));
		infos.push_back("q_[alfa-beta] [futas vegen]:"); infos.push_back(convert_DtoS(getQAB(spin1,spin2)));
		infos.push_back(""); infos.push_back("");
		infos.push_back("Alap rendszer:"); infos.push_back("");
		infos.push_back("<m>:"); infos.push_back(convert_DtoS(getVectorsAVG(s_avg1)));
		infos.push_back("q_EA [futas vegen]:"); infos.push_back(convert_DtoS(getQEA_sAVG(s_avg1)));

		//Korrelacio szamitasa es ezzel kapcsolatos eredmenyek kimentese
		cout << "	- Kapcsolt-Korrelacio szamitasa & adatmentes" << endl;
		C = getConCorrMatrix(s_avg1,s2_avg1);	//[9]
		C_CDF = getUpperTri_fromSqMatrix_doub(C);
		giveExtrLines_of_matrix(C,lineMin,lineMax,corrMin,corrMax);
		infos.push_back("Kapcsolt-Korrelacioban leggyengebb spin [i]:"); infos.push_back(convert_ItoS(lineMin+1));
		infos.push_back("sum(|c[i]|) ="); infos.push_back(convert_DtoS(corrMin));
		infos.push_back("<s_i> ="); infos.push_back(convert_DtoS(s_avg1[lineMin]));
		infos.push_back("effField(i) ="); infos.push_back(convert_DtoS(effField1[lineMin]));
		infos.push_back("Kapcsolt-Korrelacioban legerosebb spin [j]:"); infos.push_back(convert_ItoS(lineMax+1));
		infos.push_back("sum(|c[j]|) ="); infos.push_back(convert_DtoS(corrMax));
		infos.push_back("<s_j> ="); infos.push_back(convert_DtoS(s_avg1[lineMax]));
		infos.push_back("effField(j) ="); infos.push_back(convert_DtoS(effField1[lineMax]));
	
		saveMatDoub(C,"C" + strN + strT,"Connected-Correlation C Matrix");
		saveMatDoub3d(C,"C" + strN + strT + "_3D","Connected-Correlation Matrix for 3D plot");
		saveVecDoub(C_CDF,"C" + strN + strT + "_CDF","C Matrix's UpperTriangular for CDF");
		savePDDF_vec(C_CDF,bin_C,-1,1,"C" + strN + strT + "_PDDF","C matrix Vsz.-suruseg eo.-fv.");
		savePDDF_vec_alter(C_CDF,bin_C,-1,1,"C" + strN + strT + "_PDDF_alter","C matrix Vsz.-suruseg eo.-fv.");
	
		//Kovariancia szamitasa es ezzel kapcsolatos eredmenyek kimentese
		cout << "	- Kovariancia szamitasa & adatmentes" << endl;
		C = getCorrMatrix(s2_avg1);		//[10]
		C_CDF = getUpperTri_fromSqMatrix_doub(C);
		giveExtrLines_of_matrix(C,lineMin,lineMax,corrMin,corrMax);
		infos.push_back("Korrelacioban leggyengebb spin [k]:"); infos.push_back(convert_ItoS(lineMin+1));
		infos.push_back("sum(|c[k]|) ="); infos.push_back(convert_DtoS(corrMin));
		infos.push_back("<s_k> ="); infos.push_back(convert_DtoS(s_avg1[lineMin]));
		infos.push_back("effField(k) ="); infos.push_back(convert_DtoS(effField1[lineMin]));
		infos.push_back("Korrelacioban legerosebb spin [l]:"); infos.push_back(convert_ItoS(lineMax+1));
		infos.push_back("sum(|c[l]|) ="); infos.push_back(convert_DtoS(corrMax));
		infos.push_back("<s_l> ="); infos.push_back(convert_DtoS(s_avg1[lineMax]));
		infos.push_back("effField(l) ="); infos.push_back(convert_DtoS(effField1[lineMax]));

		saveMatDoub(C,"K" + strN + strT,"Correlation [K] Matrix");
		saveMatDoub3d(C,"K" + strN + strT + "_3D","Correlation Matrix for 3D plot");
		saveVecDoub(C_CDF,"K" + strN + strT + "_CDF","Correlation Matrix's UpperTriangular for CDF");
		savePDDF_vec(C_CDF,bin_C,-1,1,"K" + strN + strT + "_PDDF","Correlation matrix Vsz.-suruseg eo.-fv.");
		savePDDF_vec_alter(C_CDF,bin_C,-1,1,"K" + strN + strT + "_PDDF_alter","Correlation matrix Vsz.-suruseg eo.-fv.");

		cout << "Alternativ rendszer adatainak kimentese" << endl;

		saveMatDoub(E2,"E" + strN + strT + "_inhomField","Energy values, sampling freq="+convert_ItoS(eSampling));
		saveMatDoub(m_norm2,"m" + strN + strT + "_inhomField", "Magnetism values, sampling freq="+convert_ItoS(mSampling));
		saveMatDoub(m_PDDF_norm2,"m" + strN + strT + "_inhomField" + "_PDDF","Magnesezettseg Vsz.-suruseg eo.-fv.");
		saveMatDoub(q_EA2,"q_EA" + strN + strT + "_inhomField","q_EA, Edward-Anderson order-parameter, sampling freq="+convert_ItoS(qEASampling));
		saveVecDoub(s_avg2,"s_avg" + strN + strT + "_inhomField","<si> vector");
		saveVecDoub(effField2,"effField" + strN + strT + "_inhomField","");
		
		infos.push_back(""); infos.push_back("");
		infos.push_back("Alternativ rendszer:"); infos.push_back("inhomogen ter");
		infos.push_back("<m>:"); infos.push_back(convert_DtoS(getVectorsAVG(s_avg2)));
		infos.push_back("q_EA [futas vegen]:"); infos.push_back(convert_DtoS(getQEA_sAVG(s_avg2)));

		cout << "	- Kapcsolt-Korrelacio szamitasa & adatmentes" << endl;
		C = getConCorrMatrix(s_avg2,s2_avg2);	//[9]
		C_CDF = getUpperTri_fromSqMatrix_doub(C);
		giveExtrLines_of_matrix(C,lineMin,lineMax,corrMin,corrMax);
		infos.push_back("Kapcsolt-Korrelacioban leggyengebb spin [i]:"); infos.push_back(convert_ItoS(lineMin+1));
		infos.push_back("sum(|c[i]|) ="); infos.push_back(convert_DtoS(corrMin));
		infos.push_back("<s_i> ="); infos.push_back(convert_DtoS(s_avg2[lineMin]));
		infos.push_back("effField(i) ="); infos.push_back(convert_DtoS(effField2[lineMin]));
		infos.push_back("Kapcsolt-Korrelacioban legerosebb spin [j]:"); infos.push_back(convert_ItoS(lineMax+1));
		infos.push_back("sum(|c[j]|) ="); infos.push_back(convert_DtoS(corrMax));
		infos.push_back("<s_j> ="); infos.push_back(convert_DtoS(s_avg2[lineMax]));
		infos.push_back("effField(j) ="); infos.push_back(convert_DtoS(effField2[lineMax]));

		saveMatDoub(C,"C" + strN + strT + "_inhomField","Connected-Correlation C Matrix");
		saveMatDoub3d(C,"C" + strN + strT + "_3D" + "_inhomField","Connected-Correlation Matrix for 3D plot");
		saveVecDoub(C_CDF,"C" + strN + strT + "_CDF" + "_inhomField","C Matrix's UpperTriangular for CDF");
		savePDDF_vec(C_CDF,bin_C,-1,1,"C" + strN + strT + "_PDDF" + "_inhomField","C matrix Vsz.-suruseg eo.-fv.");
		savePDDF_vec_alter(C_CDF,bin_C,-1,1,"C" + strN + strT + "_PDDF_alter" + "_inhomField","C matrix Vsz.-suruseg eo.-fv.");

		//Kovariancia szamitasa es ezzel kapcsolatos eredmenyek kimentese
		cout << "	- Korrelacio szamitasa & adatmentes" << endl;
		C = getCorrMatrix(s2_avg2);		//[10]
		C_CDF = getUpperTri_fromSqMatrix_doub(C);
		giveExtrLines_of_matrix(C,lineMin,lineMax,corrMin,corrMax);
		infos.push_back("Korrelacioban leggyengebb spin [k]:"); infos.push_back(convert_ItoS(lineMin+1));
		infos.push_back("sum(|c[k]|) ="); infos.push_back(convert_DtoS(corrMin));
		infos.push_back("<s_k> ="); infos.push_back(convert_DtoS(s_avg2[lineMin]));
		infos.push_back("effField(k) ="); infos.push_back(convert_DtoS(effField2[lineMin]));
		infos.push_back("Korrelacioban legerosebb spin [l]:"); infos.push_back(convert_ItoS(lineMax+1));
		infos.push_back("sum(|c[l]|) ="); infos.push_back(convert_DtoS(corrMax));
		infos.push_back("<s_l> ="); infos.push_back(convert_DtoS(s_avg2[lineMax]));
		infos.push_back("effField(l) ="); infos.push_back(convert_DtoS(effField2[lineMax]));

		saveMatDoub(C,"K" + strN + strT + "_inhomField","Correlation K Matrix");
		saveMatDoub3d(C,"K" + strN + strT + "_3D" + "_inhomField","Correlation Matrix for 3D plot");
		saveVecDoub(C_CDF,"K" + strN + strT + "_CDF" + "_inhomField","Correlation Matrix's UpperTriangular for CDF");
		savePDDF_vec(C_CDF,bin_C,-1,1,"K" + strN + strT + "_PDDF" + "_inhomField","Correlation matrix Vsz.-suruseg eo.-fv.");
		savePDDF_vec_alter(C_CDF,bin_C,-1,1,"K" + strN + strT + "_PDDF_alter" + "_inhomField","Correlation matrix Vsz.-suruseg eo.-fv.");

		cout << "qAB es info-file kimentese" << endl;
		saveMatDoub(q_AB,"q_AB" + strN + strT,"q_(Alpha-Betha) parameter" + convert_ItoS(qABSampling));
		savePDDF_mat(q_AB,1,bin_qAB,-1.0,+1.0,"q_AB" + strN + strT + "_PDDF","q_(Alpha-Betha) parameter" + convert_ItoS(qABSampling));
		saveLS_infos(infos,"infos" + strN + strT);

		cout << endl;
	}
			
	cout << endl;
	cout << "Homerseklet fuggetlen adatok kimentese" << endl << endl;		

	double fi_J = getJfrustrate_int(J);		//[11] J frusztraltsaganak kiszamitasa
	saveMatInt(J,"J" + strN + "_inhomField","J_Matrix, Frustration Rate: " + convert_DtoS(fi_J));
	saveMatInt3d(J,"J" + strN + "_3D","J Matrix 3D plothoz");
	saveMatInt_meshPGM_forJ(J,"J" + strN,20);

	system("pause");
}
