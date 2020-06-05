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

void spinFlop_MMC_pin(vector<vector<bool> >  J_pos, vector<vector<bool> >  J_val, vector<bool>& s, double betha, int chosenSpin, MTRand& MTrnd)
{
	int N = s.size();
	double dE = 0;
	double expBdE;
	
	//rnd kivalasztunk egy spint, ami nem a megfogott
	int i = int(floor(MTrnd() * N));	//rnd kivalasztjuk az i-edik spint
	if (i == chosenSpin)
		while (i == chosenSpin)
			i = int(floor(MTrnd() * N));

	dE = 0;
		
	for (int j=0; j<N; j++)
	{
		if(J_pos[i][j] == true)
		{
			if((J_val[i][j] ^ (s[i] ^ s[j])) == false)			//'==false', mert igy vesszuk figyelembe, hogy az i-edig spint megforditottuk
				dE-=2;
			else
				dE+=2;
		}
	}

	dE /= sqrt(double(N));
									
		//Metropolis-Monte-carlo
	expBdE = exp((-1.0)*betha*dE);
	
	if((dE < 0) || (MTrnd() < expBdE)) 
		s[i] = s[i] ^ true;				//A spin atforditasa (bool negalasa ekvivalens ezzel: s_tmp[i] = (s_tmp[i] == true) ? false : true;)

}


int main()
{
	///--------------Parameterek beolvasasa-----------------///
	cout << "Parameterek beolvasasa" << endl << endl;
	int N;
	string fileP = "03_SG_Param_MMC_PinSpin.dat";
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
	int chosenSpin = loadParameter_int(fileP,21)-1;
	string chosenSpinIndex = loadParameter_string(fileP,21);
	bool chosenSpinVal = loadParameter_bool(fileP,24,"+1");
	double kB = loadParameter_double(fileP,27);
	int timeMax = loadParameter_int(fileP,30);
	double bin_C = loadParameter_double(fileP,33);
	int eSampling = loadParameter_int(fileP,36);
	int mSampling = loadParameter_int(fileP,39);	
	int qEASampling = loadParameter_int(fileP,42);
	int snapshotFreq = loadParameter_int(fileP,45);
	int snapshotX = loadParameter_int(fileP,48);
	int snapshotY = loadParameter_int(fileP,51);
	int spinPixelSize = loadParameter_int(fileP,54);

	int eSize =  int(floor((double(timeMax)/double(eSampling))))+1;		//energia vektor merete; +1 a 0. idopillanat miatt
	int mSize = int(floor((double(timeMax)/double(mSampling))))+1;		//magnesezettseg vektor merete; +1 a 0. idopillanat miatt
	int qEAsize = int(floor((double(timeMax)/double(qEASampling))))+1;	//qEA vektor merete; +1 a 0. idopillanat miatt

	Tmax += dT/2.0;

	list<string> infos;

	//*** Vektorok deklaralasa ***//
	vector<double> vecZeros_D(N,0.0);
	vector<int> vecZeros_I(N,0);
	vector<bool> vecFalses(N,false);

	vecZeros_D.resize(eSize,0.0); vector<vector<double> > E(2,vecZeros_D);
	vecZeros_I.resize(mSize,0); vector<vector<int> > m(2,vecZeros_I);
	vecZeros_D.resize(mSize,0.0); vector<vector<double> > m_norm(2,vecZeros_D);
	vecZeros_D.resize(qEAsize,0.0); vector<vector<double> > q_EA(2,vecZeros_D);
	
	vecZeros_D.resize(N,0.0); vecZeros_I.resize(N,0);
	vector<double> hField(N,0.0);
	vector<vector<double> >  C(N,vecZeros_D);
	vector<double> C_CDF(int((N*(N-1))/2),0.0);
	vector<bool> spin(N,false);
	vector<bool> spin_1st(N,false);
	vector<vector<int> >  J(N,vecZeros_I);
	vector<vector<bool> >  J_pos(N,vecFalses);
	vector<vector<bool> >  J_val(N,vecFalses);
	vector<int> s_sum(N,0);
	vector<vector<int> >  s2_sum(N,vecZeros_I);
	vector<double> s_avg(N,0.0);
	vector<vector<double> >  s2_avg(N,vecZeros_D);
	vector<vector<double> >  sisj(N,vecZeros_D);
	vector<double> effField(N,0.0);				//belso ter IDO-atlagolva

	int lineMin; int lineMax;
	double corrMin; double corrMax;

	int mBinNum = N+1;
	vecZeros_I.resize(mBinNum,0); vector<vector<int> >  m_PDDF(2,vecZeros_I);
	vecZeros_D.resize(mBinNum,0.0);  vector<vector<double> >  m_PDDF_norm(2,vecZeros_D);	//A vegen ebben kerul m_PDDF_sum normalva
	vecZeros_D.resize(N,0.0); vecZeros_I.resize(N,0);

	setXaxis_forVecInMat_D(E,0.0,double(eSampling));
	setXaxis_forVecInMat_D(m_norm,0.0,double(mSampling));
	setXaxis_forVecInMat_I(m,0,mSampling);
	setXaxis_forVecInMat_D(q_EA,0.0,double(qEASampling));
	setXaxis_forVecInMat_I(m_PDDF,(-1)*N,2);
	setXaxis_forVecInMat_D(m_PDDF_norm,-1.0,(2.0/double(N)));

	double betha;
	int sweep10percent = int(floor(double(timeMax/10.0)));

	string chosenSpinVal_ = (chosenSpinVal == true) ? "+1" : "-1";

	string strT;
	string strN = "_N"+convert_ItoS(N);
	string frTriang;

	MTRand MTrnd;		//Mersenne Twister rnd szam generator
	
	//J beolvasasa / generalasa
	if (fileJ != "nem")
		J = load_Jint(N,fileJ);
	else
		J = create_rndJint_CG(N,MTrnd);

	frTriang = fr_arany_J(J);		//bruteforce megkeressük a frusztrált/nem ftrusztrált haromszogeket es osszeszamoljuk oket

	//J felbontasa poziciora (J_pos[i][j]: 1, ha van kapcsolat; 0, ha nincs) es ertekre (J_val[i][j]: 1, ha +1-es a kapcsolat; 0, ha -1-es)
	giveVal_and_Pos_matrixes(J,J_pos,J_val);

	//kiindulasi Spin-vektor beolvasasa / generalasa
	if (fileSpin != "nem")
		spin_1st = load_SpinBool_CG(N,fileSpin);
	else
		spin_1st = create_rndSpinBool_CG(N,MTrnd);

	if(fileSpin == "nem") saveVecBool(spin_1st,"spinIn_N" + strN,"");		//kiindulasi vektor kimentese

	cout << "MMC modszer inditasa N=" + convert_ItoS(N) + " eseten" << endl;
	cout << "Kivalasztott spin indexe: " << chosenSpin+1 << " ; Iranya: " << chosenSpinVal_ << endl;
	cout << endl;
		
	for (double t=Tmin; t<=Tmax; t+=dT)
	{
		strT = "_T"+convert_DtoFix2S(t);		//T homerseklet atirasa a file nevekhez
		betha = 1/(kB*t);

	///----------inicializalas-----------///
		cout << "(Inicializalas, elokeszites...)" << endl;
		infos.clear();
		s_sum.clear(); s_sum.resize(N,0);	s2_sum.clear(); s2_sum.resize(N,vecZeros_I);
		effField.clear(); effField.resize(N,0.0);
		for (int i=0; i<mBinNum; i++)
			m_PDDF[1][i] = 0;
	
		//Minten adott N-re es T-re ugyanabbol a spinallapotbol indulunk
		spin = spin_1st;
		spin[chosenSpin] = chosenSpinVal;				//Kivalasztott spin beallitesa iranyra

		//Elso 10%-nak megfelelo lepes kihagyasa a statisztikabol
		for (int st=0; st<sweep10percent; st++)			
			spinFlop_MMC_pin(J_pos,J_val,spin,betha,chosenSpin,MTrnd);
			
	///-------------ido fejlesztes----------------///
		cout << "Rendszer idofejlesztese T = " << t << " eseten" << endl << "\"|\" = 100 000 SpinFlop ~ " + convert_DtoS(100000.0/double(N)) + " sweep: ";
		for (int time=0; time<=timeMax; time++)	
		{
			if (time % 100000 == 0) cout << "|";
			spinFlop_MMC_pin(J_pos,J_val,spin,betha,chosenSpin,MTrnd);	//[12] Metropolis Monte Carlo
			giveSpinSum_MMC(spin,s_sum,s2_sum);					//[2][3] atlagok szamitasahoz osszegzes
			giveEffFieldSum_MMC(effField,J_pos,J_val,spin);		//[13] belso eff.terhez osszagzes

			if (time % eSampling == 0) E[1][int(time/eSampling)] = getEnergy(J_pos,J_val,spin,hField);	//[1] energia szamitasa
			if (time % mSampling == 0) m[1][int(time/mSampling)] = getMagnetism(spin);					//[4] magnesseg szamitasa 
			if (time % qEASampling == 0) q_EA[1][int(time/qEASampling)] = getQEA_sSUM(time,s_sum);		//[16]
			if (time % snapshotFreq == 0)
			{
				saveVecBool(spin,"spinSnapshot" + strN + strT + "_" + convert_Ito3str(int(time/snapshotFreq)),"");
				saveVecBool_meshPBM(spin,snapshotX,snapshotY,spinPixelSize,"spinSnapshot" + strN + strT + "_" + convert_Ito3str(int(time/snapshotFreq)));
			}
		}
			
		cout << endl;
		giveSpinAVG_MMC(timeMax,s_sum,s2_sum,s_avg,s2_avg);
		giveEffFieldAvg_MMC(effField,timeMax);

		for (int i=0; i<mSize; i++)
			m_norm[1][i] = double(m[1][i])/double(N);

		getPDDF_forM(m,m_PDDF);
			
		for (int i=0; i<mBinNum; i++)
			m_PDDF_norm[1][i] = double(m_PDDF[1][i])/double(mSize);			

	///-----korrelációk és távolságuk számítása-----///
		cout << "Adatok kimentese:" << endl;
			
		sisj = getSiSj(s_avg);
		saveMatDoub(E,"E" + strN + strT,"Energy values, sampling freq="+convert_ItoS(eSampling));
		savePDDF_mat(E,1,0.1,0,0,"E" + strN + strT + "_PDDF","Energia Vsz.-suruseg eo.-fv.");
		saveMatDoub(m_norm,"m" + strN + strT, "Magnetism values, sampling freq="+convert_ItoS(mSampling));
		saveMatDoub(m_PDDF_norm,"m" + strN + strT + "_PDDF","Magnesezettseg Vsz.-suruseg eo.-fv.");
		saveMatDoub(q_EA,"q_EA" + strN + strT,"q_EA, Edward-Anderson order-parameter, sampling freq="+convert_ItoS(qEASampling));
		saveMatDoub(sisj,"sisj" + strN + strT,"<si><sj> Matrix");
		saveVecDoub(s_avg,"s_avg" + strN + strT,"<si> vector");
		saveVecDoub(effField,"effField" + strN + strT,"");

		infos.push_back("Rogzitett spin: " + convert_ItoS(chosenSpin+1)); infos.push_back("Iranya: " + chosenSpinVal_);
		infos.push_back("<m>:"); infos.push_back(convert_DtoS(getVectorsAVG(s_avg)));
		infos.push_back("q_EA [futas vegen]:"); infos.push_back(convert_DtoS(getQEA_sAVG(s_avg)));
		infos.push_back(""); infos.push_back("");

		//Korrelacio szamitasa es ezzel kapcsolatos eredmenyek kimentese
		cout << "	- Kapcsolt-Korrelacio szamitasa & adatmentes" << endl;
		C = getConCorrMatrix(s_avg,s2_avg);	//[9]
		C_CDF = getUpperTri_fromSqMatrix_doub(C);

		giveExtrLines_of_matrix(C,lineMin,lineMax,corrMin,corrMax);
		infos.push_back("Kapcsolt-Korrelacioban leggyengebb spin [i]:"); infos.push_back(convert_ItoS(lineMin+1));
		infos.push_back("sum(|c[i]|) ="); infos.push_back(convert_DtoS(corrMin));
		infos.push_back("<s_i> ="); infos.push_back(convert_DtoS(s_avg[lineMin]));
		infos.push_back("effField(i) ="); infos.push_back(convert_DtoS(effField[lineMin]));
		infos.push_back("Kapcsolt-Korrelacioban legerosebb spin [j]:"); infos.push_back(convert_ItoS(lineMax+1));
		infos.push_back("sum(|c[j]|) ="); infos.push_back(convert_DtoS(corrMax));
		infos.push_back("<s_j> ="); infos.push_back(convert_DtoS(s_avg[lineMax]));
		infos.push_back("effField(j) ="); infos.push_back(convert_DtoS(effField[lineMax]));
		infos.push_back(""); infos.push_back("");

		saveMatDoub(C,"C" + strN + strT,"Connected Correlation [C] Matrix");
		saveMatDoub3d(C,"C" + strN + strT + "_3D","Connected Correlation Matrix for 3D plot");
		saveVecDoub(C_CDF,"C" + strN + strT + "_CDF","C Matrix's UpperTriangular for CDF");
		savePDDF_vec(C_CDF,bin_C,-1,1,"C" + strN + strT + "_PDDF","C matrix Vsz.-suruseg eo.-fv.");
		savePDDF_vec_alter(C_CDF,bin_C,-1,1,"C" + strN + strT + "_PDDF_alter","C matrix Vsz.-suruseg eo.-fv.");

		//Kovariancia szamitasa es ezzel kapcsolatos eredmenyek kimentese
		cout << "	- Korrelacio szamitasa & adatmentes" << endl;
		C = getCorrMatrix(s2_avg);		//[10]
		C_CDF = getUpperTri_fromSqMatrix_doub(C);

		giveExtrLines_of_matrix(C,lineMin,lineMax,corrMin,corrMax);
		infos.push_back("Korrelacioban leggyengebb spin [k]:"); infos.push_back(convert_ItoS(lineMin+1));
		infos.push_back("sum(|c[k]|) ="); infos.push_back(convert_DtoS(corrMin));
		infos.push_back("<s_k> ="); infos.push_back(convert_DtoS(s_avg[lineMin]));
		infos.push_back("effField(k) ="); infos.push_back(convert_DtoS(effField[lineMin]));
		infos.push_back("Korrelacioban legerosebb spin [l]:"); infos.push_back(convert_ItoS(lineMax+1));
		infos.push_back("sum(|c[l]|) ="); infos.push_back(convert_DtoS(corrMax));
		infos.push_back("<s_l> ="); infos.push_back(convert_DtoS(s_avg[lineMax]));
		infos.push_back("effField(l) ="); infos.push_back(convert_DtoS(effField[lineMax]));

		saveMatDoub(C,"K" + strN + strT,"Correlation [K] Matrix");
		saveMatDoub3d(C,"K" + strN + strT + "_3D","Correlation Matrix for 3D plot");
		saveVecDoub(C_CDF,"K" + strN + strT + "_CDF","Correlation Matrix's UpperTriangular for CDF");
		savePDDF_vec(C_CDF,bin_C,-1,1,"K" + strN + strT + "_PDDF","Correlation matrix Vsz.-suruseg eo.-fv.");
		savePDDF_vec_alter(C_CDF,bin_C,-1,1,"K" + strN + strT + "_PDDF_alter","K matrix Vsz.-suruseg eo.-fv.");

		saveLS_infos(infos,"infos" + strN + strT);

		cout << endl;
	}
			
	cout << endl;
	cout << "T fuggetlen adatok kimentese" << endl << endl;		

	double fi_J = getJfrustrate_int(J);		//[11] J frusztraltsaganak kiszamitasa
	saveMatInt(J,"J" + strN,"J matrix, N:" + convert_ItoS(N) + "; Frustration Rate: " + convert_DtoS(fi_J) + frTriang);
	saveMatInt3d(J,"J" + strN + "_3D","J Matrix 3D plothoz");
	saveMatInt_meshPGM_forJ(J,"J" + strN,20);

	system("pause");
}

