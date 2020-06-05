//#include "stdafx.h"
#include <fstream>
#include <list>
#include <vector>
#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
//#include <math.h>
//***Mersenne Twister Random Number Generator
#include "MersenneTwister.h"
//***Sajat fileok-ek
#include "saveFiles.cpp"
#include "readFromFiles.cpp"
#include "variousFunctions.cpp"
#include "equationFunctions.cpp"
#include "spec_cWPS_Functions.cpp"

using namespace std;

int main()
{
	///--------------Parameterek beolvasasa-----------------///
	cout << "Parameterek beolvasasa" << endl << endl;
	int N;
	string fileP = "01_SG_Param_cWPS.dat";
	string fileJ = loadParameter_string(fileP,3);
	string fileH = "";
	
	if (fileJ == "nem") 
		N = loadParameter_int(fileP,6);
	else
		N = load_N_fromJfile(fileJ);

	double Tmin = loadParameter_double(fileP,9);
	double Tmax = loadParameter_double(fileP,12);
	double dT = loadParameter_double(fileP,15);
	bool homogenField = loadParameter_bool(fileP,18,"nem");
	if (homogenField == false) fileH = loadParameter_string(fileP,18);
	double Hmin = loadParameter_double(fileP,21);
	double Hmax = loadParameter_double(fileP,24);
	double dH = loadParameter_double(fileP,27);
	double kB = loadParameter_double(fileP,30);
	double bin_C = loadParameter_double(fileP,33);
	
	int stateNum = int(pow(2.0,N));

	MTRand MTrnd;		//Mersenne Twister rnd szam generator

	Tmax += dT/2.0; Hmax += dH/2.0;

	list<string> infos;

	//*** Vektorok deklaralasa ***//
	vector<double> vecZeros(N,0.0);
	vector<int> vecZeros_I(N,0);
	vector<bool> vecFalses(N,false);

	vector<double> E(stateNum,0.0);
	vector<double> hField(N,0.0);
	vector<vector<double> >  C(N,vecZeros);
	vector<double> C_CDF(int((N*(N-1))/2),0.0);
	vector<bool> spin(N,false);
	vector<vector<int> >  J(N,vecZeros_I);
	vector<vector<bool> >  J_pos(N,vecFalses);
	vector<vector<bool> >  J_val(N,vecFalses);
	vector<double> s_avg(N,0.0);
	vector<vector<double> >  s2_avg(N,vecZeros);
	vector<vector<double> >  m;
	vector<double> effField(N,0.0);		//allapotokra kiatlagolt

	int H_num = int(((Hmax-Hmin)/dH)+1);
	int T_num = int(((Tmax-Tmin)/dT)+1);

	vector<vector<double> >  C_0;		//H=0 teres eloszlasok mentese d_C0Ck-hoz
	vector<vector<double> >  K_0;
	vecZeros.resize(H_num,0.0);
	vector<vector<double> >  d_C0Ck(T_num,vecZeros);
	vector<vector<double> >  d_CkJ(T_num,vecZeros);
	vector<vector<double> >  d_K0Kk(T_num,vecZeros);
	vector<vector<double> >  d_KkJ(T_num,vecZeros);
	vector<vector<double> >  C_Var(T_num,vecZeros);		//szorasok
	vector<vector<double> >  K_Var(T_num,vecZeros);
	vecZeros.resize(N,0.0);
	int Tindex = -1;
	int Hindex = -1;
	int lineMin;
	int lineMax;
	double corrMin; double corrMax;
	double fi_J;
	
	double betha;

	string strN = "_N"+convert_ItoS(N);
	string strT;
	string strH;
	string frTriang;

	//J beolvasasa / generalasa
	if (fileJ != "nem")
		J = load_Jint(N,fileJ);
	else
		J = create_rndJint_CG(N,MTrnd);

	frTriang = fr_arany_J(J);		//bruteforce megkeressük a frusztrált/nem ftrusztrált haromszogeket es osszeszamoljuk oket

	//J felbontasa poziciora (J_pos[i][j]: 1, ha van kapcsolat; 0, ha nincs) es ertekre (J_val[i][j]: 1, ha +1-es a kapcsolat; 0, ha -1-es)
	giveVal_and_Pos_matrixes(J,J_pos,J_val);

	if (homogenField == false) 
	{
		hField = load_inhomogenH(N, fileH);
		Hmin = 0; Hmax = 1; dH = 2;		//A H-ciklus csak 1x fusson le
	}
	
	for (double h=Hmin; h<=Hmax; h+=dH)
	{
		Hindex++; Tindex = -1;
		if (homogenField == true)					//H ter beallitasa
		{
			strH = "_H"+convert_DtoFix2S(h);		//H ter atirasa a file nevekhez
			for (int i=0; i<N; i++)					//homogen H ter
				hField[i] = h;		
		}
		else
		{
			strH = "_inhomH";		//H ter atirasa a file nevekhez
		}

		cout << "Teljes Leszamlalas N = " << N << " es H = " << h << " eseten" << endl << endl;

		///--------------Energiak Szamitasa-------------------///
		cout << "Energia szamitasa" << endl << endl;
			
		spin.clear(); spin.resize(N,false);
		E[0] = getEnergy(J_pos,J_val,spin,hField);

		for (int st=1; st<stateNum; st++)					//st -> az allapotokon megy vegig
		{
			giveNextSpin(spin);								//s leptetese (alapbol csupa 0 (false))
			E[st] = getEnergy(J_pos,J_val,spin,hField);		//energia számítása (Az allapotokhoz tartozo energia szamitasa T fuggetlen)
		}

		///-------------Kulonbozo homersekletekre szamitas----///
		for (double t=Tmin; t<=Tmax; t+=dT)
		{
			Tindex++;
			strT = "_T"+convert_DtoFix2S(t);
			cout << "Tmin=" << Tmin << "	T=" << t << "	Tmax=" << Tmax << '\t' << "(N=" << N << ";H=" << h << ")" << endl;
				
			infos.clear();

			betha = 1/(kB*t);

			//0. lepesek: minden t-nel "friss" s(2)_avg és spin vektorokkal indulunk
			spin.clear(); spin.resize(N,false);				//ezt is, mert a "giveAVGs" -ben fogjuk leptetni
			s_avg.clear(); s_avg.resize(N,0.0);
			s2_avg.clear(); s2_avg.resize(N,vecZeros);
			effField.clear(); effField.resize(N,0.0);

			cout << "Fazister vegigjarasa:" << endl <<"   - Atlagok szamitasa" << endl;
			giveSpinAVG_cWPS(spin,E,betha,stateNum,s_avg,s2_avg);	//[5][6][7]: átlagok számítása korrelációkhoz

			cout << "   - Magnesezettseg szamitasa" << endl;
			spin.clear(); spin.resize(N,false);					//itt is ujra leptetjuk a spineket
			m = giveMagn_cWPS(spin,E,betha,stateNum);			//[8]: magnesseg szemitasa

			cout << "   - Belso effektiv ter szamitasa" << endl;
			spin.clear(); spin.resize(N,false);
			giveEffFielf_cWPS(effField,stateNum,J_pos,J_val,spin);	//[13]

			cout << "Adatok kimentese" << endl;

			saveMatDoub(m,"m" + strN + strT + strH + "_PDDF","Magnesezettseg Vsz.-suruseg eo.-fv.");

			cout << "	- Kapcsolt-Korrelacio szamitasa & adatmentes" << endl;
			C = getConCorrMatrix(s_avg,s2_avg);		//[9]
			C_CDF = getUpperTri_fromSqMatrix_doub(C);

			if (h == Hmin) C_0 = C;
			d_C0Ck[Tindex][Hindex] = get_dC0Ck(C_0,C);		//[14]
			d_CkJ[Tindex][Hindex] = get_dCkJ(C,J);			//[15]
			C_Var[Tindex][Hindex] = get_SqrtVar_ofVec(C_CDF);
		
			giveExtrLines_of_matrix(C,lineMin,lineMax,corrMin,corrMax);
			infos.push_back("Kapcsolt-Korrelacioban leggyengebb spin [i]:"); infos.push_back(convert_ItoS(lineMin+1));
			infos.push_back("sum(|c[i]|) ="); infos.push_back(convert_DtoS(corrMin));
			infos.push_back("effField(i) ="); infos.push_back(convert_DtoS(effField[lineMin]));
			infos.push_back("Kapcsolt-Korrelacioban legerosebb spin [j]:"); infos.push_back(convert_ItoS(lineMax+1));
			infos.push_back("sum(|c[j]|) ="); infos.push_back(convert_DtoS(corrMax));
			infos.push_back("effField(j) ="); infos.push_back(convert_DtoS(effField[lineMax]));
			infos.push_back(""); infos.push_back("");
			
			saveMatDoub(C,"C" + strN + strT + strH,"Connected Correlation [C] Matrix");
			saveMatDoub3d(C,"C" + strN + strT + strH + "_3D","Connected Correlation Matrix for 3D plot");
			saveVecDoub(C_CDF,"C" + strN + strT + strH + "_CDF","C Matrix's UpperTriangular for CDF");
			savePDDF_vec(C_CDF,bin_C,-1,1,"C" + strN + strT + strH + "_PDDF","C matrix Vsz.-suruseg eo.-fv.");
			savePDDF_vec_alter(C_CDF,bin_C,-1,1,"C" + strN + strT + strH + "_PDDF_alter","C matrix Vsz.-suruseg eo.-fv.");

			cout << "	- Korrelacio szamitasa & adatmentes" << endl;
			C = getCorrMatrix(s2_avg);		//[9]
			C_CDF = getUpperTri_fromSqMatrix_doub(C);

			if (h == Hmin) K_0 = C;
			d_K0Kk[Tindex][Hindex] = get_dC0Ck(K_0,C);		//[14]
			d_KkJ[Tindex][Hindex] = get_dCkJ(C,J);			//[15]
			K_Var[Tindex][Hindex] = get_SqrtVar_ofVec(C_CDF);
			
			giveExtrLines_of_matrix(C,lineMin,lineMax,corrMin,corrMax);
			infos.push_back("Korrelacioban leggyengebb spin [i]:"); infos.push_back(convert_ItoS(lineMin+1));
			infos.push_back("sum(|c[i]|) ="); infos.push_back(convert_DtoS(corrMin));
			infos.push_back("effField(i) ="); infos.push_back(convert_DtoS(effField[lineMin]));
			infos.push_back("Korrelacioban legerosebb spin [j]:"); infos.push_back(convert_ItoS(lineMax+1));
			infos.push_back("sum(|c[j]|) ="); infos.push_back(convert_DtoS(corrMax));
			infos.push_back("effField(j) ="); infos.push_back(convert_DtoS(effField[lineMax]));
			infos.push_back(""); infos.push_back("");
			
			saveMatDoub(C,"K" + strN + strT + strH,"Correlation [K] Matrix");						
			saveMatDoub3d(C,"K" + strN + strT + strH +"_3D","Correlation Matrix for 3D plot");
			saveVecDoub(C_CDF,"K" + strN + strT + strH + "_CDF","Correlation Matrix's UpperTriangular for CDF");
			savePDDF_vec(C_CDF,bin_C,-1,1,"K" + strN + strT + strH + "_PDDF","Vsz.-suruseg eo.-fv.");
			savePDDF_vec_alter(C_CDF,bin_C,-1,1,"K" + strN + strT + strH + "_PDDF_alter","Vsz.-suruseg eo.-fv.");
			
			saveLS_infos(infos,"infos" + strN + strT + strH);

			cout << endl;
		}
			
		cout << endl;
		cout << "H-hoz tartozo energia kimentese" << endl << endl;
		saveVecDoub(E,"E" + strN + strH,"E");
		savePDDF_vec(E,0.1,0,0,"E" + strN + strH + "_PDDF","Energia Vsz.-suruseg eo.-fv.");
	}

	cout << "(T,H) fuggetlen adatok kimentese" << endl << endl;
	fi_J = getJfrustrate_int(J);		//[11] J frusztraltsaganak kiszamitasa
	saveMatInt(J,"J" + strN,"J matrix, N:" + convert_ItoS(N) + "; Frustration Rate: " + convert_DtoS(fi_J) + frTriang);
	saveMatInt3d(J,"J" + strN + "_3D","J Matrix 3D plothoz");
	saveMatInt_meshPGM_forJ(J,"J" + strN,20);

	saveMatDoub(d_C0Ck,"dC0Ck" + strN,"d_C0Ck; sor ~ T tengely; oszlop ~ H tengely");
	saveMatDoub3d_wAxes(d_C0Ck,Tmin,Tmax,dT,Hmin,Hmax,dH,"dC0Ck" + strN + "_3D","x tengely = T; y tengely = H");
	saveMatDoub(d_CkJ,"dCkJ" + strN, "d_CkJ; sor ~ T tengely; oszlop ~ H tengely");
	saveMatDoub3d_wAxes(d_CkJ,Tmin,Tmax,dT,Hmin,Hmax,dH,"dCkJ" + strN + "_3D","x tengely = T; y tengely = H");
	saveMatDoub(C_Var,"C_Var" + strN, "Kapcsolt-Korrelacios eloszlasok szorasa; sor ~ T tengely; oszlop ~ H tengely");
	saveMatDoub3d_wAxes(C_Var,Tmin,Tmax,dT,Hmin,Hmax,dH,"C_Var" + strN + "_3D","A szorasok szorasa= "+ convert_DtoS(get_SqrtVar_ofMat(C_Var)) +"; x tengely = T; y tengely = H");

	saveMatDoub(d_K0Kk,"dK0Kk" + strN,"d_K0Kk; sor ~ T tengely; oszlop ~ H tengely");
	saveMatDoub3d_wAxes(d_K0Kk,Tmin,Tmax,dT,Hmin,Hmax,dH,"dK0Kk" + strN + "_3D","x tengely = T; y tengely = H");
	saveMatDoub(d_KkJ,"dKkJ" + strN, "d_KkJ; sor ~ T tengely; oszlop ~ H tengely");
	saveMatDoub3d_wAxes(d_KkJ,Tmin,Tmax,dT,Hmin,Hmax,dH,"dKkJ" + strN + "_3D","x tengely = T; y tengely = H");
	saveMatDoub(K_Var,"K_Var" + strN, "Korrelacios eloszlasok szorasa; sor ~ T tengely; oszlop ~ H tengely");
	saveMatDoub3d_wAxes(K_Var,Tmin,Tmax,dT,Hmin,Hmax,dH,"K_Var" + strN + "_3D","A szorasok szorasa= "+ convert_DtoS(get_SqrtVar_ofMat(K_Var)) +"; x tengely = T; y tengely = H");

	system("pause");
}
