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

void saveMatBool_alt(vector<vector<bool> >  M, string fileName, string firstLine)
{
	firstLine = (firstLine == "") ? "# True = 1, False = 0" : "# True = 1, False = 0" + firstLine;
	fileName += ".dat";
	int row = M[0].size(); /*N*/ int col = M.size(); /*st*/
	int tmp;
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << firstLine << endl;

	for(int i=0; i<col; i++)
	{
		for(int j=0; j<row; j++)
		{
			tmp = (M[i][j] == true) ? 1 : 0;
			save << tmp << '\t';
		}

		save << endl;
	}
	save.close();	
}


//kiemeljuk kulon vektorba a nem rogzitett spineket, leptetjuk oket, majd visszatesszuk a helyukre...
void giveNextSpin_Pin (vector<bool> & spin, int chosenSpin)
{
	int N = spin.size();
	vector<bool> s(N-1,false);
	
	int j=0;
	for (int i=0; i<N; i++)
	{
		if (i==chosenSpin) i++;
		s[j] = spin[i];
		j++;
	}
	
	bool check = false;
	int k = 0;

	while (check == false)
	{
		if (s[k] == false)
		{
			check = true;
			s[k] = true;
		}
		else	//s[i] = true
		{
			s[k] = false;
		}
		k++;
	}

	j=0;
	for (int i=0; i<N; i++)
	{
		if (i==chosenSpin) i++;
		spin[i] = s[j];
		j++;
	}
}

void giveSpinAVG_cWPS_PIN(vector<bool> s, vector<double> E, double betha, int stateNum, int chosenSpin, vector<double>& s_avg, vector<vector<double> > & s2_avg)
{
	int N = s.size();
	double Z = 0.0;
	double expBE = 0.0;
	
	for (int st=0; st<stateNum; st++)
	{
		expBE = exp((-1) * betha * E[st]);
		Z += expBE;		//[7]

		for (int i=0; i<N; i++)
		{
			s_avg[i] += (s[i] == true) ? expBE : -expBE;

			for (int j=0; j<N; j++)
				s2_avg[i][j] += ((s[i] ^ s[j]) != true) ? expBE : -expBE;		//[6]; a feltetel ekvivalens ezzel: (s_i * s_j == 1)
		}

		if (st < stateNum-1) giveNextSpin_Pin(s,chosenSpin);
	}

	//---Normalas---//
	for (int i=0; i<N; i++)
	{
		s_avg[i] /= Z;

		for (int j=0; j<N; j++)
			s2_avg[i][j] /= Z;
	}
}

vector<vector<double> > giveMagn_cWPS_PIN(vector<bool> s, vector<double> E, double betha, int stateNum, int chosenSpin)
{
	int N = s.size();
	double Z = 0.0;
	double expBE = 0.0;
	double avgMagn = 0.0;
	int mNum = N + 1;		//lehetseges m ertekek szama
	int mNumAll = 2*N + 1;
	int mVal = 0;			//aktualis m ertek

	vector<double> vecZeros(mNum,0.0);		
	vector<vector<double> >  magn(2,vecZeros);

	setXaxis_forVecInMat_D(magn,-1.0,(2.0/double(N)));
	
	for (int st=0; st<stateNum; st++)
	{
		mVal = 0;

		expBE = exp((-1) * betha * E[st]);
		Z += expBE;

		for (int i=0; i<N; i++)
			mVal += (s[i] == true) ? 1 : -1;

		magn[1][int((mVal+N)/2)] += expBE;			//a masodik oszlopba szummazzuk a Boltzman-faktorokat

		if (st < stateNum-1) giveNextSpin_Pin(s,chosenSpin);
	}

	//---Normalas---//
	for (int st=0; st<mNum; st++)
		magn[1][st] /= Z;
		
	return magn;
}

void giveEffFielf_cWPS_PIN(vector<double>& effField, int stateNum, int chosenSpin, vector<vector<bool> >  J_pos, vector<vector<bool> >  J_val, vector<bool> spin)
{
	int N = spin.size();

	for (int st=0; st<stateNum; st++)
	{
		for (int i=0; i<N; i++)
		{
			for (int j=0; j<N; j++)
			{
				if ((i!=j) && (J_pos[i][j] == true))		//onmagara nem szamolunk teret
				{
					if ((J_val[i][j] ^ spin[j]) == false)		//ekvivalens: effField[i] += J_ij * s_j; 
						effField[i] ++;		
					else
						effField[i] --;
				}
			}
		}
		if (st < stateNum-1) giveNextSpin_Pin(spin,chosenSpin);
	}

	for (int i=0; i<N; i++)
		effField[i] /= (double(N-1) * double(stateNum));
}

int main()
{
	///--------------Parameterek beolvasasa-----------------///
	cout << "Parameterek beolvasasa" << endl << endl;
	int N;
	string fileP = "03_SG_Param_cWPS_PinSpin.dat";
	string fileJ = loadParameter_string(fileP,3);
	string fileH = "";
	
	if (fileJ == "nem") 
		N = loadParameter_int(fileP,6);
	else
		N = load_N_fromJfile(fileJ);

	double Tmin = loadParameter_double(fileP,9);
	double Tmax = loadParameter_double(fileP,12);
	double dT = loadParameter_double(fileP,15);
	int chosenSpin = loadParameter_int(fileP,18)-1;
	string chosenSpinIndex = loadParameter_string(fileP,18);
	bool chosenSpinVal = loadParameter_bool(fileP,21,"+1");
	double kB = loadParameter_double(fileP,24);
	double bin_C = loadParameter_double(fileP,27);
	
	int stateNum = int(pow(2.0,N)/2.0);

	MTRand MTrnd;		//Mersenne Twister rnd szam generator

	Tmax += dT/2.0;

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

	vecZeros.resize(N,0.0);
	int lineMin;
	int lineMax;
	double corrMin; double corrMax;
	double fi_J;
	
	double betha;

	string chosenSpinVal_ = (chosenSpinVal == true) ? "+1" : "-1";

	string strN = "_N"+convert_ItoS(N);
	string strT;
	string frTriang;

	//J beolvasasa / generalasa
	if (fileJ != "nem")
		J = load_Jint(N,fileJ);
	else
		J = create_rndJint_CG(N,MTrnd);

	frTriang = fr_arany_J(J);		//bruteforce megkeressük a frusztrált/nem ftrusztrált haromszogeket es osszeszamoljuk oket

	//J felbontasa poziciora (J_pos[i][j]: 1, ha van kapcsolat; 0, ha nincs) es ertekre (J_val[i][j]: 1, ha +1-es a kapcsolat; 0, ha -1-es)
	giveVal_and_Pos_matrixes(J,J_pos,J_val);
	
	cout << "Teljes Leszamlalas N = " << N << " eseten" << endl << endl;
	cout << "Kivalasztott spin indexe: " << chosenSpin+1 << " ; Iranya: " << chosenSpinVal_ << endl;

	///--------------Energiak Szamitasa-------------------///
	cout << "Energia szamitasa" << endl << endl;
			

	spin.clear(); spin.resize(N,false);
	spin[chosenSpin] = chosenSpinVal;				//Kivalasztott spin beallitesa iranyra

	for (int st=0; st<stateNum; st++)					//st -> az allapotokon megy vegig
	{
		E[st] = getEnergy(J_pos,J_val,spin,hField);		//energia számítása (Az allapotokhoz tartozo energia szamitasa T fuggetlen)
		if (st < stateNum-1) giveNextSpin_Pin(spin,chosenSpin);
	}

	///DELETE///
	vecFalses.resize(N,false);
	vector<vector<bool> >  spinStates(stateNum,vecFalses);
	spin.clear(); spin.resize(N,false);
	spin[chosenSpin] = chosenSpinVal;				//Kivalasztott spin beallitesa iranyra

	for (int st=0; st<stateNum; st++)					//st -> az allapotokon megy vegig
	{
		spinStates[st] = spin;
		if (st < stateNum-1) giveNextSpin_Pin(spin,chosenSpin);
	}

	saveMatBool_alt(spinStates,"STATES","");
	///DELETE///

	///-------------Kulonbozo homersekletekre szamitas----///
	for (double t=Tmin; t<=Tmax; t+=dT)
	{
		strT = "_T"+convert_DtoFix2S(t);
		cout << "Tmin=" << Tmin << "	T=" << t << "	Tmax=" << Tmax << endl;
				
		infos.clear();

		betha = 1/(kB*t);

		//0. lepesek: minden t-nel "friss" s(2)_avg és spin vektorokkal indulunk
		spin.clear(); spin.resize(N,false);				//ezt is, mert a "giveAVGs" -ben fogjuk leptetni
		spin[chosenSpin] = chosenSpinVal;				//Kivalasztott spin beallitesa iranyra
		s_avg.clear(); s_avg.resize(N,0.0);
		s2_avg.clear(); s2_avg.resize(N,vecZeros);
		effField.clear(); effField.resize(N,0.0);

		cout << "Fazister vegigjarasa:" << endl <<"   - Atlagok szamitasa" << endl;
		giveSpinAVG_cWPS_PIN(spin,E,betha,stateNum,chosenSpin,s_avg,s2_avg);	//[5][6][7]: átlagok számítása korrelációkhoz

		cout << "   - Magnesezettseg szamitasa" << endl;
		spin.clear(); spin.resize(N,false);					//itt is ujra leptetjuk a spineket
		spin[chosenSpin] = chosenSpinVal;					//Kivalasztott spin beallitesa iranyra
		m = giveMagn_cWPS_PIN(spin,E,betha,stateNum,chosenSpin);			//[8]: magnesseg szamitasa

		cout << "   - Belso effektiv ter szamitasa" << endl;
		spin.clear(); spin.resize(N,false);
		spin[chosenSpin] = chosenSpinVal;					//Kivalasztott spin beallitesa iranyra
		giveEffFielf_cWPS_PIN(effField,stateNum,chosenSpin,J_pos,J_val,spin);	//[13]

		cout << "Adatok kimentese" << endl;

		saveMatDoub(m,"m" + strN + strT  + "_PDDF","Magnesezettseg Vsz.-suruseg eo.-fv.");

		cout << "	- Kapcsolt-Korrelacio szamitasa & adatmentes" << endl;
		C = getConCorrMatrix(s_avg,s2_avg);		//[9]
		C_CDF = getUpperTri_fromSqMatrix_doub(C);
		
		infos.push_back("Rogzitett spin: " + convert_ItoS(chosenSpin+1)); infos.push_back("Iranya: " + chosenSpinVal_);

		giveExtrLines_of_matrix(C,lineMin,lineMax,corrMin,corrMax);
		infos.push_back("Kapcsolt-Korrelacioban leggyengebb spin [i]:"); infos.push_back(convert_ItoS(lineMin+1));
		infos.push_back("sum(|c[i]|) ="); infos.push_back(convert_DtoS(corrMin));
		infos.push_back("effField(i) ="); infos.push_back(convert_DtoS(effField[lineMin]));
		infos.push_back("Kapcsolt-Korrelacioban legerosebb spin [j]:"); infos.push_back(convert_ItoS(lineMax+1));
		infos.push_back("sum(|c[j]|) ="); infos.push_back(convert_DtoS(corrMax));
		infos.push_back("effField(j) ="); infos.push_back(convert_DtoS(effField[lineMax]));
		infos.push_back(""); infos.push_back("");
			
		saveMatDoub(C,"C" + strN + strT ,"Connected Correlation [C] Matrix");
		saveMatDoub3d(C,"C" + strN + strT  + "_3D","Connected Correlation Matrix for 3D plot");
		saveVecDoub(C_CDF,"C" + strN + strT  + "_CDF","C Matrix's UpperTriangular for CDF");
		savePDDF_vec(C_CDF,bin_C,-1,1,"C" + strN + strT  + "_PDDF","C matrix Vsz.-suruseg eo.-fv.");
		savePDDF_vec_alter(C_CDF,bin_C,-1,1,"C" + strN + strT  + "_PDDF_alter","C matrix Vsz.-suruseg eo.-fv.");

		cout << "	- Korrelacio szamitasa & adatmentes" << endl;
		C = getCorrMatrix(s2_avg);		//[9]
		C_CDF = getUpperTri_fromSqMatrix_doub(C);
			
		giveExtrLines_of_matrix(C,lineMin,lineMax,corrMin,corrMax);
		infos.push_back("Korrelacioban leggyengebb spin [i]:"); infos.push_back(convert_ItoS(lineMin+1));
		infos.push_back("sum(|c[i]|) ="); infos.push_back(convert_DtoS(corrMin));
		infos.push_back("effField(i) ="); infos.push_back(convert_DtoS(effField[lineMin]));
		infos.push_back("Korrelacioban legerosebb spin [j]:"); infos.push_back(convert_ItoS(lineMax+1));
		infos.push_back("sum(|c[j]|) ="); infos.push_back(convert_DtoS(corrMax));
		infos.push_back("effField(j) ="); infos.push_back(convert_DtoS(effField[lineMax]));
		infos.push_back(""); infos.push_back("");
			
		saveMatDoub(C,"K" + strN + strT ,"Correlation [K] Matrix");						
		saveMatDoub3d(C,"K" + strN + strT  +"_3D","Correlation Matrix for 3D plot");
		saveVecDoub(C_CDF,"K" + strN + strT  + "_CDF","Correlation Matrix's UpperTriangular for CDF");
		savePDDF_vec(C_CDF,bin_C,-1,1,"K" + strN + strT  + "_PDDF","Vsz.-suruseg eo.-fv.");
		savePDDF_vec_alter(C_CDF,bin_C,-1,1,"K" + strN + strT  + "_PDDF_alter","Vsz.-suruseg eo.-fv.");
			
		saveLS_infos(infos,"infos" + strN + strT );
		
		cout << endl;
		cout << "Energia kimentese" << endl << endl;
		saveVecDoub(E,"E" + strN ,"E");
		savePDDF_vec(E,0.1,0,0,"E" + strN  + "_PDDF","Energia Vsz.-suruseg eo.-fv.");
	}

	cout << "T fuggetlen adatok kimentese" << endl << endl;
	fi_J = getJfrustrate_int(J);		//[11] J frusztraltsaganak kiszamitasa
	saveMatInt(J,"J" + strN,"J matrix, N:" + convert_ItoS(N) + "; Frustration Rate: " + convert_DtoS(fi_J) + frTriang);
	saveMatInt3d(J,"J" + strN + "_3D","J Matrix 3D plothoz");

	system("pause");
}
