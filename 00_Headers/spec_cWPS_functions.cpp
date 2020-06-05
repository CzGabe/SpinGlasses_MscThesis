//#ifndef SAVEFILES_H
//#define SAVEFILES_H

//*********************************************//
//count Whole Phase Space specifikus fuggvenyek//
//*********************************************//

using namespace std;

//"spin-lepteto" (gyakorlatilag ugy mukodik, mint egy binaris szamlalo)
void giveNextSpin (vector<bool> & s)
{
	bool check = false;
	int i = 0;

	while (check == false)
	{
		if (s[i] == false)
		{
			check = true;
			s[i] = true;
		}
		else	//s[i] = true
		{
			s[i] = false;
		}
		i++;
	}
}

//[6]-[7] kepletek: Teljes leszamlalasnal Boltzman-faktorokkal szamoljuk az atlagokat (az allapotok fuggetlenek egymastol, nincs "idofejlodes")
void giveSpinAVG_cWPS(vector<bool> s, vector<double> E, double betha, int stateNum, vector<double>& s_avg, vector<vector<double> > & s2_avg)
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

		if (st < stateNum-1) giveNextSpin(s);
	}

	//---Normalas---//
	for (int i=0; i<N; i++)
	{
		s_avg[i] /= Z;

		for (int j=0; j<N; j++)
			s2_avg[i][j] /= Z;
	}
}

//[8] keplet: Magnesseg szamitasa - kicsit trukkos: egybol az eloszlast szamolja.
vector<vector<double> > giveMagn_cWPS(vector<bool> s, vector<double> E, double betha, int stateNum)
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

		if (st < stateNum-1) giveNextSpin(s);
	}

	//---Normalas---//
	for (int st=0; st<mNum; st++)
		magn[1][st] /= Z;
		
	return magn;
}

//[13] Belos effektiv ter szamitasa
void giveEffFielf_cWPS(vector<double>& effField, int stateNum, vector<vector<bool> >  J_pos, vector<vector<bool> >  J_val, vector<bool> spin)
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
		if (st < stateNum-1) giveNextSpin(spin);
	}

	for (int i=0; i<N; i++)
		effField[i] /= (double(N-1) * double(stateNum));
}

//#endif