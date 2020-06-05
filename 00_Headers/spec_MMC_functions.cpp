//#ifndef SAVEFILES_H
//#define SAVEFILES_H

//**********************************************//
//Metropolis - Monte-Carlo specifikus fuggvenyek//
//**********************************************//

using namespace std;

//[12] keplet: energia-különbseget szamol a következo lepeshez kepest es csinal egy spin-flop-ot (leptet), ha szukseges
void spinFlop_MMC(vector<vector<bool> >  J_pos, vector<vector<bool> >  J_val, vector<bool>& s, double betha, vector<double> h, MTRand& MTrnd)
{
	int N = s.size();
	double dE = 0;
	double expBdE;
	
	int i = int(floor(MTrnd() * (N)));	//rnd kivalasztjuk az i-edik spint
	//if (i==N) i=N-1;

	//Energiaváltozás			Az MMC reszig ekvivalens ezzel: for (int j=0; j<N; j++) {dE += 1/N(J[i][j] * s[i] * s[j]);}  dE += h * s[i];
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

	if (s[i] == false)
		dE -= 2*h[i];
	else
		dE += 2*h[i];
									
	//Metropolis-Monte-carlo
	expBdE = exp((-1.0)*betha*dE);
	
	if((dE < 0) || (MTrnd() < expBdE)) 
		s[i] = s[i] ^ true;				//A spin atforditasa (bool negalasa ekvivalens ezzel: s_tmp[i] = (s_tmp[i] == true) ? false : true;)
}

//[2]-[3] keplet: Idobeli sulyozott atlagok szamitasahoz oszegzi a spin-ertekeket (a giveAVGs_MMC fuggvennyel egyutt!)
void giveSpinSum_MMC(vector<bool> s, vector<int>& s_sum, vector<vector<int> >& s2_sum)	
{
	int N = s.size();
	int tmp_si, tmp_sj;
	for (int i=0; i<N; i++)
	{
		tmp_si = (s[i] == true) ? 1 : -1;
		
		s_sum[i] += tmp_si;		
		
		for (int j=0; j<N; j++)
		{
			tmp_sj = (s[j] == true) ? 1 : -1;
			s2_sum[i][j] += tmp_si*tmp_sj;
		}
	}
}

//[2]-[3] keplet: Az idoben oszegzett a spin-ertekekbol (s_sum) visszaadja az Idobeli sulyozott atlagokat
void giveSpinAVG_MMC(int time, vector<int> s_sum, vector<vector<int> > s2_sum, vector<double>& s_avg, vector<vector<double> >& s2_avg)	
{
	int N = s_sum.size();
	for (int i=0; i<N; i++)
	{		
		s_avg[i] = double(s_sum[i])/double(time+1);		//+1 a 0. idopillanat miatt
		
		for (int j=0; j<N; j++)
			s2_avg[i][j] = double(s2_sum[i][j])/double(time+1);	//+1 a 0. idopillanat miatt
	}
}

////[2]-[3] keplet: Idobeli sulyozott atlagok szamitasa
//void giveAVGs_MMC(int time, vector<bool> s, vector<double>& s_avg, vector<vector<double> >& s2_avg)	
//{
//	//megjegyzes: time 0-tol indul, de az egyenletben 1-tol, ezert van az elcsusztatas
//	int N = s.size();
//	double tmp_si; double tmp_sj;
//	for (int i=0; i<N; i++)
//	{
//		tmp_si = (s[i] == true) ? 1.0 : -1.0;
//		
//		s_avg[i]=(s_avg[i]*(double(time)) + tmp_si)/double(time+1);		
//		
//		for (int j=0; j<N; j++)
//		{
//			tmp_sj = (s[j] == true) ? 1.0 : -1.0;
//			s2_avg[i][j]=(s2_avg[i][j]*(time) + tmp_si*tmp_sj)/double(time+1);
//		}
//	}
//}

//[13] keplet: idoatlagolt belso effektiv ter szamitasahoz szummazza a belso tereket
void giveEffFieldSum_MMC(vector<double>& effField, vector<vector<bool> >  J_pos, vector<vector<bool> >  J_val, vector<bool> spin)
{
	int N = spin.size();
	double actualEffField = 0.0;

	for (int i=0; i<N; i++)
	{
		actualEffField = 0.0;			//inicializalas
		for (int j=0; j<N; j++)
		{
			if ((i!=j) && (J_pos[i][j] == true))		//onmagara nem szamolunk teret
			{
				if ((J_val[i][j] ^ spin[j]) == false)		//ekvivalens: effField[i] += J_ij * s_j; 
					actualEffField ++;		
				else
					actualEffField --;
			}
		}
		effField[i] += actualEffField;
	}
}

//[13] keplet: idoatlagolt belso effektiv terhez normalja es atlagolja idore a szummazott tereket
void giveEffFieldAvg_MMC(vector<double>& effField, int time)
{
	int N = effField.size();

	for (int i=0; i<N; i++)
		effField[i] /= double(time+1) * double(N-1);	//+1 a 0. idopillanat miatt
}

//[16] keplet: Edwards-Anderson rendparameter szamitasa az idoatlagolt spinertekek (<si>) alajan
double getQEA_sAVG(vector<double> s_avg)
{
	double qEA = 0.0;
	int N = s_avg.size();

	for (int i=0; i<N; i++)
		qEA += s_avg[i] * s_avg[i];

	qEA /= double(N);

	return qEA;
}

//[16] keplet: Edwards-Anderson rendparameter szamitasa az ido-osszegzett spinertekek (s_sum) alajan
double getQEA_sSUM(int actualTime, vector<int> s_sum)
{
	double qEA = 0.0;

	int N = s_sum.size();
	double tmp = 0.0;

	for (int i=0; i<N; i++)
	{
		tmp = double(s_sum[i])/double(actualTime+1);
		qEA += tmp * tmp;
	}

	qEA /= double(N);

	return qEA;
}

//[17] keplet: q_(alpha-betha) rendparameter szamitasa
double getQAB(vector<bool> spin1, vector<bool> spin2)
{
	double qAB = 0.0;
	int N = spin1.size();

	for (int i=0; i<N; i++)
		qAB += ((spin1[i] ^ spin2[i]) == false) ? +1.0 : -1.0;

	qAB /= double(N);

	return qAB;
}


//#endif