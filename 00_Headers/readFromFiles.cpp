//#ifndef READFROMFILES_H
//#define READFROMFILES_H

using namespace std;

string loadParameter_string(string fileName, int lineNum)
{
	string answer;

	char *fileName_char = (char*)fileName.c_str();
	ifstream paramFile(fileName_char);

	if (paramFile.is_open())
	{
		for(int i=0; i<lineNum; i++)
			getline (paramFile,answer);
	}
	else
	{
		cout << fileName << " nem megnyithato"; 
		system("pause");
	}

	paramFile.close();

	return answer;
}

double loadParameter_double(string fileName, int lineNum)
{
	string tmp;
	double answer;

	char *fileName_char = (char*)fileName.c_str();
	ifstream paramFile(fileName_char);

	if (paramFile.is_open())
	{
		for(int i=0; i<lineNum; i++)
			getline (paramFile,tmp);

		answer = atof(tmp.c_str());
	}
	else
	{
		cout << fileName << " nem megnyithato"; 
		system("pause");
	}

	paramFile.close();

	return answer;
}

int loadParameter_int(string fileName, int lineNum)
{
	string tmp;
	int answer;

	char *fileName_char = (char*)fileName.c_str();
	ifstream paramFile(fileName_char);

	if (paramFile.is_open())
	{
		for(int i=0; i<lineNum; i++)
			getline (paramFile,tmp);

		answer = atoi(tmp.c_str());
	}
	else
	{
		cout << fileName << " nem megnyithato"; 
		system("pause");
	}

	paramFile.close();

	return answer;
}

bool loadParameter_bool(string fileName, int lineNum, string compare)
{
	string tmp;
	bool answer;
	
	char *fileName_char = (char*)fileName.c_str();
	ifstream paramFile(fileName_char);
	
	if (paramFile.is_open())
	{
		for(int i=0; i<lineNum; i++)
			getline (paramFile,tmp);

		answer = (tmp == compare) ? true : false;
	}
	else
	{
		cout << fileName << " nem megnyithato"; 
		system("pause");
	}

	paramFile.close();

	return answer;
}

int load_N_fromJfile(string fileName)
{
	string tmp = ""; string val = "";
	int N;

	char *fileName_char = (char*)fileName.c_str();
	ifstream Jfile(fileName_char);

	if (Jfile.is_open())
	{
		while (tmp != ":") 
			tmp = Jfile.get();

		while (tmp != ";")
		{
			tmp = Jfile.get();
			if (tmp != ";") val += tmp;
		}

		N = atoi(val.c_str());
	}
	else
	{
		cout << "Gond van a " << fileName << " file-lal. Kerem kapcsolja ki!" << endl;
		system("pause");
	}
	
	return N;
}

void load_XandY_fromJfile(string fileName, int& X, int& Y)
{
	string tmp = ""; string val = "";

	char *fileName_char = (char*)fileName.c_str();
	ifstream Jfile(fileName_char);

	if (Jfile.is_open())
	{
		while (tmp != "=") 
			tmp = Jfile.get();

		while (tmp != "|")
		{
			tmp = Jfile.get();
			if (tmp != "|") val += tmp;
		}

		X = atoi(val.c_str());

		tmp = ""; val = "";

		while (tmp != "=")
			tmp = Jfile.get();

		while (tmp != "|")
		{
			tmp = Jfile.get();
			if (tmp != "|") val += tmp;
		}

		Y = atoi(val.c_str());
	}
	else
	{
		cout << "Gond van a " << fileName << " file-lal. Kerem kapcsolja ki!" << endl;
		system("pause");
	}
}

//J beolvasasa
vector<vector<bool> >  load_Jbool_CG(int N, string fileName)
{
	vector<bool> vecFalse(N,false);
	vector<vector<bool> >  J(N,vecFalse);
	
	char *fileName_char = (char*)fileName.c_str();

	FILE *Jfile;
	Jfile = fopen(fileName_char,"r");
	if (Jfile == NULL) 
	{
		cout << "Gond van a " << fileName << " file-lal. Kerem kapcsolja ki!" << endl;
		system("pause");
	}
	
	while((fgetc(Jfile))!='\n'){}		//Elso sor atugrasa

	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			J[i][j] = ((fgetc(Jfile)) == '1') ? true : false;
			fgetc(Jfile);
		}
	}

	fclose(Jfile);
	return J;
}

//J beolvasasa
vector<vector<int> >  load_Jint(int N, string fileName)
{
	vector<int> vecZeros_I(N,0);
	vector<vector<int> >  J(N,vecZeros_I);
	string tmp = "";
	string val = "";

	char *fileName_char = (char*)fileName.c_str();
	
	ifstream Jfile(fileName_char);
	if (Jfile.is_open())
	{
		while (tmp != "\n") {tmp = Jfile.get();}	//1. sor atugrasa

		for (int i=0; i<N; i++)
		{
			for (int j=0; j<N; j++)
			{
				tmp ="";
				val = "";
				while ((tmp != "\t") && (tmp != "\n"))
				{
					val += tmp;
					tmp = Jfile.get();
					if ((j==0) && (tmp == "\n")) tmp = Jfile.get();
				}
				J[i][j] = atoi(val.c_str());
			}
			//while (tmp != "\n") {tmp = Jfile.get();}
		}
	}
	else
	{
		cout << "Gond van a " << fileName << " file-lal. Kerem kapcsolja ki!" << endl;
		system("pause");
	}

	return J;
}

//s beolvasasa
vector<bool> load_SpinBool_CG(int N, string fileName)
{
	vector<bool> s(N);
	string tmp = "";
	
	char *fileName_char = (char*)fileName.c_str();

	ifstream Sfile(fileName_char);
	if (Sfile.is_open())
	{
		while (tmp != "\n") {tmp = Sfile.get();}	//1. sor atugrasa

		for (int i=0; i<N; i++)
		{
			while (tmp != "\t") {tmp = Sfile.get();}
			tmp = Sfile.get();
			s[i] = (tmp == "1") ? true : false;
		}
	}
	else
	{
		cout << "Gond van a " << fileName << " file-lal. Kerem kapcsolja ki!" << endl;
		system("pause");
	}	
	
	return s;
}

//inhomogen h beolvasasa
vector<double> load_inhomogenH(int N, string fileName)
{
	vector<double> h(N);
	string tmp = ""; string val = "";
	
	char *fileName_char = (char*)fileName.c_str();

	ifstream Hfile(fileName_char);
	if (Hfile.is_open())
	{
		while (tmp != "\n") {tmp = Hfile.get();}	//1. sor atugrasa

		for (int i=0; i<N; i++)
		{
			tmp = "";
			val = "";
			while (tmp != "\t") {tmp = Hfile.get();}

			while (tmp != "\n")
			{
				val += tmp;
				tmp = Hfile.get();
			}

			h[i] = atof(val.c_str());
		}
	}
	else
	{
		cout << "Gond van a " << fileName << " file-lal. Kerem kapcsolja ki!" << endl;
		system("pause");
	}	
	
	return h;
}


//#endif