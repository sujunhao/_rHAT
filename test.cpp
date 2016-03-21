#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <queue>
#include <unistd.h>
#include <stdint.h>
// #include "global_alignment.h"
using namespace std;

int32_t gdb[100][100];

string s1, s2;


void get_global_alignment(size_t i, size_t j, vector<int> &v)
{
	if (i && j && gdb[i][j] == (gdb[i-1][j-1] + (s1[i-1] == s2[j-1]))) 
	{
		get_global_alignment(i-1, j-1, v);
		v.push_back(0);
		return;
	}
	if (i && gdb[i][j] == gdb[i-1][j])
	{
		get_global_alignment(i-1, j, v);
		v.push_back(1);
		return;
	}
	if (j && gdb[i][j] == gdb[i][j-1])
	{
		get_global_alignment(i, j-1, v);
		v.push_back(2);
		return;
	}
	if (i) 
	{
		get_global_alignment(i-1, j, v);
		v.push_back(1);
		return;
	}
	if (j) 
	{
		get_global_alignment(i, j-1, v);
		v.push_back(2);
		return;
	}
	v.push_back(0);
	return;

}

void Global_Alignment(const std::string s1, const std::string s2)
{
	int32_t tmp = 0;

	for (size_t i = 0; i <= s1.size(); ++i)
	{
		for (size_t j = 0; j <= s2.size(); ++j)
		{
			tmp = 0;
			if (i) tmp = std::max(tmp, gdb[i-1][j]);
			if (j) tmp = std::max(tmp, gdb[i][j-1]);
			if (i > 0 && j > 0)	tmp = std::max(tmp, gdb[i-1][j-1] + (s1[i-1] == s2[j-1]));
			gdb[i][j] = tmp;
			printf("%d ", gdb[i][j]);
		}
		printf("\n");
	}
}

class global_alignment
{
	public:
		global_alignment(string s11, string s22, string &s3, string &s4)
		{
			s1 = s11;
			s2 = s22;
			s_1 = s3;
			s_2 = s4;
		}
		global_alignment()
		{
			
		}
		void set_string(string s11, string s22, string &s3, string &s4)
		{
			s1 = s11;
			s2 = s22;
			s_1 = s3;
			s_2 = s4;
		}

};


int main()
{
	s1 = "AATTCCCCTAAAA";
	s2 = "TATTGCAGCTTAATTTTTT";
	cout << s1 << s2 << endl;
	Global_Alignment(s1, s2);
	string s[2];

	vector<int> v;
	get_global_alignment(s1.size(), s2.size(), v);

	for (size_t i=0, m=0, n=0; i<v.size(); ++i) 
	{
		switch(v[i])
		{
			case 0:
				s[0].append(s1, m++, 1);
				s[1].append(s2, n++, 1);
				break;
			case 1:
				s[0].append(s1, m++, 1);

				s[1].append("_");
				break;
			case 2:
				s[1].append(s2, n++, 1);
				s[0].append("_");
				break;

		}
		printf("%d ", v[i]);
	}
	cout << endl << s[0] << endl << s[1] << endl;
	return 0;
}