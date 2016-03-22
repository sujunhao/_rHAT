#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <queue>
#include <map>
#include <unistd.h>
#include <stdint.h>
// #include "global_alignment.h"
using namespace std;

const double cy = 3, cn = -1, cg = -2;
const double INF = 0xFFFFFFFF;

class alignment
{
	public:
		alignment()
		{
		}
		~alignment()
		{
		}


		void alignment_semi(string s11, string s22, string &s3, string &s4)
		{

			s1 = s11;
			s2 = s22;
			// cout << s1 << endl << s2 << endl;
			set_semi_dp();
			v.clear();
			size_t a = s1.size(), b = s2.size();
			double k = dp[a-1][b-1];
			for (size_t i = 0; i < s1.size(); ++i)
			{
				if (dp[i][s2.size() - 1] > k) 
				{
					k = dp[i][s2.size() - 1];
					a = i + 1;
					b = s2.size();
				}
			}
			for (size_t i = 0; i < s2.size(); ++i)
			{
				if (dp[s1.size() - 1][i] > k) 
				{
					k = dp[s1.size() - 1][i];
					a = s1.size();
					b = i + 1;
				}
			}
			get_semi_dp(a, b, v);

			s_1.clear();
			s_2.clear();

			for (size_t i=0, m=0, n=0; i<v.size(); ++i) 
			{
				switch(v[i])
				{
					case 0:
						s_1.append(s1, m++, 1);
						s_2.append(s2, n++, 1);
						break;
					case 1:
						s_1.append(s1, m++, 1);
						s_2.append("_");
						break;
					case 2:
						s_1.append("_");
						s_2.append(s2, n++, 1);
						break;

				}
				printf("%d ", v[i]);
			}
			s3 = s_1;
			s4 = s_2;
			cout << endl << s_1 << endl << s_2 << endl;
		}

		void alignment_g(string s11, string s22, string &s3, string &s4)
		{
			s1 = s11;
			s2 = s22;
			// cout << s1 << endl << s2 << endl;
			set_g_dp();
			v.clear();
			get_g_dp(s1.size(), s2.size(), v);

			s_1.clear();
			s_2.clear();

			for (size_t i=0, m=0, n=0; i<v.size(); ++i) 
			{
				switch(v[i])
				{
					case 0:
						s_1.append(s1, m++, 1);
						s_2.append(s2, n++, 1);
						break;
					case 1:
						s_1.append(s1, m++, 1);
						s_2.append("_");
						break;
					case 2:
						s_1.append("_");
						s_2.append(s2, n++, 1);
						break;

				}
				printf("%d ", v[i]);
			}
			s3 = s_1;
			s4 = s_2;
			cout << endl << s_1 << endl << s_2 << endl;
		}
		
		void set_semi_dp()
		{
			double tmp = 0.0;
			for (size_t i = 0; i <= s1.size(); ++i)
			{
				while (dp.size() <= i)
				{
					vector<double> v;
					v.clear();
					dp.push_back(v);
				}
				for (size_t j = 0; j <= s2.size(); ++j)
				{
					while (dp[i].size() <= j) dp[i].push_back(-INF);
					if (i == 0 || j == 0)
					{
						dp[i][j] = 0;
						continue;
					}
					if (i + j == 0)	dp[0][0] = 0.0;
					tmp = dp[i][j];
					if (i) tmp = std::max(tmp, dp[i-1][j] + cg);
					if (j) tmp = std::max(tmp, dp[i][j-1] + cg);
					if (i > 0 && j > 0)	tmp = std::max(tmp, dp[i-1][j-1] + ((s1[i-1]==s2[j-1]) ? cy : cn));
					dp[i][j] = tmp;
					printf("%3.0lf ", dp[i][j]);
				}
				printf("\n");
			}
		}

		void set_g_dp()
		{
			double tmp = 0.0;
			for (size_t i = 0; i <= s1.size(); ++i)
			{
				while (dp.size() <= i)
				{
					vector<double> v;
					v.clear();
					dp.push_back(v);
				}
				for (size_t j = 0; j <= s2.size(); ++j)
				{
					while (dp[i].size() <= j) dp[i].push_back(-INF);
					if (i + j == 0)	dp[0][0] = 0.0;
					tmp = dp[i][j];
					if (i) tmp = std::max(tmp, dp[i-1][j] + cg);
					if (j) tmp = std::max(tmp, dp[i][j-1] + cg);
					if (i > 0 && j > 0)	tmp = std::max(tmp, dp[i-1][j-1] + ((s1[i-1]==s2[j-1]) ? cy : cn));
					dp[i][j] = tmp;
					printf("%3.0lf ", dp[i][j]);
				}
				printf("\n");
			}
		}

		void get_semi_dp(size_t i, size_t j, vector<int> &v)
		{
			if (i && j && dp[i][j] == (dp[i-1][j-1] + ((s1[i-1]==s2[j-1]) ? cy : cn))) 
			{
				get_semi_dp(i-1, j-1, v);
				v.push_back(0);
				return;
			}
			if (j && dp[i][j] == dp[i][j-1] + cg)
			{
				get_semi_dp(i, j-1, v);
				v.push_back(2);
				return;
			}
			if (i && dp[i][j] == dp[i-1][j] + cg)
			{
				get_semi_dp(i-1, j, v);
				v.push_back(1);
				return;
			}
			// if (i) 
			// {
			// 	get_dp(i-1, j, v);
			// 	v.push_back(1);
			// 	return;
			// }
			// if (j) 
			// {
			// 	get_dp(i, j-1, v);
			// 	v.push_back(2);
			// 	return;
			// }
			// v.push_back(0);
			return;
		}

		void get_g_dp(size_t i, size_t j, vector<int> &v)
		{
			if (i && j && dp[i][j] == (dp[i-1][j-1] + ((s1[i-1]==s2[j-1]) ? cy : cn))) 
			{
				get_g_dp(i-1, j-1, v);
				v.push_back(0);
				return;
			}
			if (j && dp[i][j] == dp[i][j-1] + cg)
			{
				get_g_dp(i, j-1, v);
				v.push_back(2);
				return;
			}
			if (i && dp[i][j] == dp[i-1][j] + cg)
			{
				get_g_dp(i-1, j, v);
				v.push_back(1);
				return;
			}
			if (i) 
			{
				get_g_dp(i-1, j, v);
				v.push_back(1);
				return;
			}
			if (j) 
			{
				get_g_dp(i, j-1, v);
				v.push_back(2);
				return;
			}
			// v.push_back(0);
			return;
		}
	private:
		vector< vector<double> > dp;
		vector<int> v;
		string s1, s2, s_1, s_2;
};






int main()
{

	string s1 = "AAAAAAAAAAAAAAAATCGGG", s2 = "GCATCGGGGGGG";
	
	string s3, s4;
	alignment ga;
	ga.alignment_g(s1, s2, s3, s4);
	ga.alignment_semi(s1, s2, s3, s4);
	// cout << s1 << endl << s2 << endl;
	// cout << endl <<s3 << endl << s4 << endl;
	// Global_Alignment(s1, s2);
	// string s[2];

	// vector<int> v;
	// get_global_alignment(s1.size(), s2.size(), v);

	// for (size_t i=0, m=0, n=0; i<v.size(); ++i) 
	// {
	// 	switch(v[i])
	// 	{
	// 		case 0:
	// 			s[0].append(s1, m++, 1);
	// 			s[1].append(s2, n++, 1);
	// 			break;
	// 		case 1:
	// 			s[0].append(s1, m++, 1);

	// 			s[1].append("_");
	// 			break;
	// 		case 2:
	// 			s[1].append(s2, n++, 1);
	// 			s[0].append("_");
	// 			break;

	// 	}
	// 	printf("%d ", v[i]);
	// }
	// cout << endl << s[0] << endl << s[1] << endl;
	return 0;
}