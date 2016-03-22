#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <queue>
#include <unistd.h>
#include <stdint.h>
using namespace std;


class ALIGNMENT
{
	public:
		ALIGNMENT()
		{
			cy = 3; cn = -1; cg = -2;
			INF = 0xFFFFFFFF;
		}
		~ALIGNMENT()
		{

		}

		void get_local(string s11, string s22, string &s3, string &s4)
		{
			s1 = s11;
			s2 = s22;
			// cout << s1 << endl << s2 << endl;
			alignment_state = LOCAL;
			set_dp();
			v.clear();
			size_t a = s1.size(), b = s2.size();
			double k = dp[a][b];
			for (size_t i = 0; i <= s1.size(); ++i)
			{
				for (size_t j = 0; j<=s2.size(); ++j)
				{
					if (dp[i][j] > k)
					{
						k = dp[i][j];
						a = i;
						b = j;
					}
				}
			}
			// cout << endl << k << endl;

			get_dp(a, b, v);

			s_1.clear();
			s_2.clear();

			for (size_t i=0, m=L[0], n=L[1]; i<v.size(); ++i) 
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
				// printf("%d ", v[i]);
			}
			s3 = s_1;
			s4 = s_2;
			// cout << endl << s_1 << endl << s_2 << endl;
		}

		void get_semi(string s11, string s22, string &s3, string &s4)
		{

			s1 = s11;
			s2 = s22;
			// cout << s1 << endl << s2 << endl;
			alignment_state = SEMI;
			set_dp();
			v.clear();
			size_t a = s1.size(), b = s2.size();
			double k = dp[a][b];
			for (size_t i = 0; i < s1.size(); ++i)
			{
				if (dp[i][s2.size()] > k) 
				{
					k = dp[i][s2.size()];
					a = i;
					b = s2.size();
				}
			}
			for (size_t i = 0; i < s2.size(); ++i)
			{
				if (dp[s1.size()][i] > k) 
				{
					k = dp[s1.size()][i];
					a = s1.size();
					b = i;
				}
			}
			// cout << endl << k << endl;
			get_dp(a, b, v);

			s_1.clear();
			s_2.clear();

			for (size_t i=0, m=L[0], n=L[1]; i<v.size(); ++i) 
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
				// printf("%d ", v[i]);
			}
			s3 = s_1;
			s4 = s_2;
			// cout << endl << s_1 << endl << s_2 << endl;
		}

		void get_global(string s11, string s22, string &s3, string &s4)
		{
			s1 = s11;
			s2 = s22;
			// cout << s1 << endl << s2 << endl;
			alignment_state = GLOBAL;
			set_dp();
			v.clear();
			// cout << endl << dp[s1.size()][s2.size()] << endl;
			get_dp(s1.size(), s2.size(), v);

			s_1.clear();
			s_2.clear();

			for (size_t i=L[0], m=L[1], n=0; i<v.size(); ++i) 
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
				// printf("%d ", v[i]);
			}
			s3 = s_1;
			s4 = s_2;
			// cout << endl << s_1 << endl << s_2 << endl;
		}
		
		void set_dp()
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
					if (alignment_state == GLOBAL)
					{
						dp[i][j] = -INF;
					}
					else if (alignment_state == SEMI)
					{
						if (i == 0) dp[i][j] = 0;
						else if (j == 0) dp[i][j] = 0;
						else dp[i][j] = -INF;
					}
					else
					{
						dp[i][j] = 0;
					}
					if (i==0 && j==0) dp[i][j] = 0;
					tmp = dp[i][j];
					if (i) tmp = std::max(tmp, dp[i-1][j] + cg);
					if (j) tmp = std::max(tmp, dp[i][j-1] + cg);
					if (i > 0 && j > 0)	tmp = std::max(tmp, dp[i-1][j-1] + ((s1[i-1]==s2[j-1]) ? cy : cn));
					dp[i][j] = tmp;
					// printf("%3.0lf ", dp[i][j]);
				}
				// printf("\n");
			}
		}

		void get_dp(size_t i, size_t j, vector<int> &v)
		{
			if (alignment_state == GLOBAL)
			{
				if (i==0 && j==0)
				{
					L[0] = i;
					L[1] = j;
					return;
				}
			}
			else if (alignment_state == SEMI)
			{
				if (i == 0  || j ==0)
				{
					L[0] = i;
					L[1] = j;
					return;
				}
			}
			else
			{
				if (dp[i][j] == 0) 
				{
					L[0] = i;
					L[1] = j;
					return;
				}
			}

			if (i && j && dp[i][j] == (dp[i-1][j-1] + ((s1[i-1]==s2[j-1]) ? cy : cn))) 
			{
				get_dp(i-1, j-1, v);
				v.push_back(0);
				return;
			}
			if (j && dp[i][j] == dp[i][j-1] + cg)
			{
				get_dp(i, j-1, v);
				v.push_back(2);
				return;
			}
			if (i && dp[i][j] == dp[i-1][j] + cg)
			{
				get_dp(i-1, j, v);
				v.push_back(1);
				return;
			}
		}

	private:
		double cy, cn, cg, INF;
		size_t L[2];
		vector< vector<double> > dp;
		enum{LOCAL, GLOBAL, SEMI} alignment_state;
		vector<int> v;
		string s1, s2, s_1, s_2;
};
