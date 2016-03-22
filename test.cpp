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


class global_alignment
{
	public:
		global_alignment()
		{
		}
		~global_alignment()
		{
			
		}
		void set(string s11, string s22, string &s3, string &s4)
		{
			s1 = s11;
			s2 = s22;
			cout << s1 << endl << s2 << endl;
			set_dp();
			v.clear();
			get_dp(s1.size(), s2.size(), v);

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
				// printf("%d ", v[i]);
			}
			s3 = s_1;
			s4 = s_2;
			// cout << endl << s_1 << endl << s_2 << endl;
		}

		void set(string s11, string s22)
		{
			s1 = s11;
			s2 = s22;
			cout << s1 << endl << s2 << endl;
			set_dp();
			v.clear();
			get_dp(s1.size(), s2.size(), v);

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
				// printf("%d ", v[i]);
			}
			cout << endl << s_1 << endl << s_2 << endl;
		}


		void set_dp()
		{
			int32_t tmp = 0;

			for (size_t i = 0; i <= s1.size(); ++i)
			{
				while (dp.size() <= i)
				{
					vector<int32_t> v;
					v.clear();
					dp.push_back(v);
				}
				for (size_t j = 0; j <= s2.size(); ++j)
				{
					tmp = 0;
					while (dp[i].size() <= j) dp[i].push_back(0);
					if (i) tmp = std::max(tmp, dp[i-1][j]);
					if (j) tmp = std::max(tmp, dp[i][j-1]);
					if (i > 0 && j > 0)	tmp = std::max(tmp, dp[i-1][j-1] + (s1[i-1] == s2[j-1]));
					dp[i][j] = tmp;
					// printf("%d ", dp[i][j]);
				}
				// printf("\n");
			}
		}


		void get_dp(size_t i, size_t j, vector<int> &v)
		{
			if (i && j && dp[i][j] == (dp[i-1][j-1] + (s1[i-1] == s2[j-1]))) 
			{
				get_dp(i-1, j-1, v);
				v.push_back(0);
				return;
			}
			if (i && dp[i][j] == dp[i-1][j])
			{
				get_dp(i-1, j, v);
				v.push_back(1);
				return;
			}
			if (j && dp[i][j] == dp[i][j-1])
			{
				get_dp(i, j-1, v);
				v.push_back(2);
				return;
			}
			if (i) 
			{
				get_dp(i-1, j, v);
				v.push_back(1);
				return;
			}
			if (j) 
			{
				get_dp(i, j-1, v);
				v.push_back(2);
				return;
			}
			v.push_back(0);
			return;
		}



	private:
		vector< vector<int32_t> > dp;
		vector<int> v;
		string s1, s2, s_1, s_2;
};



int main()
{
	string s1 = "AAAAGGGTTGTGGGGGTACCCGTTTGGTTGTG";
	string s2 = "ATAAGGTTGTTACCCAAAAGGGTGTCG";
	string s3, s4;
	global_alignment ga;
	ga.set(s1, s2, s3, s4);
	cout << endl <<s3 << endl << s4 << endl;
	// cout << s1 << s2 << endl;
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