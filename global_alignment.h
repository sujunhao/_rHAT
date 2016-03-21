#include <iostream>
#include <string>
#include <fstream>
#include <bitset>
#include <cstring>
#include <cmath>
#include <algorithm> 
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <vector>


int32_t gdb[100][100];

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
			if (i > 0 && j > 0 && s1[i-1] == s2[j-1])	tmp = std::max(tmp, gdb[i-1][j-1] + 1);
			gdb[i][j] = tmp;
		}
	}

	size_t i = s1.size(), j = s2.size();
	std::string s[2];
	while (i + j)
	{
		if (i && j && gdb[i][j] == gdb[i-1][j-1] + 1)
		{
			s[1][i - 1] = s1[i - 1];
			s[2][j - 1] = s2[j - 1];
			--i;
			--j;
		}
		else if (i && gdb[i][j] == gdb[i-1][j])
		{
			s[1][i - 1] = s1[i - 1];
			s[2][j - 1] = ' ';
			--i;
		}
		else if (j && gdb[i][j] == gdb[i][j-1])
		{
			s[1][i - 1] = ' ';
			s[2][j - 1] = s2[j - 1];
			--j;
		}
	}

	std::cout << s[1] << std::endl << s[2] << std::endl;

}