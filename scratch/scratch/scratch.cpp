// scratch.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <string>
int factorial(int n)
{
	int prod = 1;
	for(int i = 1; i <= n; ++i)
	{
		prod = prod*i;
	}
	return prod;
}

//int solution(std::vector<int> &A) {
//	// write your code in C++14 (g++ 6.2.0)
//	std::unordered_map<int, int> potential_pairs;
//	for (int i = 0, max = A.size(); i < max; ++i)
//	{
//		if (potential_pairs.count(A[i]) == 0)
//		{
//			potential_pairs[A[i]] = 0;
//		}
//		potential_pairs[A[i]] = potential_pairs[A[i]] + 1;
//	}
//
//	long num_of_pairs = 0;
//	for (auto it : potential_pairs)
//	{
//		int n = it.second;
//		num_of_pairs += (1.0 / 2.0)*(n)*(n - 1);
//		if (num_of_pairs > 1e9)
//		{
//			return 1e9;
//		}
//	}
//
//	return num_of_pairs;
//
//}

//int solution(std::string &S) {
//	// write your code in C++14 (g++ 6.2.0)
//	int max_len = 0;
//	for(int i = 0, n = S.length(); i < n; ++i)
//	{
//		std::string prefix = S.substr(0, i);
//		std::string suffix = S.substr(n-i, n);
//		if(prefix == suffix && i > max_len)
//		{
//			max_len = i;
//		}
//	}
//
//	return max_len;
//}

void dfs(int root, std::vector<int>& T, std::vector<int>& TR, std::unordered_set<int>& visited, int level, std::vector<int>& ans)
{
	ans[level] = ans[level] + 1;
	level = level + 1;
	visited.insert(root);
	if(visited.count(T[root]) == 0)
	{
		root = T[root];
		dfs(root, T, TR, visited, level, ans);
		//visited.erase(root);
	}
	else if(visited.count(TR[root]) == 0)
	{
		root = TR[root];
		dfs(root, T, TR, visited, level, ans);
		//visited.erase(root);
	}
}


std::vector<int> solution(std::vector<int> &T) {
	// write your code in C++14 (g++ 6.2.0)
	std::vector<int> TR(T.size());
	int capital = 0;
	int n = T.size() - 1;
	for (int i = 0; i < T.size(); ++i)
	{
		TR.insert(TR.begin() + T[i], i);
		if (T[i] == i) 
		{
			capital = i;
		}
	}

	std::unordered_set<int> visited;
	int level = 0;
	std::vector<int> ans(T.size(),0);

	dfs(capital, T, TR, visited, level, ans);

	return ans;


}

int main()
{
	std::vector<int> vec;
	vec.push_back(0);
	vec.push_back(1);
	vec.push_back(4);
	vec.push_back(9);
	vec.push_back(0);
	vec.push_back(5);
	//std::string s = "abbabba";
	solution(vec);
	std::cout << "Press any key" << std::endl;
	getchar();
	return 0;
}

