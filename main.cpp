
#include <algorithm>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <numeric>

using namespace std;

// 三角形结构数组的最小路径和问题
//class Solution {
//public:
//    int minimumTotal(vector<vector<int>>& triangle) {
//        if (triangle.size() == 0 || triangle[0].size() == 0)    return 0;
//        if (triangle.size() == 1) return triangle[0][0];
//        vector<int> start = triangle[0];
//        int min_num = 0;
//        vector<vector<int>> dp;
//        dp.push_back(start);
//        for (int i = 1; i < triangle.size(); i++) {
//            vector<int> temp(triangle[i].size(), 0);
//            dp.push_back(temp);
//            for (int j = 0; j < triangle[i].size(); j++) {
//                if (j == 0) {
//                    dp[i][j] = triangle[i][j] + dp[i-1][j];
//                    min_num = dp[i][j];
//                }
//                else if (j == triangle[i].size() - 1) {
//                    dp[i][j] = triangle[i][j] + dp[i - 1][j - 1];
//                    if (dp[i][j] < min_num)  min_num = dp[i][j];
//                }
//                else {
//                    dp[i][j] = triangle[i][j] + min(dp[i-1][j-1], dp[i-1][j]);
//                    if (dp[i][j] < min_num)  min_num = dp[i][j];
//                }
//            }
//        }
//        return min_num;
//    }
//
//};


// 首尾相接的小偷问题，第一个元素和最后一个元素取一个，因此分别去掉，求两种情况下的最大值即可
//class Solution {
//public:
//    int rob(vector<int>& nums) {
//        vector<int> nums1 = nums;
//        vector<int> nums2 = nums;
//        nums1.erase(nums1.begin());
//        nums2.erase(nums2.end() - 1);
//        int max_1 = 0, max_2 = 0;
//        int dp_1[100] = { 0 }, dp_2[100] = { 0 };
//
//        dp_1[0] = nums1[0];
//        dp_1[1] = max(nums1[0], nums1[1]);
//
//        dp_2[0] = nums2[0];
//        dp_2[1] = max(nums2[0], nums2[1]);
//
//        for (int i = 2; i < nums1.size(); i++) {
//            dp_1[i] = max(nums1[i] + dp_1[i - 2], dp_1[i - 1]);
//        }
//        max_1 = max(dp_1[nums1.size() - 1], dp_1[nums1.size() - 2]);
//
//        for (int i = 2; i < nums2.size(); i++) {
//            dp_2[i] = max(nums2[i] + dp_2[i - 2], dp_2[i - 1]);
//        }
//        max_2 = max(dp_2[nums2.size() - 1], dp_2[nums2.size() - 2]);
//
//        return max(max_1, max_2);
//    }
//};

// n个1加减运算为目标值的组合数问题，非leetcode题目，纯属是对原题目理解出现错误
//class Solution {
//public:
//    int fact(int n) {
//        if (n == 0 || 1) return 1;
//        return n * fact(n - 1);
//    }
//    int findTargetSumWays(vector<int>& nums, int S) {
//        int size = nums.size();
//        if ((size - S) % 2 == 1) return 0;
//        int pos_nums = (size + S) / 2;
//        int neg_nums = (size - S) / 2;
//        int temp = fact(pos_nums) * fact(neg_nums);
//        return fact(size) / temp;
//    }
//};

// 众数
//class Solution {
//public:
//    int majorityElement(vector<int>& nums) {
//        if (nums.empty())    return 0;
//        unordered_map<int, int> um;
//        int temp = 0, num = 0, max = 0;
//        for (int i : nums) {
//            if (um.find(i) == um.end()) {
//                um.insert({ i,1 });
//            }
//            else {
//                temp = um[i] + 1;
//                if (temp > max) { max = temp; num = i; }
//                um.erase(i);
//                um.insert({ i,temp });
//            }
//        }
//        return num;
//    }
//};

// 重量、价值相同的0-1背包问题
// lintcode 92
// 初始化dp[0] = 0, 保证了背包时容量为0时，已经放满，无法继续放入，即可放入价值为0
// dp数组为背包中价值
//class Solution {
//public:
//    int backPack(int target, vector<int> nums) {
//        vector<int> dp(target + 1,0);
//        if (nums.size() == 0) return 0;
//        for (int i : nums) {
//            for (int j = target; j > 0; j--) {
//                if (j >= i)
//                    dp[j] = max(dp[j], dp[j - i] + i);
//                else
//                    break;
//            }
//        }
//        return dp[target];
//    }
//};

// 部分和问题，转化为重量、价值相同的0-1背包问题
// lintcode 563
// 初始化dp[0] = 1, 保证了背包时容量为0时，已经放满，无法继续放入，即满足部分和条件，是问题的一个解
// dp数组为组合数
//class Solution {
//public:
//    int findTargetSumWays(vector<int>& nums, int S) {
//        int size = nums.size();
//        int sum = accumulate(nums.begin(), nums.end(), 0);
//        if ((sum - S) % 2 == 1 || sum < S) return 0;
//        int target = (sum + S) / 2;
//        vector<int> dp(target + 1, 0);
//        dp[0] = 1;
//        for (int num : nums) {
//            int i = target;
//            while (i >= num) {
//                dp[i] = dp[i] + dp[i - num];
//                i = i - 1;
//            }
//        }
//        return dp[target];
//    }
//};


// 有条件满足即可，可以是取当前元素满足，亦可以是上一次取值满足，因此取或操作
//class Solution {
//public:
//    bool canPartition(vector<int>& nums) {
//        int sum = accumulate(nums.begin(), nums.end(), 0);
//        if (sum % 2)   return false;
//        int target = sum / 2;
//        if (find(nums.begin(), nums.end(), target) != nums.end())  return true;
//        vector<int> dp(target+1,0);
//        dp[0] = 1;
//        for (int i : nums) {
//            for (int j = target; j > 0; j--) {
//                if (j >= i)
//                    dp[j] = dp[j] || dp[j - i];
//            }
//        }
//        return dp[target];
//    }
//};


// 整数拆分，求拆分后能得到的最大乘积
// 完全背包问题思路，难点在于dp的取值，作为分解数时无法取到本身，因此不能作为递推项
// n > 4 时，显然拆分后的最大积大于本身，n < 4 时本身大于拆分后的最大积，因此可以对小于4的情况特别处理
// 本代码使用特别处理方式，注释两行为统一处理
//class Solution {
//public:
//    int integerBreak(int n) {
//        if (n < 4 && n >= 2) return n - 1;
//        vector<int> dp(n + 1, 1);
//        dp[0] = 0;
//        dp[1] = 1;
//        dp[2] = 2;
//        dp[3] = 3;
//        // dp[1] = 0;
//        int k = 0;
//        for (int i = 4; i <= n; i++) {
//            for (int j = 1; j <= i; j++) {
//                k = i - j;
//                //dp[i] = max(dp[i], max(dp[j], j) * k);
//                dp[i] = max(dp[i], dp[j] * k);
//            }
//        }
//        return dp[n];
//    }
//};


// 与部分和问题区别，部分和问题结果无排列，本题结果有排列，循环部分不同
// 本题循环部分先取值再计数，不同取值顺序都可计算到，对比部分和问题的先取不同target值再取值，不区分排列
// 另一难点为判断target值是否能由所有值不同组合计算得到(太难了，就不优化了，暴力算，用 unsigned long long防止溢出）
class Solution {
public:
    unsigned long long combinationSum4(vector<int>& nums, int target) {
        vector<unsigned long long> dp(target + 1, 0);
        dp[0] = 1;
        for (int j = 1; j <= target; j++) {
            for (int i : nums) {
                if (j >= i)
                    dp[j] = dp[j] + dp[j - i];
            }
        }
        return dp[target];
    }
};

int main() {
    vector<int> data = { {1,2,5}};
    Solution solution;
    cout << solution.combinationSum4(data,6);
}