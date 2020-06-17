
#include <algorithm>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <numeric>

using namespace std;

// 二维网格路径问题
// leetcode 62
//class Solution {
//public:
//    int uniquePaths(int m, int n) {
//        if (m == 1 || n == 1) return 1;
//        int dp[m][n];
//        dp[0][0] = 1;
//        for (int i = 1; i < m; i++) {
//            dp[i][0] = 1;
//        }
//        for (int i = 1; i < n; i++) {
//            dp[0][i] = 1;
//        }
//        for (int i = 1; i < m; i++) {
//            for (int j = 1; j < n; j++) {
//                dp[i][j] = dp[i - 1][j] + dp[i][j - 1];
//            }
//        }
//        return dp[m - 1][n - 1];
//    }
//};


// 二维网格路径问题,带障碍物
// leetcode 63
//class Solution {
//public:
//    int uniquePathsWithObstacles(vector<vector<int>>& obstacleGrid) {
//        if (obstacleGrid[0][0] == 1) return 0;
//        if (obstacleGrid.size() == 0) return 0;
//        if (obstacleGrid.size() == 1) {
//            int sum = 0;
//            for (int i = 0; i < obstacleGrid[0].size(); i++) {
//                sum = sum + obstacleGrid[0][i];
//            }
//            if (sum > 0) return 0;
//            return 1;
//        }
//        if (obstacleGrid[0].size() == 1) {
//            int sum = 0;
//            for (int i = 0; i < obstacleGrid.size(); i++) {
//                sum = sum + obstacleGrid[i][0];
//            }
//            if (sum > 0) return 0;
//            return 1;
//        }
//        int m = obstacleGrid.size();
//        int n = obstacleGrid[0].size();
//        long dp[m][n];
//        dp[0][0] = 1;
//        for (int i = 1; i < m; i++) {
//            if (obstacleGrid[i][0] == 0)  dp[i][0] = 1;
//            else {
//                for (int k = i; k < m; k++) {
//                    dp[k][0] = 0;
//                }
//                break;
//            }
//        }
//        for (int i = 1; i < n; i++) {
//            if (obstacleGrid[0][i] == 0)  dp[0][i] = 1;
//            else {
//                for (int k = i; k < n; k++) {
//                    dp[0][k] = 0;
//                }
//                break;
//            }
//        }
//        for (int i = 1; i < m; i++) {
//            for (int j = 1; j < n; j++) {
//                if (obstacleGrid[i][j] == 1)
//                    dp[i][j] = 0;
//                else
//                    dp[i][j] = dp[i - 1][j] + dp[i][j - 1];
//            }
//        }
//        return dp[m - 1][n - 1];
//    }
//};

// 二维网格路径问题,最小路径和
// leetcode 64
//class Solution {
//public:
//    int minPathSum(vector<vector<int>>& grid) {
//        if (grid.size() == 0 || grid[0].size() == 0) return 0;
//        if (grid.size() == 1) return accumulate(grid[0].begin(), grid[0].end(), 0);
//        if (grid[0].size() == 1) {
//            int sum = 0;
//            for (int i = 0; i < grid.size(); i++) {
//                sum = sum + grid[i][0];
//            }
//            return sum;
//        }
//        int m = grid.size();
//        int n = grid[0].size();
//        int dp[m][n];
//        dp[0][0] = grid[0][0];
//        for (int i = 1; i < m; i++) {
//            dp[i][0] = grid[i][0] + dp[i - 1][0];
//        }
//        for (int i = 1; i < n; i++) {
//            dp[0][i] = grid[0][i] + dp[0][i - 1];
//        }
//        for (int i = 1; i < m; i++) {
//            for (int j = 1; j < n; j++) {
//                dp[i][j] = grid[i][j] + min(dp[i - 1][j], dp[i][j - 1]);
//            }
//        }
//        return dp[m - 1][n - 1];
//    }
//};

// 三角形结构数组的最小路径和问题
// leetcode 120
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


// 首尾相接的打家劫舍问题，第一个元素和最后一个元素取一个，因此分别去掉，求两种情况下的最大值即可
// leetcode 213
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
// leetcode 416
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
// leetcode 343
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
// leetcode 377
//class Solution {
//public:
//    unsigned long long combinationSum4(vector<int>& nums, int target) {
//        vector<unsigned long long> dp(target + 1, 0);
//        dp[0] = 1;
//        for (int j = 1; j <= target; j++) {
//            for (int i : nums) {
//                if (j >= i)
//                    dp[j] = dp[j] + dp[j - i];
//            }
//        }
//        return dp[target];
//    }
//};


// 中规中矩的动态规划，这里使用自顶向下，非最优时间
// leetcode 931
//class Solution {
//public:
//
//    int min_in_three(int a, int b, int c) {
//        int min = 0;
//        min = a < b ? a : b;
//        min = min < c ? min : c;
//        return min;
//    }
//    int minFallingPathSum(vector<vector<int>>& A) {
//        if (A.size() == 0 || A[0].size() == 0)   return 0;
//        int size = A.size();
//        int min_ = INT_MAX;
//        vector<vector<int>> dp(size, A[0]);
//        for (int i = 1; i < size; i++) {
//            for (int j = 0; j < size; j++) {
//                if (j == 0)
//                    dp[i][j] = min(A[i][j] + dp[i - 1][j], A[i][j] + dp[i - 1][j + 1]);
//                else if (j == size - 1)
//                    dp[i][j] = min(A[i][j] + dp[i - 1][j], A[i][j] + dp[i - 1][j - 1]);
//                else
//                    dp[i][j] = min_in_three(A[i][j] + dp[i - 1][j - 1], A[i][j] + dp[i - 1][j], A[i][j] + dp[i - 1][j + 1]);
//            }
//        }
//        for (int r : dp[size - 1]) {
//            min_ = min(min_, r);
//        }
//        return min_;
//    }
//};


// leetcode 1155
//class Solution {
//public:
//    int numRollsToTarget(int d, int f, int target) {
//        vector<unsigned long long> temp(target + 1, 0);
//        vector<vector<unsigned long long>> dp(d + 1, temp);
//        for (int i = 1; i <= min(f, target); i++)
//            dp[1][i] = 1;
//        if (d == 1 && target <= f)   return  1;
//        if (d == 1 && target > f)  return 0;
//        for (int i = d; i <= target; i++) {
//            for (int j = 2; j <= d; j++) {
//                for (int k = 1; k <= f; k++) {
//                    if (i > k)
//                        dp[j][i] = dp[j][i] + dp[j - 1][i - k];
//                }
//            }
//        }
//        return dp[d][target] % 1000000007;
//    }
//};

// leetcode 16.11 跳水板
class Solution {
public:
    vector<int> divingBoard(int shorter, int longer, int k) {
        vector<int> result;
        if (k == 0)
            return result;
        if (shorter == longer) {
            result.push_back(shorter);
            return result;
        }
        for (size_t i = shorter * k; i <= longer * k; i += (longer - shorter)){
            result.push_back(i);
        }
        return result;
    }
};

int main() {
    vector<int> data = {};
    Solution solution;
    data = solution.divingBoard(1,2,3);
}