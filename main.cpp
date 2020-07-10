
#include <algorithm>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <numeric>
#include <random>

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


// 骰子同时扔，以dp[i][j]表示i个骰子扔出和为j的组合数，显然，dp[i][j] = dp[i-1][j-1] + dp[i-1][j-2] +...+ dp[i-1][j-f]
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
//        for (int i = 1; i <= target; i++) {
//            for (int j = 2; j <= d; j++) {
//                for (int k = 1; k <= f; k++) {
//                    if (i > k)
//                        dp[j][i] = (dp[j][i] + dp[j - 1][i - k]) % 1000000007;
//                }
//            }
//        }
//        // return dp[d][target] % 1000000007;   // 会溢出导致结果错误
//        // 1E9+7 与 1E+9 都为质数，结果过大时一般都使用 mod 的结果，对质数求 mod，减少碰撞且相加不溢出int，相乘不溢出 unsigned long long
//        return dp[d][target];
//    }
//};


// 完全背包问题变形，若m = n + k*k，则d[m] = d[n] + 1；即有状态转移方程d[m] = d[m - k * k] + 1;
// 由1，2，3等无法分解为例，平方数最多时为全部取1，所以初始化时初始化为大于等于本身的值
// leetcode 279
//class Solution {
//public:
//    int numSquares(int n) {
//        vector<int> dp(n + 1, n);
//        dp[0] = 0;
//        for (int i = 1; i <= n; i++) {
//            for (int j = 1; j * j <= i; j++) {
//                dp[i] = min(dp[i], dp[i - j * j] + 1);
//            }  
//        }
//        return dp[n];
//    }
//};


// 与 leetcode 279 类似，关键在判断所给的数组所有元素组合是否能凑出目标值，这里将初始值初始化为比目标值大的数
// 当所有组合都无法凑出目标值，结果不变，仍然为初始化时的值；当数组存在1，则有可能使得最终结果为目标值本身，因此必须将结果数组
// 初始化为比目标值大的值，这样在结果返回时更好处理
// leetcode 322
//class Solution {
//public:
//    int coinChange(vector<int>& coins, int amount) {
//        vector<int> dp(amount + 1, amount + 1);
//        dp[0] = 0;
//        for (int i = 1; i <= amount; i++) {
//            for (int j : coins) {
//                if (i >= j)
//                    dp[i] = min(dp[i], dp[i - j] + 1);
//            }
//        }
//        return (dp[amount] > amount) ? -1 : dp[amount];
//    }
//};


// leetcode 576 
// 状态转移方程不难列出来，[x,y]的点可以移动到[x-1,y],[x,y-1],[x+1,y],[x,y+1]，移动后剩余移动步数减一
// 设dp[x,y,s]表示[x,y]的点移动s步可能产生的出界路径数，则dp[x,y,s] = [x-1,y,s-1] + [x,y-1,s-1] + [x+1,y,s-1] + [x,y+1,s-1]
// 这里使用加两行(0行、m+2行)、两列(0列、m+2列)作为出界判定，这部分初始化为1，其他为0
// vs使用的编译器数组不能参数初始化，太不方便了

//class Solution {
//public:
//    int findPaths(int m, int n, int N, int i, int j) {
//        int x = i + 1, y = j + 1, M = m + 2, n_ = n + 2;
//        vector<unsigned long long> temp(N + 1, 0);
//        vector<vector<unsigned long long>> temp_(n_, temp);
//        vector<vector<vector<unsigned long long>>> dp(M, temp_);
//        for (int a = 0; a < M; a++)
//            for (int b = 0; b < n_; b++)
//                for (int c = 0; c < N + 1; c++)
//                    if (a == 0 || a == (m + 1) || b == 0 || (b == n + 1))
//                        dp[a][b][c] = 1;
//
//        for (int c = 1; c <= N; c++)
//            for (int a = 1; a <= m; a++)
//                for (int b = 1; b <= n; b++)
//                    dp[a][b][c] = (dp[a - 1][b][c - 1] + dp[a][b - 1][c - 1] + dp[a + 1][b][c - 1] + dp[a][b + 1][c - 1]) % 1000000007;
//
//        return dp[x][y][N];
//    }
//};


// leetcode 338 天秀！
// 以二进制来讲，右移一位后1的位数与原来的数至多只有一位的差别，原来的数最右位是否为1
// 可以得到状态转换方程 dp[n] = dp[n >> 1] + (n & 1)

//class Solution {
//public:
//    vector<int> countBits(int num) {
//        vector<int> dp(num + 1, 0);
//        for (int i = 0; i <= num; i++)
//            dp[i] = dp[i >> 1] + (i & 1);
//        return dp;
//    }
//};

// leetcode 16.11 k个数，只有shorter和longer两个值，求所有可能和，由小到大排列
//class Solution {
//public:
//    vector<int> result;
//    vector<int> divingBoard(int shorter, int longer, int k) {
//        vector<int> result;
//        if (k == 0)
//            return result;
//        if (shorter == longer) {
//            result.push_back(shorter * k);
//            return result;
//        }
//        for (size_t i = shorter * k; i <= longer * k; i += (longer - shorter)) {
//            result.push_back(i);
//        }
//        return result;
//    }
//};


// leetcode 172 阶乘末尾0个数，只需要求所有乘数中5的因子数量，如25为2个，125为3个
//class Solution {
//public:
//    int trailingZeroes(int n) {
//        if (n < 5)   return 0;
//        int sum = 0;
//        int k = n / 5, c = 0;
//        for (int i = 1; i <= k; i++) {
//            if (i % 5 == 0) {
//                c = i;
//                while (c % 5 == 0) {
//                    c = c / 5;
//                    sum += 1;
//                }
//                sum += 1;
//            }
//            else
//                sum += 1;
//        }
//        return sum;
//    }
//};
//class Solution {
//public:
//    int trailingZeroes(int n) {
//        if (n < 5)   return 0;
//        return n / 5 + trailingZeroes(n / 5);
//    }
//};

// leetcode 470 拒绝采样
//class Solution {
//public:
//    // leetcode 本身提供rand7()，这里自己实现
//    int rand7() {
//        return 1 + rand() % 7;
//    }
//
//    int rand10() {
//        int m = 0, n = 0, result = 0;
//        do {
//            m = rand7();
//            n = rand7();
//            // 此处若使用 m*n ，在[1,m*n]范围内取值概率不均匀，且无法取到质数，映射无效
//            // 使用(n - 1) * 7在[1,m*n]范围均匀随机取0、7、14、21、28、35、42
//            // 再加上[1，7]范围内均匀随机生成的 m，即可以实现在[1,49]范围取得均匀随机数
//            result = (n - 1) * 7 + m;
//        } while (result > 40);
//        return (result - 1) % 10 + 1;
//    }
//};

// leetcode 292, 两人拿石块游戏，拿到最后一块胜利，判断是否可以胜利
//class Solution {
//public:
//    bool canWinNim(int n) {
//        if (n != 0 && n % 4 == 0)
//            return false;
//        return true;
//    }
//};

// leetcode 1137, 三递归会超时，使用动态规划思想，使用循环
//class Solution {
//public:
//    int tribonacci(int n) {
//        // if(n == 0){
//        //     return 0;
//        // }
//        // if(n == 1 || n == 2){
//        //     return  1;
//        // }
//        // return tribonacci(n - 1) + tribonacci(n - 2) + tribonacci(n - 3);
//        vector<int> result(40, 0);
//        result[0] = 0;
//        result[1] = 1;
//        result[2] = 1;
//        if (n < 3)   return result[n];
//        for (int i = 3; i <= n; i++) {
//            result[i] = result[i - 1] + result[i - 2] + result[i - 3];
//        }
//        return result[n];
//    }
//};

// leetcode面试题10-I Fibonacci数列
//class Solution {
//public:
//    int fib(int n) {
//        vector<int> result(128, 0);
//        result[0] = 0;
//        result[1] = 1;
//        if (n < 2)   return result[n];
//        for (int i = 2; i <= n; i++) {
//            result[i] = (result[i - 1] + result[i - 2]) % 1000000007;
//        }
//        return result[n];
//    }
//};

// leetcode 面试题08.05
//class Solution {
//public:
//    int multiply(int A, int B) {
//        // 指数加速
//        // A*B=(A*2)*(B/2)=(A*2)*(floor(B/2) + 2)=(A*2)*floor(B/2)+A*(B%2)
//        if (B == 1)  return A;
//        if (B == 0)  return 0;
//        if (B & 1)
//            return multiply(A << 1, B >> 1) + A;
//        else
//            return multiply(A << 1, B >> 1);
//    }
//};

// leetcode 1313
//class Solution {
//public:
//    vector<int> decompressRLElist(vector<int>& nums) {
//        vector<int> r;
//        for (int i = 1; i <= nums.size() / 2; i++) {
//            vector<int> temp(nums[2 * (i - 1)], nums[2 * (i - 1) + 1]);
//            r.insert(r.end(), temp.begin(), temp.end());
//        }
//        return r;
//    }
//};

// leetcode 260 只出现一次的数字
// 本题主要思想为分组，将只出现一次的两个数字分别分到两个组中做异或操作，变成 leetcode 136 题
//class Solution {
//public:
//    vector<int> singleNumber(vector<int>& nums) {
//        vector<int> result(2, 0);
//        int s(0);
//        int div(0);
//        for (int num : nums) {
//            s ^= num;
//        }
//        div = s & (-s);
//        for (int num : nums) {
//            if (num & div) {
//                result[0] ^= num;
//            }
//            else
//                result[1] ^= num;
//        }
//        return result;
//    }
//};

// leetcode 112
//class Solution {
//public:
// 
//    struct TreeNode {
//    int val;
//    TreeNode *left;
//    TreeNode *right;
//    TreeNode(int x) : val(x), left(NULL), right(NULL) {}
//    };
//
//    bool hasPathSum(TreeNode* root, int sum) {
//        if (root == NULL)    return false;
//        if (root->left == NULL && root->right == NULL)    return sum == root->val;
//        return hasPathSum(root->left, sum - root->val) || hasPathSum(root->right, sum - root->val);
//    }
//};

// leetcode 17.04
//class Solution {
//public:
//    int missingNumber(vector<int>& nums) {
//        int sum = accumulate(nums.begin(), nums.end(), 0);
//        int max = *max_element(nums.begin(), nums.end());
//        if (nums.size() == max + 1 || max == 0)  return max + 1;
//        return (1 + max) * max / 2 - sum == 0 ? 0 : (1 + max) * max / 2 - sum;
//    }
//};

// leetcode 1502
//class Solution {
//public:
//    bool canMakeArithmeticProgression(vector<int>& arr) {
//        if (arr.size() == 2) return true;
//        sort(arr.begin(), arr.end());
//        int diff = arr[1] - arr[0];
//        for (int i = 2; i < arr.size(); i++) {
//            if (arr[i] - arr[i - 1] != diff)
//                return false;
//        }
//        return true;
//    }
//};

// leetcode 16.07
//class Solution {
//public:
//    int maximum(int a, int b) {
//        // gnu/clang 编译器
//        //int x = a ^ b;
//        //int bit = 0;
//        //__asm__(
//        //    "bsr %1, %0"
//        //    : "=r" (bit)
//        //    : "r" (x)
//        //);
//        //int num = 1 << bit;
//        //return num & (1 << 31) ? (a & (1 << 31) ? b : a) : (a & num ? a : b);
//
//        // vs
//        int x = a ^ b;
//        int bit = 0;
//        __asm {
//            mov eax, x
//            bsr ebx, eax
//            mov bit, ebx
//        }
//        int num = 1 << bit;
//        return num & (1 << 31) ? (a & (1 << 31) ? b : a) : (a & num ? a : b);
//    }
//};

// leetcode 342
//class Solution {
//public:
//    bool isPowerOfFour(int num) {
//        if ((num & 0x80000000) == 0x80000000)    return false;
//        if ((num & (num - 1)) == 0 && num) {
//            while (num > 2) {
//                num >>= 2;
//            }
//            return num == 2 ? false : true;
//        }
//        return false;
//    }
//};

// leetcode 633
class Solution {
public:
    bool judgeSquareSum(int c) {
        int half = int(sqrt(c));
        for (unsigned int i = 0, j = half; i <= j;) {
            unsigned int sum = i * i + j * j;
            if (sum == c) return true;
            if (sum > c) j--;
            if (sum < c) i++;
        }
        return false;
    }
};

int main() {
    //vector<int> data = {3,5};
    Solution solution;
    cout << solution.judgeSquareSum(2147482647);
}