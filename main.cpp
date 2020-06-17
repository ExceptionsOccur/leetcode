
#include <algorithm>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <numeric>

using namespace std;

// ��ά����·������
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


// ��ά����·������,���ϰ���
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

// ��ά����·������,��С·����
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

// �����νṹ�������С·��������
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


// ��β��ӵĴ�ҽ������⣬��һ��Ԫ�غ����һ��Ԫ��ȡһ������˷ֱ�ȥ��������������µ����ֵ����
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

// n��1�Ӽ�����ΪĿ��ֵ����������⣬��leetcode��Ŀ�������Ƕ�ԭ��Ŀ�����ִ���
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

// ����
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

// ��������ֵ��ͬ��0-1��������
// lintcode 92
// ��ʼ��dp[0] = 0, ��֤�˱���ʱ����Ϊ0ʱ���Ѿ��������޷��������룬���ɷ����ֵΪ0
// dp����Ϊ�����м�ֵ
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

// ���ֺ����⣬ת��Ϊ��������ֵ��ͬ��0-1��������
// lintcode 563
// ��ʼ��dp[0] = 1, ��֤�˱���ʱ����Ϊ0ʱ���Ѿ��������޷��������룬�����㲿�ֺ��������������һ����
// dp����Ϊ�����
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


// ���������㼴�ɣ�������ȡ��ǰԪ�����㣬���������һ��ȡֵ���㣬���ȡ�����
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


// ������֣����ֺ��ܵõ������˻�
// ��ȫ��������˼·���ѵ�����dp��ȡֵ����Ϊ�ֽ���ʱ�޷�ȡ��������˲�����Ϊ������
// n > 4 ʱ����Ȼ��ֺ���������ڱ���n < 4 ʱ������ڲ�ֺ����������˿��Զ�С��4������ر���
// ������ʹ���ر���ʽ��ע������Ϊͳһ����
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


// �벿�ֺ��������𣬲��ֺ������������У������������У�ѭ�����ֲ�ͬ
// ����ѭ��������ȡֵ�ټ�������ͬȡֵ˳�򶼿ɼ��㵽���ԱȲ��ֺ��������ȡ��ͬtargetֵ��ȡֵ������������
// ��һ�ѵ�Ϊ�ж�targetֵ�Ƿ���������ֵ��ͬ��ϼ���õ�(̫���ˣ��Ͳ��Ż��ˣ������㣬�� unsigned long long��ֹ�����
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


// �й��оصĶ�̬�滮������ʹ���Զ����£�������ʱ��
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


// ����ͬʱ�ӣ���dp[i][j]��ʾi�������ӳ���Ϊj�����������Ȼ��dp[i][j] = dp[i-1][j-1] + dp[i-1][j-2] +...+ dp[i-1][j-f]
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
//        // return dp[d][target] % 1000000007;   // ��������½������
//        // 1E9+7 �� 1E+9 ��Ϊ�������������ʱһ�㶼ʹ�� mod �Ľ������������ mod��������ײ����Ӳ����int����˲���� unsigned long long
//        return dp[d][target];
//    }
//};


// ��ȫ����������Σ���m = n + k*k����d[m] = d[n] + 1������״̬ת�Ʒ���d[m] = d[m - k * k] + 1;
// ��1��2��3���޷��ֽ�Ϊ����ƽ�������ʱΪȫ��ȡ1�����Գ�ʼ��ʱ��ʼ��Ϊ���ڵ��ڱ����ֵ
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


// �� leetcode 279 ���ƣ��ؼ����ж���������������Ԫ������Ƿ��ܴճ�Ŀ��ֵ�����ｫ��ʼֵ��ʼ��Ϊ��Ŀ��ֵ�����
// ��������϶��޷��ճ�Ŀ��ֵ��������䣬��ȻΪ��ʼ��ʱ��ֵ�����������1�����п���ʹ�����ս��ΪĿ��ֵ������˱��뽫�������
// ��ʼ��Ϊ��Ŀ��ֵ���ֵ�������ڽ������ʱ���ô���
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
// ״̬ת�Ʒ��̲����г�����[x,y]�ĵ�����ƶ���[x-1,y],[x,y-1],[x+1,y],[x,y+1]���ƶ���ʣ���ƶ�������һ
// ��dp[x,y,s]��ʾ[x,y]�ĵ��ƶ�s�����ܲ����ĳ���·��������dp[x,y,s] = [x-1,y,s-1] + [x,y-1,s-1] + [x+1,y,s-1] + [x,y+1,s-1]
// ����ʹ�ü�����(0�С�m+2��)������(0�С�m+2��)��Ϊ�����ж����ⲿ�ֳ�ʼ��Ϊ1������Ϊ0
// vsʹ�õı��������鲻�ܲ�����ʼ����̫��������

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


// leetcode 338 ���㣡
// �Զ���������������һλ��1��λ����ԭ����������ֻ��һλ�Ĳ��ԭ����������λ�Ƿ�Ϊ1
// ���Եõ�״̬ת������ dp[n] = dp[n >> 1] + (n & 1)

//class Solution {
//public:
//    vector<int> countBits(int num) {
//        vector<int> dp(num + 1, 0);
//        for (int i = 0; i <= num; i++)
//            dp[i] = dp[i >> 1] + (i & 1);
//        return dp;
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

// leetcode 16.11 ��ˮ��
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
        for (size_t i = shorter * k; i <= longer * k; i += (longer - shorter)) {
            result.push_back(i);
        }
        return result;
    }
};

int main() {
    vector<int> data = {};
    Solution solution;
    data = solution.divingBoard(1, 2, 3);
}