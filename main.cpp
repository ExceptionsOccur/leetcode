
#include <algorithm>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <numeric>

using namespace std;

// �����νṹ�������С·��������
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


// ��β��ӵ�С͵����
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

// n��1�Ӽ�����ΪĿ��ֵ�����������
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
class Solution {
public:
    bool canPartition(vector<int>& nums) {
        int sum = accumulate(nums.begin(), nums.end(), 0);
        if (sum % 2)   return false;
        int target = sum / 2;
        if (find(nums.begin(), nums.end(), target) != nums.end())  return true;
        vector<int> dp(target+1,0);
        dp[0] = 1;
        for (int i : nums) {
            for (int j = target; j > 0; j--) {
                if (j >= i)
                    dp[j] = dp[j] || dp[j - i];
            }
        }
        return dp[target];
    }
};

int main() {
    vector<int> data = { {1,2,5} };
    Solution solution;
    cout << solution.canPartition(data);
}