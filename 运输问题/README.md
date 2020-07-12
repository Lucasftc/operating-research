## README

本方案实现了运输问题的解决，其中traffic是主函数，接受价格矩阵，各地的生产量向量和各地的销售量向量作为输入，输出为调度矩阵和最终的最优总代价。其余函数为辅助函数：

+ greedy（）：利用贪心算法找到一个可行初始解
+ dual-variable 利用势差法计算对偶变量的值
+ check_sigma:计算各个空格的检验数
+ find_close_path：如果出现了负的检验数，寻找闭回路
+ adjust-distribution：跟据找到的闭回路调整调度
+ calculate：计算最终的总代价