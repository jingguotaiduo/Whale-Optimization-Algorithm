import random
import numpy as np
import math
from matplotlib import pyplot as plt

class City(object):  # 城市类
    def __init__(self, ID, CoordinateX, CoordinateY):
        self.ID = ID
        self.CoordinateX = CoordinateX
        self.CoordinateY = CoordinateY

    def shuchu(self):
        print(self.ID, self.CoordinateX, self.CoordinateY)

def getDistanceOfTwoCity(city1, city2):  # 计算两个城市距离的函数
    return math.sqrt((pow((city1.CoordinateX - city2.CoordinateX), 2) + pow((city1.CoordinateY - city2.CoordinateY), 2)) / 10.0)

def ReplaceContinueousSpace(str): # 替换字符串中的空格，并以逗号分隔符连接
    n = len(str)
    newstr = ''
    for i in range(n):
        if str[i] != ' ':
            newstr += str[i]
            if i+1 < n and str[i+1] == ' ':
                newstr += ','
    return newstr


Distances = []  # 城市距离矩阵

class CityGraph(object):  # 城市群类
    def __init__(self, filename):
        self.Citys = []  # 所有的城市
        self.n = 0  # 城市的数量
        self.readDataFile(filename)
        self.computeDistances()

    def readDataFile(self, filename):  # 读取城市数据文件函数
        f = open(filename)
        line = f.readline()
        while line:
            line_str = line.strip('\n')
            city_linestr = ReplaceContinueousSpace(line_str).split(',')
            city_i = City(int(city_linestr[0]), float(city_linestr[1]), float(city_linestr[2]))
            self.Citys.append(city_i)
            line = f.readline()
        f.close()
        self.n = len(self.Citys)

    def computeDistances(self):  # 得到城市距离矩阵的函数
        for i in range(0, self.n):
            dist_i = []
            for j in range(0, self.n):
                dist_ij = getDistanceOfTwoCity(self.Citys[i], self.Citys[j])
                dist_i.append(dist_ij)
            Distances.append(dist_i)

    def shuchu(self):  # 输出城市群信息的函数
        print('城市ID-----城市横坐标X-----城市纵坐标Y')
        for i in range(self.n):
            self.Citys[i].shuchu()
        print('----------各城市间的距离矩阵---------')
        for i in range(0,self.n):
            for j in range(0, self.n):
                print(round(Distances[i][j],3), end=',')
            print('\n')

# The Whale Optimization Algorithm By jing_zhong 2022.6.2
class WOA(object):
    def __init__(self, population_size=50, dimension=29, lb=1, ub=29, Max_iteration = 100):
        self.dimension = dimension  # 每个鲸鱼个体位置的维度
        self.lb = lb  # 鲸鱼个体位置的下边界
        self.ub = ub  # 鲸鱼个体位置的上边界
        self.PopSize = population_size  # 鲸鱼种群大小
        self.Max_iter = Max_iteration  # 最大迭代次数
        self.Population = []  # 鲸鱼种群
        self.BestGeti = []  # 最优鲸鱼个体
        self.BestFitness = -100  # 最优鲸鱼个体的适应度
        self.Convergence_curve = []  # 各代鲸鱼种群最优适应度集合,用于绘制收敛曲线

    def Init_Population(self):  # 种群初始化，评价所有鲸鱼个体的适应度，获得最优个体
        self.Population = []
        for i in range(0, self.PopSize):
            geti_i = []
            for j in range(0, self.dimension):
                xj = j + 1
                geti_i.append(xj)
            random.shuffle(geti_i)
            self.Population.append(geti_i)
        self.BestGeti = self.Population[0]
        self.BestFitness = self.Evaluate(self.Population[0])
        for i in range(0,self.PopSize):
            fitness_i = self.Evaluate(self.Population[i])
            print('种群初始化第 {} 个服务链座鱼鲸个体 {} 的适应度为 {}'.format(i + 1, self.Population[i], fitness_i))
            if fitness_i < self.BestFitness:
                self.BestGeti = self.Population[i]
                self.BestFitness = fitness_i
        print('\n初始化种群最优服务链座鱼鲸个体 {} 的适应度为 {}'.format(self.BestGeti, self.BestFitness))

    def Valid_Geti(self,geti):
        for j in range(0, self.dimension):
            if geti[j] > self.lb and geti[j] < self.ub:
                geti[j] = int(math.floor(geti[j]))
            else:
                geti[j] = random.randint(int(self.lb), int(self.ub))
        # 调整鲸鱼个体解的有效性
        n = len(geti)
        std_geti = []  # 标准的一个全排列
        std_flag = []  # 标记是否需要修改
        Count = []  # 统计出现的次数
        for i in range(0, n):
            std_geti.append(i + 1)
            Count.append(0)
            std_flag.append(True)
        for i in range(0, n):
            findNum = std_geti[i]
            count = 0
            for j in range(0, n):
                if (findNum == geti[j]):
                    count = count + 1
                    if count > 1:
                        std_flag[j] = False
            Count[i] = count
        for i in range(0, n):
            if Count[i] == 0:
                for j in range(0, n):
                    if std_flag[j] == False:
                        geti[j] = std_geti[i]
                        std_flag[j] = True
                        break

    def Evaluate(self, geti):  # 鲸鱼个体适应度评价函数
        sumDist = 0
        n = len(geti)
        for i in range(0, n-1):
            sumDist += Distances[geti[i]-1][geti[i+1]-1]
        sumDist += Distances[geti[n-1]-1][geti[0]-1]
        return int(sumDist)

    def Start(self):
        self.Init_Population() # 初始化鲸鱼种群
        t = 0
        while t < self.Max_iter:
            a = 2 - t * (2 / self.Max_iter)  # a decreases linearly fron 2 to 0
            a2 = -1 + t * ((-1) / self.Max_iter)  # a2 linearly dicreases from -1 to -2 to calculate t
            for i in range(0, self.PopSize):
                r1 = random.random()  # r1 is a random number in [0,1]
                r2 = random.random()  # r2 is a random number in [0,1]
                A = 2 * a * r1 - a
                C = 2 * r2
                b = 1
                l = (a2 - 1) * random.random() +1
                p = random.random()
                for j in range(0,self.dimension):
                    if p < 0.5:
                        if abs(A) >= 1:
                            rand_leader_index = math.floor(self.PopSize * random.random())
                            X_rand = self.Population[rand_leader_index]
                            D_X_rand = abs(C * X_rand[j] - self.Population[i][j])
                            self.Population[i][j] = X_rand[j] - A * D_X_rand
                        else:
                            D_Leader = abs(C * self.BestGeti[j] - self.Population[i][j])
                            self.Population[i][j] = self.BestGeti[j] - A * D_Leader
                    else:
                        distance2Leader = abs(self.BestGeti[j] - self.Population[i][j])
                        self.Population[i][j] = distance2Leader * np.exp(b * l) * np.cos(l * 2 * math.pi) + self.BestGeti[j]
            t = t + 1
            self.Convergence_curve.append(self.BestFitness)
            for i in range(0,self.PopSize):
                self.Valid_Geti(self.Population[i]) # Return back the search agents that go beyond the boundaries of the search space
                fitness_i = self.Evaluate(self.Population[i])
                if fitness_i < self.BestFitness:
                    self.BestGeti = self.Population[i]
                    self.BestFitness = fitness_i
                print('第 {} 代种群第 {} 个座鱼鲸个体 {} 的适应度为 {}'.format(t + 1, i+1, self.Population[i], fitness_i))
            print('第 {} 代种群 最优座鱼鲸个体 {} 的适应度为 {}'.format(t + 1, self.BestGeti, self.BestFitness))
        print('{} 次迭代后WOA种群最优座鱼鲸个体 {} 的适应度值为 {}'.format(self.Max_iter, self.BestGeti, self.BestFitness))

if __name__ == "__main__":
    filename = 'bayg29.tsp'
    citygraph = CityGraph(filename)
    citygraph.shuchu()
    woa = WOA(50, 29, 1, 29, 500)
    woa.Start()
    epochs = range(len(woa.Convergence_curve))
    plt.figure()
    plt.plot(epochs, woa.Convergence_curve, "b", label="woa")
    plt.title('WOA_ServiceChain')
    plt.xlabel("Epochs")
    plt.ylabel("Fitness")
    plt.legend()
    plt.show()
