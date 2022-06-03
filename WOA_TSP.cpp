#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <ctime>

using namespace std;

const int citycount = 29;//城市的数量
#define PI       3.14159265358979323846   // pi
#define N  999 //精度为小数点后面3位
double round(double r)
{
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

class City//城市类
{
    public:
        string ID;//城市名称ID
        double x, y;//城市点的二维坐标x和y
        void shuchu()
        {
            std::cout <<ID+":"<<"("<< x << "," << y <<")"<< endl;
        }
};

class CityGraph//城市图
{
    public:
        City city[citycount];//城市数组
        double distance[citycount][citycount];//城市间的距离矩阵
        void Readcoordinatetxt(string txtfilename)//读取城市坐标文件的函数
        {
            ifstream myfile(txtfilename, ios::in);
            double x = 0, y = 0;
            int z = 0;
            if (!myfile.fail())
            {
                int i = 0;
                while (!myfile.eof() && (myfile >> z >> x >> y))
                {
                    city[i].ID = to_string(z);//城市名称转化为字符串
                    city[i].x = x; city[i].y = y;
                    i++;
                }
            }
            else
                cout << "文件不存在";
            myfile.close();//计算城市距离矩阵
            for (int i = 0; i < citycount; i++)
                for (int j = 0; j < citycount; j++)
                {
                    distance[i][j] = sqrt((pow((city[i].x - city[j].x), 2) + pow((city[i].y - city[j].y), 2))/10.0);//计算城市ij之间的伪欧式距离
                    if (round(distance[i][j] < distance[i][j]))distance[i][j] = round(distance[i][j]) + 1;
                    else distance[i][j] = round(distance[i][j]);
                }
        }
        void shuchu()
        {
            cout << "城市名称 " << "坐标x" << " " << "坐标y" << endl;
            for (int i = 0; i < citycount; i++)
                city[i].shuchu();
            cout << "距离矩阵： " << endl;
            for (int i = 0; i < citycount; i++)
            {
                for (int j = 0; j < citycount; j++)
                {
                    if (j == citycount - 1)
                       cout << distance[i][j] << endl;
                    else
                        cout << distance[i][j] << "  ";
                }
            }
        }
};


CityGraph Map_City;//全局对象城市图

int * Random_N(int n)// 随机生成一个个体解的函数
{
	int *geti;
	geti = new int[n];
	int j = 0;
	while(j<n)
	{
		while (true)
		{
			int flag = -1;
			int temp = rand() % n + 1;
			if (j > 0)
			{
				int k = 0;
				for(; k < j; k++)
				{
					if (temp == *(geti + k))break;
				}
				if (k == j)
				{
					*(geti + j) = temp;
					flag = 1;
				}
			}
			else
			{
				*(geti + j) = temp;
				flag = 1;
			}
			if (flag == 1)break;
		}
		j++;
	}
	return geti;
}

class JingYu//鲸鱼类，一个座头鲸代表一个个体解
{
	public:
		int *X;//鲸鱼位置向量
		double fitness;//鲸鱼的适应度
		void Init()
		{
			X = new int[citycount];
			int *M = Random_N(citycount);
			for (int j = 0; j < citycount; j++)
				X[j] = *(M + j);
			fitness=0;
		}
		void shuchu()
		{
			for (int j = 0; j < citycount; j++)
			{
				if (j == citycount -1) std::cout << X[j] << " "<<fitness<<endl;
				else std::cout << X[j] << "->";
			}
		}
};
class WOA
{
    public:
        JingYu *Population;//下一代种群
        int dimension;
        double LB;//位置向量每个分量的下确界
        double UB;//位置向量每个分量的上确界
        int Pop_Size;//种群大小
        int Itetime;//迭代次数
        JingYu Bestgeti;//最优鲸鱼个体
        double BestFitness;//最优鲸鱼个体的适应度


        void Init(int popsize,double dimen,double lb,double ub,int itetime,string filename)
        {
            Map_City.Readcoordinatetxt(filename);
            Map_City.shuchu();
            Pop_Size = popsize;
            dimension = dimen;
            LB = lb;
            UB = ub;
            Itetime = itetime;
            Population = new JingYu[Pop_Size];
            for (int i = 0; i < Pop_Size; i++)
            {
                Population[i].Init();
                Population[i].fitness = Evaluate(Population[i]);
            }
            Bestgeti.Init();
            for (int j = 0; j < citycount; j++)
                Bestgeti.X[j] = Population[0].X[j];
            Bestgeti.fitness = Evaluate(Bestgeti);
            for (int i = 0; i < Pop_Size; i++)
            {
                if (Population[i].fitness < Bestgeti.fitness)
                {
                    for (int j = 0; j < citycount; j++)
                        Bestgeti.X[j] = Population[i].X[j];
                    Bestgeti.fitness = Evaluate(Bestgeti);
                }
            }
            cout<<"初始化鲸鱼种群如下："<<endl;
            for (int i = 0; i < Pop_Size; i++)
                Population[i].shuchu();
            cout<<"初始化最优的鲸鱼个体如下："<<endl;
            Bestgeti.shuchu();
        }

        void Adjuxt_validGeti(JingYu jy)//调整鲸鱼个体位置有效性的函数
        {
            for(int j=0;j<dimension;j++)
            {

                if(jy.X[j] > LB && jy.X[j] < UB)
                    jy.X[j] = (int)jy.X[j];
                else
                    jy.X[j] = int(rand() % dimension);
            }
            int route[citycount];
            bool flag[citycount];
            int biaoji[citycount];
            for (int j = 0; j < citycount; j++)
            {
                route[j] = j + 1;
                flag[j] = false;
                biaoji[j] = 0;
            }

            for (int j = 0; j < citycount; j++)
            {
                int num = 0;
                for (int k = 0; k < citycount; k++)
                {
                    if (jy.X[k] == route[j])
                    {
                        biaoji[k] = 1;
                        num++; break;
                    }
                }
                if (num == 0) flag[j] = false;
                else if (num == 1) flag[j] = true;
            }
            for (int k = 0; k < citycount; k++)
            {
                if (flag[k] == false)
                {
                    int i = 0;
                    for (; i < citycount; i++)
                    {
                        if (biaoji[i] != 1)break;
                    }
                    jy.X[i] = route[k];
                    biaoji[i] = 1;
                }
            }
        }
        double Evaluate(JingYu jy)//计算鲸鱼个体适应值的函数
        {
            double Sum_dist=0;
            for (int j = 0; j < citycount-1; j++)
                Sum_dist += Map_City.distance[jy.X[j]][jy.X[j + 1]];
            Sum_dist += Map_City.distance[jy.X[citycount -1]][jy.X[0]];
            jy.fitness = Sum_dist;
            return Sum_dist;
        }

        void UpdateBestgeti()//更新鲸鱼最优个体的函数
        {
            for (int i = 0; i < Pop_Size; i++)
            {
                if (Population[i].fitness < Bestgeti.fitness)
                {
                    for (int j = 0; j < citycount; j++)
                        Bestgeti.X[j] = Population[i].X[j];
                }
            }
            Bestgeti.fitness = Evaluate(Bestgeti);
        }

        void ShuchuPopulation(int ite)
        {
            for (int i = 0; i < Pop_Size; i++)
            {
                cout << "第"<<ite<<"代种群中第"<<i + 1<<"个鲸鱼"<<"->";
                for (int j = 0; j < dimension; j++)
                {
                    if (j == dimension - 1) cout << Population[i].X[j] <<")对应的适应度为： "<<Evaluate( Population[i])<< endl;
                    else if(j==0)
                        cout <<"("<<  Population[i].X[j] << ", ";
                    else
                        cout <<  Population[i].X[j] << ", ";
                }
            }
        }
        void ShuchuBestgeti(int ite)
        {
            for (int j = 0; j < dimension; j++)
            {
                if (j == dimension - 1) std::cout << Bestgeti.X[j] << ")" << "对应的适应度为：" << Evaluate(Bestgeti) << endl;
                else if (j == 0) std::cout <<  "第"<<ite<<"代种群最优鲸鱼个体(" << Bestgeti.X[j] << ",";
                else std::cout << Bestgeti.X[j] << ",";
            }
        }
        void WOA_TSP(int popsize,double dimen,double lb,double ub,int Max_iter,string filename)
        {
            ofstream outfile;
            outfile.open("result.txt",ios::trunc);
            Init(popsize,dimen,lb,ub,Max_iter,filename);

            outfile<<"初始化鲸鱼种群如下："<<endl;
            for (int i = 0; i < Pop_Size; i++)
            {
                outfile << "初始化种群中第"<<i + 1<<"个鲸鱼"<<"->";
                for (int j = 0; j < dimension; j++)
                {
                    if (j == dimension - 1) outfile << Population[i].X[j] <<")对应的适应度为： "<<Evaluate( Population[i])<< endl;
                    else if(j==0)
                        outfile <<"("<<  Population[i].X[j] << ", ";
                    else
                        outfile <<  Population[i].X[j] << ", ";
                }
            }
            outfile<<"初始化最优的鲸鱼个体如下："<<endl;
            for (int j = 0; j < dimension; j++)
            {
                if (j == dimension - 1) outfile << Bestgeti.X[j] << ")" << "对应的适应度为：" << Evaluate(Bestgeti) << endl;
                else if (j == 0) outfile <<  "初始种群最优鲸鱼个体(" << Bestgeti.X[j] << ",";
                else outfile << Bestgeti.X[j] << ",";
            }

            int t=0;
            while (t < Itetime)
            {
                double a = 2 - t * (2 / Itetime);  // a decreases linearly fron 2 to 0
                double a2 = -1 + t * ((-1) / Itetime);  // a2 linearly dicreases from -1 to -2 to calculate t
                for (int i = 0; i < Pop_Size; i++)
                {
                    double r1 = rand() % (N + 1) / (float)(N + 1);  // r1 is a random number in [0,1]
                    double r2 = rand() % (N + 1) / (float)(N + 1);  // r2 is a random number in [0,1]
                    double A = 2 * a * r1 - a;
                    double C = 2 * r2;
                    double b = 1;
                    double l = (a2 - 1) * (rand() % (N + 1) / (float)(N + 1)) +1;
                    double p = rand() % (N + 1) / (float)(N + 1);
                    for(int j=0;j<dimension;j++)
                    {
                        if(p < 0.5)
                        {
                            if(abs(A) >= 1)
                            {
                                int rand_leader_index = rand() % dimension;
                                JingYu X_rand;
                                X_rand.Init();
                                for(int k=0;k<dimension;k++)
                                    X_rand.X[k] = Population[rand_leader_index].X[k];
                                double D_X_rand = abs(C * X_rand.X[j] - Population[i].X[j]);
                                Population[i].X[j] = int(X_rand.X[j] - A * D_X_rand);
                            }
                            else
                            {
                                double D_Leader = abs(C * Bestgeti.X[j] - Population[i].X[j]);
                                Population[i].X[j] = Bestgeti.X[j] - A * D_Leader;
                            }
                        }
                        else
                        {
                            double distance2Leader = abs(Bestgeti.X[j] - Population[i].X[j]);
                            Population[i].X[j] = int(distance2Leader * exp(b * l) * cos(l * 2 * PI) + Bestgeti.X[j]);
                        }
                    }
                }
                t = t + 1;
                for (int i = 0; i < Pop_Size; i++)
                {
                    Adjuxt_validGeti(Population[i]);
                    Population[i].fitness = Evaluate(Population[i]);
                }
                UpdateBestgeti();
                ShuchuPopulation(t);
                ShuchuBestgeti(t);
                //写到文本文件中
                for (int i = 0; i < Pop_Size; i++)
                {
                    outfile << "第"<<t<<"代种群中第"<<i + 1<<"个鲸鱼"<<"->";
                    for (int j = 0; j < dimension; j++)
                    {
                        if (j == dimension - 1) outfile << Population[i].X[j] <<")对应的适应度为： "<<Evaluate( Population[i])<< endl;
                        else if(j==0)
                            outfile <<"("<<  Population[i].X[j] << ", ";
                        else
                            outfile <<  Population[i].X[j] << ", ";
                    }
                }
                for (int j = 0; j < dimension; j++)
                {
                    if (j == dimension - 1) outfile << Bestgeti.X[j] << ")" << "对应的适应度为：" << Evaluate(Bestgeti) << endl;
                    else if (j == 0) outfile <<  "第"<<t<<"代种群最优鲸鱼个体(" << Bestgeti.X[j] << ",";
                    else outfile << Bestgeti.X[j] << ",";
                }
            }
            cout<<"****************迭代结束！****************"<<endl;
            for (int j = 0; j < dimension; j++)
            {
                if (j == dimension - 1) std::cout << Bestgeti.X[j] << ")" << "对应的适应度为：" << Evaluate(Bestgeti) << endl;
                else if (j == 0) std::cout <<  Max_iter<<"次迭代后最终得到的最优鲸鱼个体(" << Bestgeti.X[j] << ",";
                else std::cout << Bestgeti.X[j] << ",";
            }

            outfile<<"****************迭代结束！****************"<<endl;
            for (int j = 0; j < dimension; j++)
            {
                if (j == dimension - 1) outfile << Bestgeti.X[j] << ")" << "对应的适应度为：" << Evaluate(Bestgeti) << endl;
                else if (j == 0) outfile <<  Max_iter<<"次迭代后最终得到的最优鲸鱼个体(" << Bestgeti.X[j] << ",";
                else outfile << Bestgeti.X[j] << ",";
            }
            outfile.close();
        }
};

int main()
{
    srand((int)time(0));  // 产生随机种子，否则每次的随机结果都是一样
	std::cout << "****************鲸鱼优化算法求解旅行商问题！****************" << endl;
	WOA woa;
	woa.WOA_TSP(50,29,1,29,100,"D:\\program files\\CodeBlocks\\myworkspace\\WOA_TSP\\bayg29.tsp");
	return 0;
}
