#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <string>
#include <time.h>

using namespace std;
//-----------------------------------------------------------------------------------------------------------//
/////////////////////////预处理部分////////////////////////
//-----------------------------------------------------------------------------------------------------------//
#define PI 3.1415926535
#define u0 4*PI*1e-7
#define G0 6.67259*1e-11
//-----------------------------------------------------------------------------------------------------------//
//=====================正演函数定义部分=====================================================================//
//-----------------------------------------------------------------------------------------------------------//
//>>>>>>>>部分正演参数定义
double F_I;
double F_A;
double F_M;
double T_M;
double* B_T;
double* F_T;
double* IF_T;
double* mb;
double* G;
int xstart, xend, ystart, yend, zstart, zend;
double cx, cy, cz;
double x = 0, y = 0, z = 0;
double* Fm0;//拟合模型
double* Fdate;//拟合数据
//重磁异常正演
double ForwardT(double k1, double k2, double k3, double k4, double k5, double k6,
                double tx, double ty, double tz,
                double x1, double x2, double y1, double y2, double z1, double z2);
double ForwardG(double tx, double ty, double tz,
                double x1, double x2, double y1, double y2, double z1, double z2);
//核函数获取
void precp_getG(double F_I, double F_A,
                double tx, double ty, double tz,
                double dx, double dy, double dz,
                double* cp,
                double* G, int points, int choice);
void test(double* test, int testlong);
//-----------------------------------------------------------------------------------------------------------//
//=====================反演函数定义部分======================//
//-----------------------------------------------------------------------------------------------------------//
//部分反演参数定义
double* m0;
double* m_get;
double* G_T;
double* G_TG;
double* GG_T;
double eps;
double m_max, m_min;
double w_up, w_down;
int flag;
double* wd_vector, * wm_vector;
double miu0;
////>>>>>>>>>>>>>>反演函数>>>>>>>>>>>>>>>>>//
void inversion(double* m0, double* mb, double miu0, double eps, int cd, int wg,
               double m_min, double m_max, double* d_tra, double* F_T, double* B_T, double* wm_vector, double* wd_vector,
               int flag, double* A);
////>>>>>>>>>>>>>>矩阵相关函数>>>>>>>>>>>>>>>>>//
void array_diag(double* rarry, double* getdiag, int numlong);
//>>>>>>>>>>>>>>向量扩展对角矩阵>>>>>>>>>>>>>>>>>//
void vectodiag(double* avect, int vectlong, double* neediag);
//>>>>>>>>>>>>>>矩阵乘法>>>>>>>>>>>>>>>>>>>>>>>>>//
//void array_dot(double* ary1, int line1, int ray1, double* ary2, int line2, int ray2, double* dotary);
void array_dot(double* matrix_1, int line1, int column1, double* matrix_2, int line2, int column2, double* output_matrix);
//>>>>>>>>>>>>>>矩阵数乘>>>>>>>>>>>>>>>>>>>>>>>>>//
void array_multiply(double chengshu, double* ary3, int ary3long, double* ary4);
//>>>>>>>>>>>>>>矩阵减法>>>>>>>>>>>>>>>>>>>>>>>>>//
void array_minus(double* ary4, double* ary5, int ary4long, double* aryminus);
//>>>>>>>>>>>>>>矩阵加法>>>>>>>>>>>>>>>>>>>>>>>>>//
void array_add(double* ary4, double* ary5, int ary4long, double* aryminus);
//>>>>>>>>>>>>>>向量指数>>>>>>>>>>>>>>>>>>>>>>>>>//
void vector_index(double* vecindex, double* getindex, int veclong1, double index);
//>>>>>>>>>>>>>>向量点积>>>>>>>>>>>>>>>>>>>>>>>>>//
double vector_dot(double* vec1, double* vec2, int veclong);
//>>>>>>>>>>>>>>矩阵转置>>>>>>>>>>>>>>>>>>>>>>>>>//
void array_tpose(double* getas, double* getat, int getax, int getay);
//>>>>>>>>>>>>>>向量相乘>>>>>>>>>>>>>>>>>>>>>>>>>//
void vector_multiply(double* vector1, double* vector2, int velong, double* vector3);
//>>>>>>>>>>>>对角矩阵向量算法>>>>>>>>>>>>>>>>>>>//
void array_diagdotleft(double* yuan, double* ve, int velong, int yuanxlong, int yuanylong, double* gets);
void array_diagdotright(double* yuan, double* ve, int velong, int yuanxlong, int yuanylong, double* gets);
///>>>>>>>>>>>>> amp加速下矩阵乘法>>>>>>>>>>>>>>>>>>>>>>>>>//
void MultiplyWithAMP(double* aMatrix, double* bMatrix, double* productMatrix, int row, int inner, int col);
//-----------------------------------------------------------------------------------------------------------//
//=======================================================//
/////////////////////////主函数部分////////////////////////
//=======================================================//
//-----------------------------------------------------------------------------------------------------------//

int main()
{
    clock_t start, end;
    start = clock();
    // ======================共用参数=====================//
    double x_step;//测点测线间距
    double y_step;
    int divide_x;
    int divide_y;   //网格剖分数
    int divide_z;
    double Xscale;
    double Yscale;  // 网格剖分尺度
    double Zscale;
    int Lines; // 测线数
    int Point; // 测点数
    double* zcp;
    double* cp;//存储角点坐标
    int xy_get, xz_get, yz_get;
    double *x_save;
    double *y_save;
    F_A = 0.;
    F_I = 0.;
    //=======================================================//
    /////////////////////////正演计算部分////////////////////////
    //=======================================================//
    int upto;
    int autoget;
    int amount;
    int cd;
    int wg;
    int choice;
    int unit;
    cout << "*========================================================================================================*\n";
    cout << "*-----------------------------------------&聚焦重磁三维正反演程序&----------------------------------------*\n";
    cout << "*=============================================&成都理工大学&=============================================*\n";
    cout << "请仔细阅读以下注意事项：\n";
    cout << "（1）输入不同参数需用空格分开\n";
    cout << "（2）8G运存电脑最大剖分网格数尽量小于5000  5000左右每次迭代次数时间>=6s\n\n";
    cout << "&=========================================&正演程序说明&=================================================&\n\n";
    cout << "（1）正演测点坐标由（0,0）开始向正方向取，单个剖分长方体xyz中心位置代表该剖分单元坐标\n";
    cout << "（2）正演数据以“正演背景磁异常.txt、正演长方体磁异常.txt、正演长方体叠加异常.txt、正演模型.txt”命名存放于根目录下\n";
    cout << "（3）正演长方体坐标自动输入文件以“正演长方体坐标.txt放于根目录下，每一行格式为 x0 x1 y0 y1 z0 z1\n\n";
    cout << "&=========================================&反演程序说明&=================================================&\n\n";
    cout << "（1）反演数据以“反演数据.txt”命名存放于根目录下，数据分三列，第一列为X，第二列为Y，第三列为磁异常值\n";
    cout << "（2）请根据具体数据输入合适的dx dy dz大小\n（3）背景磁异常值将以“背景场值.txt”输出于根目录下\n";
    cout << "（4）物性上限-物性下限的差值最好在推测地质体磁化强度左右浮动\n";
    cout << "（5）聚焦因子一般取1e-6  正则化参数建议取0.0001\n";
    cout << "（6）输出剖面将以“剖面xy.txt 剖面xz.txt 剖面yz.txt”的形式输出于根目录下\n";
    cout << "（7）反演模型将以“反演模型.txt”的形式输出于根目录下\n";
    cout << "（8）反演数据与原数据拟合差以“拟合差.txt”的形式输出于根目录下\n";
    cout << "*========================================================================================================*\n";
    cout << "请选择单位（1）m （2）km\n";
    cin >> unit;
    cout << "请选择正演或反演程序：正演程序（1） 反演程序（2）\n";
    cin >> upto;
    if (upto == 1)
    {
        cout << "请选择:重力正演（1） 磁正演（2）\n";
        cin >> choice;
        if (choice == 2)
        {
            cout << "输入磁倾角\n";
            cin >> F_I;
            cout << "请输入磁偏角\n";
            cin >> F_A;
        }
        cout << "请输入网格剖分数 xdivide  ydivide zdivide\n";
        cin >> divide_x;
        cin >> divide_y;
        cin >> divide_z;
        cout << "请输入网格剖分尺度dx dy dz\n";
        cin >> Xscale;
        cin >> Yscale;
        cin >> Zscale;
        cout << "请输入测点数X\n";
        cin >> Point;
        cout << "请输入测线数Y\n";
        cin >> Lines;
        cout << "请输入测点测线间距\n";
        cin >> x_step;
        cin >> y_step;
        if (unit == 2)
        {
            Xscale = Xscale * 1000;
            Yscale = Yscale * 1000;
            Zscale = Zscale * 1000;
            x_step = x_step * 1000;
            y_step = y_step * 1000;
        }
        cout << "请输入背景场强度\n";
        cin >> F_M;
        cout << "请输入正演的长方体个数\n";
        cin >> amount;
        cp = new double[amount * 6];
        if (amount)
        {
            cout << "设置各长方体角点坐标 1.手动输入\t2.自动输入\n";
            cin >> autoget;
            switch (autoget)
            {
                case 1:
                    for (int i = 0; i < amount; i++)
                    {
                        cout << "请输入第" << i + 1 << "个长方形角点坐标 x0 x1 y0 y1 z0 z1（以空格分割）\n";
                        for (int j = 0; j < 6; j++)
                        {
                            cin >> cp[i * 6 + j];
                        }
                    }
                case 2:
                    ifstream infile("正演长方体坐标.txt");
                    string getcft;
                    for (int i = 0; i <amount; i++)//
                    {
                        getline(infile, getcft);
                        stringstream stream;
                        stream << getcft;
                        for (int getc = 0; getc < 6; getc++)
                        {

                            stream >> cp[i * 6 + getc];
                        }
                    }
            }
            //test(cp, amount * 6);
            if(choice == 1)
            {
                cout << "请输入长方体剩余密度\n";
            }
            else if (choice==2)
            {
                cout << "请输入长方体磁化强度\n";
            }
            cin >> T_M;
        }
        //======================================================================正演参数设定
        cd = Point * Lines;
        wg = divide_x * divide_y * divide_z;
        mb = new double[wg];
        for (int reset = 0; reset < wg; reset++)
        {
            mb[reset] = 1 * F_M;
        }
        double* mF;
        mF = new double[wg];
        memset(mF, 0, sizeof(double) * wg);
        zcp = new double[6];
        zcp[0] = 0;
        zcp[2] = 0;
        zcp[4] = 0;
        zcp[1] = divide_x * Xscale;
        zcp[3] = divide_y * Yscale;
        zcp[5] = divide_z * Zscale;
        ofstream outfileB_T("正演背景异常.txt");
        ofstream outfileT("正演长方体异常.txt");
        ofstream oufileIT("正演长方体叠加异常.txt");
        B_T = new double[Point * Lines];
        F_T = new double[Point * Lines];
        memset(B_T, 0, sizeof(double) * cd);
        memset(F_T, 0, sizeof(double) * cd);
        IF_T = new double[Point * Lines];
        memset(IF_T, 0, sizeof(double) * cd);
        G = new double[cd * wg];
        //======================================================================正演计算
        int save_get = 0;
        x_save = new double[cd];
        y_save = new double[cd];
        for (int xFi = 0; xFi < Point; xFi++)
        {
            for (int yFi = 0; yFi < Lines; yFi++)
            {
                x = x_step * xFi;
                y = y_step * yFi;
                x_save[save_get] = x;
                y_save[save_get] = y;
                precp_getG(F_I, F_A, x, y, 0, Xscale, Yscale, Zscale, zcp, G, save_get,choice);
                save_get += 1;
            }
        }
        //test(G, cd * wg);
        //MultiplyWithAMP(G, mb, B_T, cd, wg, 1);
        //cout << "if finish";
        array_dot(G, cd, wg, mb, wg, 1, B_T);//背景场计算
        //test(B_T,cd);
        //cout << "yes";
        //=====================================================================根据角点坐标计算模型参数
        for (int cf = 0; cf < amount; cf++)
        {
            xstart = (int)cp[cf * 6] / Xscale;
            xend = (int)cp[cf * 6 + 1] / Xscale;
            ystart = (int)cp[cf * 6 + 2] / Yscale;
            yend = (int)cp[cf * 6 + 3] / Yscale;
            zstart = (int)cp[cf * 6 + 4] / Zscale;
            zend = (int)cp[cf * 6 + 5] / Zscale;
            for (int cfx = xstart; cfx <= xend; cfx++)
            {
                for (int cfy = ystart; cfy <= yend; cfy++)
                {
                    for (int cfz = zstart; cfz <= zend; cfz++)
                    {
                        mF[cfx * divide_y * divide_z + cfy * divide_z + cfz] = T_M;
                        //cout << mF[cfx * divide_y * divide_z + cfy * divide_z + cfz];
                    }
                }
            }
        }
        array_dot(G, cd, wg, mF, wg, 1, F_T);//正演长方体计算
        array_add(F_T, B_T, cd, IF_T);//叠加背景场
        outfileB_T << "x\ty\tT\n";
        outfileT << "x\ty\tT\n";
        oufileIT << "x\ty\tT\n";
        for (int i = 0; i < cd; i++)
        {

            outfileB_T << x_save[i] << "\t" << y_save[i] << "\t" << B_T[i] << endl;
            outfileT << x_save[i] << "\t" << y_save[i] << "\t" << F_T[i] << endl;
            oufileIT << x_save[i] << "\t" << y_save[i] << "\t" << IF_T[i] << endl;
        }
        outfileB_T.close();
        oufileIT.close();
        outfileT.close();
    }
    else
    {
        //-----------------------------------------------------------------------------------------------------------//
        //================================反演计算部分============================
        //-----------------------------------------------------------------------------------------------------------//
        cout << "请选择：重力反演（1） 磁反演（2）\n";
        cin >> choice;
        cout << "请输入网格剖分数 xdivide  ydivide zdivide\n";
        cin >> divide_x;
        cin >> divide_y;
        cin >> divide_z;
        cout << "请输入网格剖分尺度dx dy dz\n";
        cin >> Xscale;
        cin >> Yscale;
        cin >> Zscale;
        cout << "请输入测点数X\n";
        cin >> cd;
        if (unit == 2)
        {
            Xscale = Xscale * 1000;
            Yscale = Yscale * 1000;
            Zscale = Zscale * 1000;
        }
        cout << "请输入聚焦因子\n";
        cin >> eps;
        cout << "请输入正则化参数\n";
        cin >> miu0;
        cout << "请输入背景场强度\n";
        cin >> F_M;
        cout << "请输入物性上限\n";
        cin >> w_up;
        cout << "请输入物性下限\n";
        cin >> w_down;
        cout << "请输入迭代次数\n";
        cin >> flag;
        int iauto;
        cout << "请输入需要的剖面 1.自定义设置 2.中心剖面\n";
        cin >> iauto;
        switch (iauto)
        {
            case 1:
                cout << "请输入需要第几个xy剖面\n";
                cin >> xy_get;
                cout << "请输入需要第几个yz剖面\n";
                cin >> yz_get;
                cout << "请输入需要第几个xz剖面\n";
                cin >> xz_get;
            case 2:
                xy_get = divide_z / 2;
                xz_get = divide_y / 2;
                yz_get = divide_x / 2;
        }
        //======================================================================反演参数设定
        double* X_get, * Y_get;
        m_max = F_M + w_up;
        m_min = F_M + w_down;
        wg = divide_x * divide_y * divide_z;
        F_T = new double[cd];
        B_T = new double[cd];
        IF_T = new double[cd];
        X_get = new double[cd];
        Y_get = new double[cd];
        zcp = new double[6];
        zcp[0] = 0;
        zcp[2] = 0;
        zcp[4] = 0;
        zcp[1] = divide_x * Xscale;
        zcp[3] = divide_y * Yscale;
        zcp[5] = divide_z * Zscale;
        //======================================================================反演数据智能录入
        //数据及对应坐标录入
        ifstream infileget("./double_cf.txt");
        string str1;
        for (int fi = 0; fi <= cd; fi++)//
        {
            getline(infileget, str1);
            int get_x = 0, get_y = 0;
            double c_x = 0., c_y = 0.;
            stringstream stream1;
            if (fi == 0)
            {
                stream1 << str1;
            }
            else
            {
                stream1 << str1;
                stream1 >> c_x >> c_y;
                X_get[fi - 1] = c_x;
                Y_get[fi - 1] = c_y;
                stream1 >> F_T[fi-1];
            }
        }
        infileget.close();
        //======================================================================核矩阵计算+背景场+叠加场计算
        mb = new double[wg];
        for (int reset = 0; reset < wg; reset++)
        {
            mb[reset] = 1 * F_M;
        }
        G = new double[cd * wg];
        for (int i = 0; i < cd; i++)
        {
            precp_getG(90, 0, X_get[i], Y_get[i], 0, Xscale, Yscale, Zscale, zcp, G,i,choice);
        }
        array_dot(G, cd, wg, mb, wg, 1, B_T);
        array_add(B_T, F_T, cd, IF_T);
        //======================================================================反演
        wm_vector = new double[wg];
        wd_vector = new double[cd];
        m0 = new double[wg];
        G_T = new double[cd * wg];
        G_TG = new double[wg * wg];
        memset(G_TG, 0, sizeof(double) * wg * wg);
        GG_T = new double[cd * cd];
        memset(GG_T, 0, sizeof(double) * cd * cd);
        for (int reset = 0; reset < wg; reset++)
        {
            m0[reset] = mb[reset];
        }
        array_tpose(G, G_T, cd, wg);
        MultiplyWithAMP(G, G_T, GG_T, cd, wg, cd);
        //array_dot(G, cd, wg, G_T, wg, cd, GG_T);
        MultiplyWithAMP(G_T, G, G_TG, wg, cd, wg);
        //array_dot(G_T, wg, cd, G, cd, wg, G_TG);
        array_diag(GG_T, wd_vector, cd);
        array_diag(G_TG, wm_vector, wg);
        vector_index(wd_vector, wd_vector, cd, 0.5);
        vector_index(wm_vector, wm_vector, wg, 0.5);
        inversion(m0, mb, miu0, eps, cd, wg, m_min, m_max, IF_T, F_T, B_T, wm_vector, wd_vector, flag, G);
        //=====================================================================减去背景模型
        m_get = new double[wg];
        array_minus(m0, mb, wg, m_get);
        //======================================================================背景场输出
        ofstream btget("背景场值.txt");
        btget << "x\ty\tB_T\n";
        for (int bi = 0; bi < cd; bi++)
        {
            btget << X_get[bi] << "\t" << Y_get[bi] << "\t" << B_T[bi] << endl;
        }
        //======================================================================反演数据输出
        //======================================================================三剖面数据输出
        ofstream xyget("xy剖面.txt");
        xyget << "x\ty\tm_T\n";
        for (int xd = 0; xd < divide_x; xd++)
        {
            for (int yd = 0; yd < divide_y; yd++)
            {

                xyget << xd * Xscale << "\t" << yd * Yscale << "\t" << m_get[xd * divide_y * divide_z + yd * divide_z + xy_get] << "\n";
            }
        }
        ofstream xzget("xz剖面.txt");
        xzget << "x\tz\tm_T\n";
        for (int xd = 0; xd < divide_x; xd++)
        {
            for (int zd = 0; zd < divide_z; zd++)
            {

                xzget << xd * Xscale << "\t" << zd * Zscale << "\t" << m_get[xd * divide_y * divide_z + xz_get * divide_z + zd] << "\n";
            }
        }
        ofstream yzget("yz剖面.txt");
        yzget << "y\tz\tm_T\n";
        for (int yd = 0; yd < divide_y; yd++)
        {
            for (int zd = 0; zd < divide_z; zd++)
            {

                yzget << yd * Yscale << "\t" << zd * Zscale << "\t" << m_get[yz_get * divide_y * divide_z + yd * divide_z + zd] << "\n";
            }
        }
        xyget.close();
        xzget.close();
        yzget.close();
        //======================================================================z拟合差输出
        ofstream rmsget("拟合差.txt");
        rmsget << "测点\t原数据\t拟合数据\t拟合差\n";
        Fdate = new double[cd];
        Fm0 = new double[wg];
        array_minus(m0, mb, wg, Fm0);
        array_dot(G, cd, wg, Fm0, wg, 1, Fdate);
        for(int i2 =0;i2<cd;i2++)
        {
            rmsget <<i2<<"\t"<< F_T[i2] << "\t" << Fdate[i2] << "\t" << F_T[i2] - Fdate[i2] << "\n";
        }
        rmsget.close();
        ofstream Fmo("反演模型.txt");
        Fmo << "x\ty\tz\taspen\n";
        for (int ix = 0; ix < divide_x; ix++)
        {
            for (int iy = 0; iy < divide_y; iy++)
            {
                for (int iz = 0; iz < divide_z; iz++)
                {
                    Fmo << ix * Xscale << "\t" << iy * Yscale << "\t" << iz * Zscale << "\t" << Fm0[ix * divide_y * divide_z + iy * divide_z + iz] << "\n";
                }
            }

        }
        Fmo.close();
    }
    end = clock();
    cout<<"运行时间"<<(double)(end-start)/CLOCKS_PER_SEC<<"s"<<endl;
    return 0;
}
/*F_I = 90;
F_A = 0;
divide_x = 11;
divide_y = 11;
divide_z = 11;
Xscale = 20;
Yscale = 20;
Zscale = 20;
Point = 10;
Lines = 10;
F_M = 1;
eps = 1e-6;
m_max = 2;
m_min = 1;
flag = 100;*/
//=======================================================//
/////////////////////////函数实现部分//////////////////////
//=======================================================//
//=======================正演部分=======================//
//>>>>>>>>>>>>>>1.预处理及核矩阵获取>>>>>>>>>>>>>>>>>//
void precp_getG(double F_I, double F_A,
                double tx, double ty, double tz,
                double dx, double dy, double dz,
                double* cp,
                double* G, int points, int choice)
{
    int divide_x, divide_y, divide_z;
    double sum = 0;
    double trans = 0;
    double L = cos(F_I * PI / 180.0) * cos(F_A * PI / 180.0);
    double M = cos(F_I * PI / 180.0) * sin(F_A * PI / 180.0);
    double N = sin(F_I * PI / 180.0);
    double k1 = 2 * M * N;
    double k2 = 2 * N * L;
    double k3 = 2 * M * L;
    double k4 = L * L;
    double k5 = M * M;
    double k6 = -1 * N * N;
    int i, j, k, p;
    double x0 = cp[0];
    double xn = cp[1];
    double y0 = cp[2];
    double yn = cp[3];
    double z0 = cp[4];
    double zn = cp[5];
    double x1, x2, y1, y2, z1, z2;
    divide_x = (int)(xn - x0) / dx;
    divide_y = (int)(yn - y0) / dy;
    divide_z = (int)(zn - z0) / dz;
    if (choice == 1)
    {
        for (i = 0; i < divide_x; i++)
        {
            x1 = x0 + (i - 0.5) * dx;
            x2 = x0 + (i + 0.5) * dx;
            for (j = 0; j < divide_y; j++)
            {
                y1 = y0 + (j - 0.5) * dy;
                y2 = y0 + (j + 0.5) * dy;
                for (k = 0; k < divide_z; k++)
                {
                    z1 = z0 + k * dz;
                    z2 = z1 + dz;
                    if (z1 == 0)
                        trans = G0 * ForwardG(tx, ty, dz / 1000, x1, x2, y1, y2, z1, z2) ;
                    else
                    {
                        trans = G0 * ForwardG(tx, ty, 0, x1, x2, y1, y2, z1, z2);
                    }
                    G[points * divide_x * divide_y * divide_z  + i * divide_y * divide_z +j * divide_z + k] = trans;

                }
            }
        }
    }
    else if (choice == 2)
    {
        for (i = 0; i < divide_x; i++)
        {
            x1 = x0 + (i - 0.5) * dx;
            x2 = x0 + (i + 0.5) * dx;
            for (j = 0; j < divide_y; j++)
            {
                y1 = y0 + (j - 0.5) * dy;
                y2 = y0 + (j + 0.5) * dy;
                for (k = 0; k < divide_z; k++)
                {
                    z1 = z0 + k * dz;
                    z2 = z1 + dz;
                    if (z1 == 0)
                        trans = u0 / (4 * PI) * ForwardT(k1, k2, k3, k4, k5, k6, tx, ty, dz / 1000, x1, x2, y1, y2, z1, z2) * 1e9;
                    else
                    {
                        trans = u0 / (4 * PI) * ForwardT(k1, k2, k3, k4, k5, k6, tx, ty, 0, x1, x2, y1, y2, z1, z2) * 1e9;
                    }
                    G[points * divide_x * divide_y * divide_z  + i * divide_y * divide_z +j * divide_z + k] = trans;

                }
            }
        }
    }
}
//>>>>>>>>>>>>>>2.1剖分单元长方体磁异常值计算>>>>>>>>>>>>>>>>>>>>//
double ForwardT(double k1, double k2, double k3, double k4, double k5, double k6,
                double tx, double ty, double tz,
                double x1, double x2, double y1, double y2, double z1, double z2)
{
    int i, j, k;
    double T = 0;
    double R = 0;
    double dltx, dlty, dltz;
    double x[] = { x1,x2 };
    double y[] = { y1,y2 };
    double z[] = { z1,z2 };
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            for (k = 0; k < 2; k++)
            {
                dltx = x[i] - tx;
                dlty = y[j] - ty;
                dltz = z[k] - tz;
                R = sqrt(dltx * dltx + dlty * dlty + dltz * dltz);
                T += pow(-1, i + j + k + 1) * (k1 * log(R + dltx) + k2 * log(R + dlty) + k3 * log(R + dltz) +
                                               k4 * atan2((dltx * dlty), (dltx * dltx + R * dltz + dltz * dltz)) +
                                               k5 * atan2((dltx * dlty), (dlty * dlty + R * dltz + dltz * dltz)) + k6 * atan2((dltx * dlty), (R * dltz)));
            }
        }
    }
    return T;
}
//>>>>>>>>>>>>>>2.2剖分单元长方体重力异常值计算>>>>>>>>>>>>>>>>>>>>//
double ForwardG(double tx, double ty, double tz,
                double x1, double x2, double y1, double y2, double z1, double z2)
{
    int i, j, k;
    double T = 0;
    double R = 0;
    double dltx, dlty, dltz;
    double x[] = { x1,x2 };
    double y[] = { y1,y2 };
    double z[] = { z1,z2 };
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            for (k = 0; k < 2; k++)
            {
                dltx = x[i] - tx;
                dlty = y[j] - ty;
                dltz = z[k] - tz;
                R = sqrt(dltx * dltx + dlty * dlty + dltz * dltz);
                T += 1e9 * pow(-1, i + j + k) * (dltx * log(dlty + R) + dlty * log(dltx + R) - dltz * atan2((dltx * dlty), (dltz * R)));
            }
        }
    }
    return T;
}
//=======================反演部分=======================//
void get_we(double* m0, int m0long, double eps, double* regetwe)
{
    for (int mi = 0; mi < m0long; mi++)
    {
        double tr = m0[mi] * m0[mi] + eps * eps;
        regetwe[mi] = 1 / pow(tr, 0.5);
    }
}
void inversion(double* m0, double* mb, double miu0, double eps, int cd, int wg,
               double m_min, double m_max, double* d_tra, double* F_T, double* B_T, double* wm_vector, double* wd_vector,
               int flag, double* A)
{
    int swit = 0;
    double* we_vector, * we_vector1, * we_vector2;
    //=======================================
    double* perfectwgwg;
    perfectwgwg = new double[wg * wg];
    double* perfectwgwg1;
    perfectwgwg1 = new double[wg * wg];
    double* perfectwgwg2;
    perfectwgwg2 = new double[wg * wg];
    double* perfectwg;
    perfectwg = new double[wg];
    double* perfectwg1;
    perfectwg1 = new double[wg];
    double* perfectcdcd;
    perfectcdcd = new double[cd * cd];
    double* perfectcd;
    perfectcd = new double[cd];
    double* perfectcdwg;
    perfectcdwg = new double[cd * wg];
    double* perfectcdwg1;
    perfectcdwg1 = new double[cd * wg];
    double cyx = 0, cyx1 = 0;
    //========================================
    double* wm_1vector;
    double* mw;
    double* Aw, * Aw_T;
    double* Q;
    double* q;
    double* f1;
    double* f0;
    double* d0;
    double t;
    double beta;
    double* m1, * T1, * T2;
    double up = 0, down = 0;
    double* fd;
    double* delta_d;
    double rms;
    double faid = 0;
    fd = new double[cd];
    delta_d = new double[cd];
    //========================================
    we_vector = new double[wg];
    we_vector1 = new double[wg];
    we_vector2 = new double[wg];
    mw = new double[wg];
    memset(mw, 0, sizeof(double) * wg);
    wm_1vector = new double[wg];
    double* dw;
    dw = new double[cd];
    f0 = new double[wg];
    Aw = new double[cd * wg];
    memset(Aw, 0, sizeof(double) * cd * wg);
    Q = new double[wg * wg];
    memset(Q, 0, sizeof(double) * wg * wg);
    Aw_T = new double[cd * wg];
    memset(Aw_T, 0, sizeof(double) * cd * wg);
    q = new double[wg];
    memset(q, 0, sizeof(double) * wg);
    f1 = new double[wg];
    memset(f1, 0, sizeof(double) * wg);
    d0 = new double[wg];
    memset(d0, 0, sizeof(double) * wg);
    m1 = new double[wg];
    memset(m1, 0, sizeof(double) * wg);
    //======================================
    //dw
    vector_multiply(wd_vector, d_tra, cd, dw);
    //test(dw, cd);
    //vectodiag(wm_vector, wg, wm);
    //vectodiag(wd_vector, cd, wd);
    vector_index(wm_vector, wm_1vector, wg, -1);
    //==========================================================================test(wd_vector, wg);
    //vectodiag(wm_1vector, wg, wm_1);
    int i = 0;
    while (i<flag)
    {
        if (i == 0)
        {
            get_we(m0, wg, eps, we_vector);
            //vectodiag(we_vector, wg, we);
            //==================================================================test(we_vector, wg);
            vector_index(we_vector, we_vector1, wg, -1);
            vector_index(we_vector, we_vector2, wg, -2);
            //vectodiag(we_vector1, wg, we_1);
            //vectodiag(we_vector2, wg, we_2);
            //mw
            array_diagdotleft(m0, wm_vector, wg, wg, 1, perfectwg);
            //====================================================================test(m0, wg);
            //array_dot(wm, wg, wg, m0, wg, 1, perfectwg);
            array_diagdotleft(perfectwg, we_vector, wg, wg, 1, mw);
            //array_dot(we, wg, wg, perfectwg, wg, 1, mw);
            //====================================================================test(mw, wg);
            //Aw
            array_diagdotleft(G, wd_vector, cd, cd, wg, perfectcdwg);
            //array_dot(wd, cd, cd, G,cd, wg, perfectcdwg);
            array_diagdotright(perfectcdwg, wm_1vector, wg, cd, wg, perfectcdwg1);
            //array_dot(perfectcdwg1, cd, wg, wm_1, wg, wg, perfectcdwg1);
            array_diagdotright(perfectcdwg1, we_vector1, wg, cd, wg, Aw);
            //array_dot(perfectcdwg1, cd, wg, we_1, wg, wg, Aw);
            //====================================================================test(Aw, cd*wg);
            //Q
            for (int qi = 0; qi < wg; qi++)
            {
                perfectwg[qi] = miu0 * eps * eps;
            }
            vectodiag(perfectwg, wg, perfectwgwg);
            array_tpose(Aw, Aw_T, cd, wg);
            MultiplyWithAMP(Aw_T, Aw,perfectwgwg1, wg, cd, wg);
            //array_dot(Aw_T, wg, cd, Aw, cd, wg, perfectwgwg1);
            array_diagdotleft(perfectwgwg1, we_vector2, wg, wg, wg, perfectwgwg2);
            //array_dot(we_2, wg, wg, perfectwgwg1, wg, wg, perfectwgwg2);
            array_add(perfectwgwg, perfectwgwg2, wg * wg, Q);
            //===================================================================test(Q, wg * wg);
            //q
            array_dot(Aw_T, wg, cd, dw, cd, 1, perfectwg);
            array_diagdotleft(perfectwg, we_vector2, wg, wg, 1, q);
            //===================================================================test(q, wg);
            //array_dot(we_2, wg, wg, perfectwg, wg, 1, q);
            //f1
            array_dot(Q, wg, wg, mw, wg, 1, perfectwg);
            array_minus(perfectwg, q, wg, f1);
            //===================================================================test(f1, wg);
            //d0
            array_multiply(-1, f1, wg, d0);
            //===================================================================test(d0, wg);
            //t
            cyx =vector_dot(d0, f1, wg);
            //cout<<"warning" << cyx;
            array_dot(Q, wg, wg, d0, wg, 1, perfectwg);
            cyx1 = vector_dot(d0, perfectwg, wg);
            t = cyx / cyx1;
            //cout << "\nthis"<<t;
            //mw
            array_multiply(t, d0, wg, perfectwg);
            array_minus(mw, perfectwg, wg, perfectwg1);
            for (int mwi = 0; mwi < wg; mwi++)
            {
                mw[mwi] = perfectwg1[mwi];
            }
            //===================================================test(mw, wg);
            //m0
            array_diagdotleft(mw, we_vector1, wg, wg, 1, perfectwg);
            //array_dot(we_1, wg, wg, mw, wg, 1, perfectwg);
            array_diagdotleft(perfectwg, wm_1vector, wg, wg, 1, m0);
            //array_dot(wm_1, wg, wg, perfectwg, wg, 1, m0);
            for (int die = 0; die < wg; die++)
            {
                if (m0[die] < m_min) { m0[die] = m_min; }
                else if (m0[die] > m_max) { m0[die] = m_max; }
            }
            //====================================================test(m0, wg);
            //we
            get_we(m0, wg, eps, we_vector);
            //vectodiag(we_vector, wg, we);
            //mw
            array_diagdotleft(m0, wm_vector, wg, wg, 1, perfectwg);
            //array_dot(wm, wg, wg, m0, wg, 1, perfectwg);
            array_diagdotleft(perfectwg, we_vector, wg, wg, 1, mw);
            //array_dot(we, wg, wg, perfectwg, wg, 1, mw);
            //test(mw, wg);
        }
        else
        {
            //we
            vector_index(we_vector, we_vector1, wg, -1);
            vector_index(we_vector, we_vector2, wg, -2);
            //vectodiag(we_vector1, wg, we_1);
            //vectodiag(we_vector2, wg, we_2);
            //Aw
            array_diagdotleft(G, wd_vector, cd, cd, wg, perfectcdwg);
            //array_dot(wd, cd, cd, G,cd, wg, perfectcdwg);
            array_diagdotright(perfectcdwg, wm_1vector, wg, cd, wg, perfectcdwg1);
            //array_dot(perfectcdwg1, cd, wg, wm_1, wg, wg, perfectcdwg1);
            array_diagdotright(perfectcdwg1, we_vector1, wg, cd, wg, Aw);
            //array_dot(perfectcdwg1, cd, wg, we_1, wg, wg, Aw);
            //====================================================================test(Aw, cd*wg);
            //Q
            for (int qi = 0; qi < wg; qi++)
            {
                perfectwg[qi] = miu0 * eps * eps;
            }
            vectodiag(perfectwg, wg, perfectwgwg);
            array_tpose(Aw, Aw_T, cd, wg);
            MultiplyWithAMP(Aw_T, Aw, perfectwgwg1, wg, cd, wg);
            //array_dot(Aw_T, wg, cd, Aw, cd, wg, perfectwgwg1);
            array_diagdotleft(perfectwgwg1, we_vector2, wg, wg, wg, perfectwgwg2);
            //array_dot(we_2, wg, wg, perfectwgwg1, wg, wg, perfectwgwg2);
            array_add(perfectwgwg, perfectwgwg2, wg * wg, Q);
            //===================================================================test(Q, wg * wg);
            //q
            array_dot(Aw_T, wg, cd, dw, cd, 1, perfectwg);
            array_diagdotleft(perfectwg, we_vector2, wg, wg, 1, q);
            //===================================================================test(q, wg);
            //array_dot(we_2, wg, wg, perfectwg, wg, 1, q);
            //f0
            for (int cx = 0; cx < wg; cx++)
            {
                f0[cx] = f1[cx];
            }
            //f1
            array_dot(Q, wg, wg, mw, wg, 1, perfectwg);
            array_minus(perfectwg, q, wg, f1);
            //beta
            cyx = vector_dot(f1, f1, wg);
            cyx1=vector_dot(f0, f0, wg);
            beta = cyx / cyx1;
            //cout << beta;
            //d0
            array_multiply(-1, f1, wg, perfectwg);
            array_multiply(beta, d0, wg, perfectwg1);
            array_add(perfectwg, perfectwg1, wg, d0);
            //t
            cyx = vector_dot(d0, f1, wg);
            array_dot(Q, wg, wg, d0, wg, 1, perfectwg);
            cyx1 = vector_dot(d0, perfectwg, wg);
            t = cyx / cyx1;
            //mw
            array_multiply(t, d0, wg, perfectwg);
            array_minus(mw, perfectwg, wg, perfectwg1);
            for (int mwi1 = 0; mwi1 < wg; mwi1++)
            {
                mw[mwi1] = perfectwg1[mwi1];
            }
            //m0
            array_diagdotleft(mw, we_vector1, wg, wg, 1, perfectwg);
            //array_dot(we_1, wg, wg, mw, wg, 1, perfectwg);
            array_diagdotleft(perfectwg, wm_1vector, wg, wg, 1, m0);
            //array_dot(wm_1, wg, wg, perfectwg, wg, 1, m0);
            for (int die = 0; die < wg; die++)
            {
                if (m0[die] < m_min) { m0[die] = m_min; }
                else if (m0[die] > m_max) { m0[die] = m_max; }
            }
            //================================================test(m0, wg);
            //we
            get_we(m0, wg, eps, we_vector);
            //vectodiag(we_vector, wg, we);
            //mw
            array_diagdotleft(m0, wm_vector, wg, wg, 1, perfectwg);
            //array_dot(wm, wg, wg, m0, wg, 1, perfectwg);
            array_diagdotleft(perfectwg, we_vector, wg, wg, 1, mw);
            //array_dot(we, wg, wg, perfectwg, wg, 1, mw);
        }
        if (swit == 0)
        {

            swit = 1;
            T1 = new double[cd];
            T2 = new double[cd];
            array_minus(m0, mb, wg, m1);
            array_dot(G, cd, wg, m1, wg, 1, perfectcd);
            array_minus(perfectcd, F_T, cd, T1);
            array_diagdotleft(T1, wd_vector, cd, cd, 1, T2);
            //array_dot(wd, cd, cd, T1, cd, 1, T2);
            up =vector_dot(T2, T2, cd);
            down = vector_dot(mw, mw, wg);
            miu0 = up / down;
            //miu0 = 2000;
            for (int reset = 0; reset < wg; reset++)
            {
                m0[reset] = mb[reset];
            }
            i = 0;
            continue;
        }
        else
        {
            array_dot(G, cd, wg, m0, wg, 1, perfectcd);
            array_minus(perfectcd, d_tra, cd, fd);
            faid = vector_dot(fd, fd, cd);
            array_minus(m0, mb, wg, m1);
            array_dot(G, cd, wg, m1, wg, 1, perfectcd);
            //test(perfectcd,cd);
            array_minus(F_T, perfectcd, cd, delta_d);
            rms = 0;
            for (int pp = 0; pp < cd; pp++)
            {
                double rmsget = delta_d[pp] / F_T[pp];
                rms += pow(rmsget, 2);
            }
            rms = sqrt(rms / cd);
            cout  << "Times=" << i << setw(18) << "Faid=" << faid  <<setw(14)<< "miu=" << miu0 << setw(18) << "Rms=" << rms << "\n";
            miu0 = miu0 * 0.9;
            i += 1;
        }
    }
    delete[]perfectcd;
    delete[]perfectcdcd;
    delete[]perfectcdwg;
    delete[]perfectcdwg1;
    delete[]perfectwg;
    delete[]perfectwg1;
    delete[]perfectwgwg1;
    delete[]perfectwgwg2;
    delete[]perfectwgwg;
    delete[]Aw;
    delete[]Aw_T;
    delete[]Q;
    delete[]q;
    delete[]d0;
    delete[]mw;
    delete[]f0;
    delete[]f1;
    delete[]wm_1vector;
    delete[]we_vector;
    delete[]we_vector1;
    delete[]we_vector2;
}


//==================矩阵函数===========================//
//>>>>>>>>>>>>>>矩阵取对角元素>>>>>>>>>>>>>>>>>//
void array_diag(double* rarry, double* getdiag, int numlong)
{
    int cplong;
    cplong = numlong;
    for (int i = 0; i < numlong; i++)
    {
        for (int j = 0; j < cplong; j++)
        {
            if (i == j)
            {
                getdiag[i] = rarry[i * cplong + j];
            }
        }
    }

}
//>>>>>>>>>>>>>>向量扩展对角矩阵>>>>>>>>>>>>>>>>>//
void vectodiag(double* avect, int vectlong, double* neediag)
{
    for (int avt = 0; avt < vectlong; avt++)
    {
        for (int jvt = 0; jvt < vectlong; jvt++)
        {
            if (avt == jvt)
            {
                neediag[avt * vectlong + jvt] = avect[avt];
            }
            else {

                neediag[avt * vectlong + jvt] = 0;
            }
        }
    }

}
//>>>>>>>>>>>>>>矩阵乘法>>>>>>>>>>>>>>>>>>>>>>>>>//
/*void array_dot(double* ary1, int line1, int ray1, double* ary2, int line2, int ray2, double* dotary)
{
	double sum;
	for (int m = 0; m < line1; m++)
	{
		for (int n = 0; n < ray2; n++)
		{
			double sum = array_dotget(m, n, ary1, ray1, ary2, ray2);
			dotary[m * ray2 + n]= sum;
		}
	}
}
*/
void array_dot(double* matrix_1, int line1, int column1, double* matrix_2, int line2, int column2, double* output_matrix)
{
    memset(output_matrix, 0, sizeof(double) * line1 * column2);
    int i, j, k;
    double sum = 0;
    for (i = 0; i < line1; i++)
        for (k = 0; k < column1; k++)
        {
            sum = matrix_1[i * column1 + k];
            for (j = 0; j < column2; j++)

                output_matrix[i * column2 + j] += sum * matrix_2[k * column2 + j];
        }
}
//>>>>>>>>>>>>>>矩阵数乘>>>>>>>>>>>>>>>>>>>>>>>>>//
void array_multiply(double chengshu, double* ary3, int ary3long, double* ary4)
{
    for (int i = 0; i < ary3long; i++)
    {
        ary4[i] = chengshu * ary3[i];
    }

}
//>>>>>>>>>>>>>>矩阵减法>>>>>>>>>>>>>>>>>>>>>>>>>//
void array_minus(double* ary4, double* ary5, int ary4long, double* aryminus)
{
    for (int i = 0; i < ary4long; i++)
    {
        aryminus[i] = ary4[i] - ary5[i];
    }

}
//>>>>>>>>>>>>>>矩阵加法>>>>>>>>>>>>>>>>>>>>>>>>>//
void array_add(double* ary4, double* ary5, int ary4long, double* aryminus)
{
    for (int i = 0; i < ary4long; i++)
    {
        aryminus[i] = ary4[i] + ary5[i];
    }

}
//>>>>>>>>>>>>>>向量指数>>>>>>>>>>>>>>>>>>>>>>>>>//
void vector_index(double* vecindex, double* getindex, int veclong1, double index)
{
    for (int i = 0; i < veclong1; i++)
    {
        getindex[i] = pow(vecindex[i], index);
    }
}
//>>>>>>>>>>>>>>向量点积>>>>>>>>>>>>>>>>>>>>>>>>>//
double vector_dot(double* vec1, double* vec2, int veclong)
{
    double fvec = 0;
    for (int vecs = 0; vecs < veclong; vecs++)
    {
        fvec += vec1[vecs] * vec2[vecs];
    }
    return fvec;
}
//>>>>>>>>>>>>>>矩阵转置>>>>>>>>>>>>>>>>>>>>>>>>>//
void array_tpose(double* getas, double* getat, int getax, int getay)
{
    for (int i = 0; i < getax; i++)
    {
        for (int j = 0; j < getay; j++)
            getat[j * getax + i] = getas[i * getay + j];
    }

}
//>>>>>>>>>>>>>>向量相乘>>>>>>>>>>>>>>>>>>>>>>>>>//
void vector_multiply(double* vector1, double* vector2, int velong, double* vector3)
{

    for (int vi = 0; vi < velong; vi++)
    {
        vector3[vi] = vector1[vi] * vector2[vi];
    }
}
///>>>>>>>>>>>>>>利用对角矩阵向量计算对角矩阵数乘>>>>>>>>>>>>>>>>>>>>>>>>>//
void array_diagdotleft(double* yuan, double* ve, int velong, int yuanxlong, int yuanylong, double* gets)
{
    memset(gets, 0, sizeof(double) * velong * yuanylong);
    for (int m = 0; m < velong; m++)
    {
        for (int n = 0; n < yuanylong; n++)
        {
            gets[m * yuanylong + n] = ve[m] * yuan[m * yuanylong + n];
        }
    }
}
void array_diagdotright(double* yuan, double* ve, int velong, int yuanxlong, int yuanylong, double* gets)
{
    memset(gets, 0, sizeof(double) * yuanxlong * velong);
    for (int m = 0; m < yuanxlong; m++)
    {
        for (int n = 0; n < velong; n++)
        {
            gets[m * velong + n] = yuan[m * velong + n] * ve[n];
        }
    }

}
//======================================================
///>>>>>>>>>>>>> amp加速下矩阵乘法>>>>>>>>>>>>>>>>>>>>>>>>>//
void MultiplyWithAMP(double* aMatrix, double* bMatrix, double* productMatrix, int row, int inner, int col)
{

    memset(productMatrix, 0, sizeof(double) * row * col);
    int i, j, k;
    double sum = 0;
    for (i = 0; i < row; i++)
        for (k = 0; k < inner; k++)
        {
            sum = aMatrix[i * inner + k];
            for (j = 0; j < col; j++)

                productMatrix[i * col + j] += sum * bMatrix[k * col + j];
        }
}
void test(double* test, int testlong)
{

    for (int i = 0; i < testlong; i++)
    {
        cout << test[i] << "\n";
    }

}
