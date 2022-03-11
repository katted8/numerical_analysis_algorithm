#include<iostream>
#include<iomanip>
#include<cmath>
using namespace std;

double aValue(double x){
    return x>=0?x:-x;
}

void printVec(double* Vec,int n){
    for(int i=0;i<n;i++){
        cout<<setw(10)<<Vec[i]<<" ";
    }
    cout<<endl;
}

void printMat(double** Mat,int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout<<setw(10)<<Mat[i][j];
        }
        cout<<endl;
    }
}

void printMat_Vec(double** Mat,double* Vec,int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout<<setw(10)<<Mat[i][j];
        }
        cout<<setw(10)<<Vec[i];
        cout<<endl;
    }
}

/*
    高斯消去法
    算法不稳定
*/
double* guass(double** a,double* b,int n){
    //消元
    for(int i=0;i<n-1;i++){
        bool flag=0;
        if(a[i][i]==0){//主元为零，换行
            for(int j=i+1;j<n;j++){
                if(a[j][i]!=0){
                    flag=1;
                    for(int k=i;k<n;k++){//换行
                        double temp=a[i][k];
                        a[i][k]=a[j][k];
                        a[j][k]=temp;
                    }
                    double temp=b[i];
                    b[i]=b[j];
                    b[j]=temp;
                    break;
                }
            }
            if(flag==0){
                throw 0;
            }
        }

        for(int j=i+1;j<n;j++){
            double l=a[j][i]/a[i][i];
            for(int k=i;k<n;k++){
                a[j][k]=a[j][k]-l*a[i][k];
            }
            b[j]=b[j]-l*b[i];
        }

    }

    for(int i=0;i<n;i++){//输出上三角方程组
        for(int j=0;j<n;j++){
            cout<<setw(5)<<a[i][j];
        }
        cout<<setw(5)<<b[i]<<endl;
    }

    for(int i=n-1;i>=0;i--){
        for(int j=n-1;j>i;j--){
            b[i]=b[i]-a[i][j]*b[j];
        }
        b[i]=b[i]/a[i][i];
    }

    return b;
}

/*
    高斯列主元消去法
    算法稳定
*/
double* guassCol(double** a,double* b,int n){
    //消元
    for(int i=0;i<n-1;i++){
        double maxNum=aValue(a[i][i]);
        int maxRow=i;
        for(int j=i+1;j<n;j++){
            if(aValue(a[j][i])>maxNum){
                maxNum=aValue(a[j][i]);
                maxRow=j;
            }
        }
        if(maxNum==0){
            throw 0;
        }
        else{
            for(int k=i;k<n;k++){//换行
                double temp=a[i][k];
                a[i][k]=a[maxRow][k];
                a[maxRow][k]=temp;
            }
            double temp=b[i];
            b[i]=b[maxRow];
            b[maxRow]=temp;
        }

        for(int j=i+1;j<n;j++){
            double l=a[j][i]/a[i][i];
            for(int k=i;k<n;k++){
                a[j][k]=a[j][k]-l*a[i][k];
            }
            b[j]=b[j]-l*b[i];
        }

    }

    for(int i=0;i<n;i++){//输出上三角方程组
        for(int j=0;j<n;j++){
            cout<<setw(5)<<a[i][j];
        }
        cout<<setw(5)<<b[i];
        cout<<endl;
    }

    for(int i=n-1;i>=0;i--){//回代
        for(int j=n-1;j>i;j--){
            b[i]=b[i]-a[i][j]*b[j];
        }
        b[i]=b[i]/a[i][i];
    }

    return b;
}

/*
    LU分解
    实际上就是普通的高斯消去法
    适用于解多个系数矩阵相同的方程
    算法不稳定
*/

//测试用例 4 -2 0 4 -2 2 -3 1 0 -3 13 -7 4 1 -7 23
double** LU_factorization(double** a,int n){
    //矩阵的第一行不变
    for(int i=0;i<n;i++){
        //行
        for(int j=i;j<n;j++){
            for(int k=0;k<i;k++){
                a[i][j]=a[i][j]-a[i][k]*a[k][j];
            }
        }

        //列
        for(int j=i+1;j<n;j++){
            for(int k=0;k<i;k++){
                a[j][i]=a[j][i]-a[j][k]*a[k][i];
            }
            a[j][i]=a[j][i]/a[i][i];
        }
    }
    return a;
}

double** LUY(double** a,double* b,int n){
    //矩阵的第一行不变
    for(int i=0;i<n;i++){
        //行
        for(int j=i;j<n;j++){
            for(int k=0;k<i;k++){
                a[i][j]=a[i][j]-a[i][k]*a[k][j];
            }
        }

        for(int k=0;k<i;k++){
            b[i]=b[i]-b[k]*a[i][k];
        }

        //列
        for(int j=i+1;j<n;j++){
            for(int k=0;k<i;k++){
                a[j][i]=a[j][i]-a[j][k]*a[k][i];
            }
            a[j][i]=a[j][i]/a[i][i];
        }
    }
    return a;
}


//测试用例 9 18 9 -27 18 45 0 -45 9 0 126 9 -27 -45 9 135
//1 2 16 8
double* guassLU(double** a,double* b,int n){
    LUY(a,b,n);
    for(int i=n-1;i>=0;i--){//回代
        for(int j=n-1;j>i;j--){
            b[i]=b[i]-a[i][j]*b[j];
        }
        b[i]=b[i]/a[i][i];
    }
    return b;
}

/*
    平方根法
    适用于对称正定矩阵
    算法稳定
*/

//LDL分解（对称矩阵）
//得到LU分解后就得到了LDL分解

//GG分解 楚列斯基分解 Cholesky分解
//得到LDL分解就得到了Cholesky分解
//但楚列斯基分解自有一套稍作简化的公式可以减少计算量，不用作LU分解
//测试用例 9 18 9 -27 18 45 0 -45 9 0 126 9 -27 -45 9 135
//1 2 16 8
double** Cholesky(double** a,int n){
    for(int i=0;i<n;i++){
        //对角
        for(int k=0;k<i;k++){
            a[i][i]=a[i][i]-a[i][k]*a[i][k];
        }
        a[i][i]=sqrt(a[i][i]);

        //列
        for(int j=i+1;j<n;j++){
            for(int k=0;k<i;k++){
                a[j][i]=a[j][i]-a[j][k]*a[i][k];
            }
            a[j][i]=a[j][i]/a[i][i];
        }
    }
    return a;
}

double** Cholesky_Y(double** a,double* b,int n){
    for(int i=0;i<n;i++){
        //对角
        for(int k=0;k<i;k++){
            a[i][i]=a[i][i]-a[i][k]*a[i][k];
            b[i]=b[i]-b[k]*a[i][k];
        }
        a[i][i]=sqrt(a[i][i]);
        b[i]=b[i]/a[i][i];

        //列
        for(int j=i+1;j<n;j++){
            for(int k=0;k<i;k++){
                a[j][i]=a[j][i]-a[j][k]*a[i][k];
            }
            a[j][i]=a[j][i]/a[i][i];
        }

        //补全矩阵
        for(int i=0;i<n-1;i++){
            for(int j=i+1;j<n;j++){
                a[i][j]=a[j][i];
            }
        }
    }
    return a;
}

double* sqrtGG(double** a,double* b,int n){
    Cholesky_Y(a,b,n);
    for(int i=n-1;i>=0;i--){//回代
        for(int j=n-1;j>i;j--){
            b[i]=b[i]-a[i][j]*b[j];
        }
        b[i]=b[i]/a[i][i];
    }
    return b;
}

/*
    改进平方根法
    基于LDL分解
    相比于平方根法少做n个开方运算
    由对称性，在算LDL分解时，可只算LU分解中的U部分
*/

//测试用例 9 18 9 -27 18 45 0 -45 9 0 126 9 -27 -45 9 135
//1 2 16 8
double** LDL(double** a,int n){
    for(int i=0;i<n;i++){
        //对角线
        for(int k=0;k<i;k++){
            a[i][i]=a[i][i]-a[k][k]*a[i][k]*a[i][k];
        }

        //列
        for(int j=i+1;j<n;j++){
            for(int k=0;k<i;k++){
                a[j][i]=a[j][i]-a[j][k]*a[i][k]*a[k][k];
            }
            a[j][i]=a[j][i]/a[i][i];
        }
    }
    return a;
}

double* LDL_Y(double** a,double* b,int n){
    for(int i=0;i<n;i++){
        //对角线
        for(int k=0;k<i;k++){
            a[i][i]=a[i][i]-a[k][k]*a[i][k]*a[i][k];
            b[i]=b[i]-b[k]*a[i][k];//顺带把y算了
        }

        //列
        for(int j=i+1;j<n;j++){
            for(int k=0;k<i;k++){
                a[j][i]=a[j][i]-a[j][k]*a[i][k]*a[k][k];
            }
            a[j][i]=a[j][i]/a[i][i];
        }
    }
    return b;
}

//测试用例 9 18 9 -27 18 45 0 -45 9 0 126 9 -27 -45 9 135
//1 2 16 8
double* improvedSqrtLDL(double** a,double* b,int n){
    LDL_Y(a,b,n);
    for(int i=n-1;i>=0;i--){
        b[i]=b[i]/a[i][i];
        for(int k=i+1;k<n;k++){
            b[i]=b[i]-b[k]*a[k][i];
        }
    }
    return b;
}

/*
    求解三对角方程组的追赶法
    严格对角占优矩阵
    本质上用的是LU分解
    算法稳定
*/

//测试用例 0 2 -3 4 -5 1 3 4 7 6 2 1 2 1 0
//1 3 4 7 6
//2 1 2 1 0
double* Chase3(double** a,double* b,int n){//n为对角线长度
    for(int i=1;i<n;i++){
        a[0][i]=a[0][i]/a[1][i-1];
        a[1][i]=a[1][i]-a[0][i]*a[2][i-1];
        b[i]=b[i]-b[i-1]*a[0][i];//顺带把y算了
    }

    b[n-1]=b[n-1]/a[1][n-1];
    for(int i=n-2;i>=0;i--){
        b[i]=(b[i]-a[2][i]*b[i+1])/a[1][i];
    }

    return b;
}









//main/////////////////////////////////////////////////////////////////////////////////////////


int main(){
    int n;
    double* result;

    cout<<"请输入方阵的阶数：";
    cin>>n;

    /*double** a=new double*[n];
    for(int i=0;i<n;i++){
        a[i]=new double[n];
    }*/

    double** a=new double*[3];//三追逐用矩阵
    for(int i=0;i<3;i++){
        a[i]=new double[n];
    }

    double* b=new double[n];

    cout<<"请输入系数矩阵的元素："<<endl;//三追逐用输入
    for(int i=0;i<3;i++){
        for(int j=0;j<n;j++){
            cin>>a[i][j];
        }
    }

    /*cout<<"请输入系数矩阵的元素："<<endl;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cin>>a[i][j];
        }
    }*/
    cout<<"请输入向量的元素："<<endl;
    for(int i=0;i<n;i++){
        cin>>b[i];
    }

    /*try{
        result=guassCol(a,b,n);
    }
    catch(int e){
        cout<<"不满足各阶顺序主子式非零"<<endl;
    }

    for(int i=0;i<n;i++){
        cout<<result[i]<<" ";
    }*/

    /*result=LU_factorization(a,n);
    printMat(result,n);*/

    /*result=guassLU(a,b,n);
    printVec(result,n);*/

    /*Cholesky_Y(a,b,n);
    printMat_Vec(a,b,n);*/

    /*sqrtGG(a,b,n);
    printVec(b,n);*/

    /*improvedSqrtLDL(a,b,n);
    printVec(b,n);*/

    Chase3(a,b,n);
    printVec(b,n);

    return 0;
}
