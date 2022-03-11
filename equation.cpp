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
    ��˹��ȥ��
    �㷨���ȶ�
*/
double* guass(double** a,double* b,int n){
    //��Ԫ
    for(int i=0;i<n-1;i++){
        bool flag=0;
        if(a[i][i]==0){//��ԪΪ�㣬����
            for(int j=i+1;j<n;j++){
                if(a[j][i]!=0){
                    flag=1;
                    for(int k=i;k<n;k++){//����
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

    for(int i=0;i<n;i++){//��������Ƿ�����
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
    ��˹����Ԫ��ȥ��
    �㷨�ȶ�
*/
double* guassCol(double** a,double* b,int n){
    //��Ԫ
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
            for(int k=i;k<n;k++){//����
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

    for(int i=0;i<n;i++){//��������Ƿ�����
        for(int j=0;j<n;j++){
            cout<<setw(5)<<a[i][j];
        }
        cout<<setw(5)<<b[i];
        cout<<endl;
    }

    for(int i=n-1;i>=0;i--){//�ش�
        for(int j=n-1;j>i;j--){
            b[i]=b[i]-a[i][j]*b[j];
        }
        b[i]=b[i]/a[i][i];
    }

    return b;
}

/*
    LU�ֽ�
    ʵ���Ͼ�����ͨ�ĸ�˹��ȥ��
    �����ڽ���ϵ��������ͬ�ķ���
    �㷨���ȶ�
*/

//�������� 4 -2 0 4 -2 2 -3 1 0 -3 13 -7 4 1 -7 23
double** LU_factorization(double** a,int n){
    //����ĵ�һ�в���
    for(int i=0;i<n;i++){
        //��
        for(int j=i;j<n;j++){
            for(int k=0;k<i;k++){
                a[i][j]=a[i][j]-a[i][k]*a[k][j];
            }
        }

        //��
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
    //����ĵ�һ�в���
    for(int i=0;i<n;i++){
        //��
        for(int j=i;j<n;j++){
            for(int k=0;k<i;k++){
                a[i][j]=a[i][j]-a[i][k]*a[k][j];
            }
        }

        for(int k=0;k<i;k++){
            b[i]=b[i]-b[k]*a[i][k];
        }

        //��
        for(int j=i+1;j<n;j++){
            for(int k=0;k<i;k++){
                a[j][i]=a[j][i]-a[j][k]*a[k][i];
            }
            a[j][i]=a[j][i]/a[i][i];
        }
    }
    return a;
}


//�������� 9 18 9 -27 18 45 0 -45 9 0 126 9 -27 -45 9 135
//1 2 16 8
double* guassLU(double** a,double* b,int n){
    LUY(a,b,n);
    for(int i=n-1;i>=0;i--){//�ش�
        for(int j=n-1;j>i;j--){
            b[i]=b[i]-a[i][j]*b[j];
        }
        b[i]=b[i]/a[i][i];
    }
    return b;
}

/*
    ƽ������
    �����ڶԳ���������
    �㷨�ȶ�
*/

//LDL�ֽ⣨�Գƾ���
//�õ�LU�ֽ��͵õ���LDL�ֽ�

//GG�ֽ� ����˹���ֽ� Cholesky�ֽ�
//�õ�LDL�ֽ�͵õ���Cholesky�ֽ�
//������˹���ֽ�����һ�������򻯵Ĺ�ʽ���Լ��ټ�������������LU�ֽ�
//�������� 9 18 9 -27 18 45 0 -45 9 0 126 9 -27 -45 9 135
//1 2 16 8
double** Cholesky(double** a,int n){
    for(int i=0;i<n;i++){
        //�Խ�
        for(int k=0;k<i;k++){
            a[i][i]=a[i][i]-a[i][k]*a[i][k];
        }
        a[i][i]=sqrt(a[i][i]);

        //��
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
        //�Խ�
        for(int k=0;k<i;k++){
            a[i][i]=a[i][i]-a[i][k]*a[i][k];
            b[i]=b[i]-b[k]*a[i][k];
        }
        a[i][i]=sqrt(a[i][i]);
        b[i]=b[i]/a[i][i];

        //��
        for(int j=i+1;j<n;j++){
            for(int k=0;k<i;k++){
                a[j][i]=a[j][i]-a[j][k]*a[i][k];
            }
            a[j][i]=a[j][i]/a[i][i];
        }

        //��ȫ����
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
    for(int i=n-1;i>=0;i--){//�ش�
        for(int j=n-1;j>i;j--){
            b[i]=b[i]-a[i][j]*b[j];
        }
        b[i]=b[i]/a[i][i];
    }
    return b;
}

/*
    �Ľ�ƽ������
    ����LDL�ֽ�
    �����ƽ����������n����������
    �ɶԳ��ԣ�����LDL�ֽ�ʱ����ֻ��LU�ֽ��е�U����
*/

//�������� 9 18 9 -27 18 45 0 -45 9 0 126 9 -27 -45 9 135
//1 2 16 8
double** LDL(double** a,int n){
    for(int i=0;i<n;i++){
        //�Խ���
        for(int k=0;k<i;k++){
            a[i][i]=a[i][i]-a[k][k]*a[i][k]*a[i][k];
        }

        //��
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
        //�Խ���
        for(int k=0;k<i;k++){
            a[i][i]=a[i][i]-a[k][k]*a[i][k]*a[i][k];
            b[i]=b[i]-b[k]*a[i][k];//˳����y����
        }

        //��
        for(int j=i+1;j<n;j++){
            for(int k=0;k<i;k++){
                a[j][i]=a[j][i]-a[j][k]*a[i][k]*a[k][k];
            }
            a[j][i]=a[j][i]/a[i][i];
        }
    }
    return b;
}

//�������� 9 18 9 -27 18 45 0 -45 9 0 126 9 -27 -45 9 135
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
    ������ԽǷ������׷�Ϸ�
    �ϸ�Խ�ռ�ž���
    �������õ���LU�ֽ�
    �㷨�ȶ�
*/

//�������� 0 2 -3 4 -5 1 3 4 7 6 2 1 2 1 0
//1 3 4 7 6
//2 1 2 1 0
double* Chase3(double** a,double* b,int n){//nΪ�Խ��߳���
    for(int i=1;i<n;i++){
        a[0][i]=a[0][i]/a[1][i-1];
        a[1][i]=a[1][i]-a[0][i]*a[2][i-1];
        b[i]=b[i]-b[i-1]*a[0][i];//˳����y����
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

    cout<<"�����뷽��Ľ�����";
    cin>>n;

    /*double** a=new double*[n];
    for(int i=0;i<n;i++){
        a[i]=new double[n];
    }*/

    double** a=new double*[3];//��׷���þ���
    for(int i=0;i<3;i++){
        a[i]=new double[n];
    }

    double* b=new double[n];

    cout<<"������ϵ�������Ԫ�أ�"<<endl;//��׷��������
    for(int i=0;i<3;i++){
        for(int j=0;j<n;j++){
            cin>>a[i][j];
        }
    }

    /*cout<<"������ϵ�������Ԫ�أ�"<<endl;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cin>>a[i][j];
        }
    }*/
    cout<<"������������Ԫ�أ�"<<endl;
    for(int i=0;i<n;i++){
        cin>>b[i];
    }

    /*try{
        result=guassCol(a,b,n);
    }
    catch(int e){
        cout<<"���������˳������ʽ����"<<endl;
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
