#include <ios>
#include<iostream>
#include<sstream>
#include <string>
#include<vector>
#include<fstream>
#include<string>
#include<iomanip>
#include<set>
using namespace std;
int m, n, slack=0, aux=0;
void show(vector<vector<double> >& A,vector<double>&b,vector<int>& c,vector<int>& Bi,vector<double>& sigma){
    int i,j;
    cout << "c:";
    for (i = 0; i < c.size(); i++)
        cout<< setw(4)<<c[i] << "  ";
    cout << endl;
    cout <<setw(4)<<"Cb"<<setw(4)<<"Xb"<<"  "<<setw(4)<<"b"<<" "<< endl;
    for (i = 0; i < A.size();i++){
        cout <<setw(4)<< c[Bi[i]] <<setw(4)<< Bi[i] <<"  "<<setw(4) << b[i];
        for (j = 0; j < A[i].size();j++)
            cout <<setprecision(4)<<setw(4)<< A[i][j]<<"  ";
        cout << endl;
    }
    cout <<"  "<<setw(12)<< "sigma"<<" ";
    for (i = 0; i < sigma.size();i++)
        cout <<setw(4)<< sigma[i]<<" ";
    cout << endl;
}

void pivot(vector<vector<double> >&A,vector<double>&b,vector<int>&c,int row,int col,vector<int>& Bi){
    cout << "row " << row << " col " << col << endl;
    int i, j;
    double key = A[row][col];
    for (i = 0; i < n + slack + aux;i++)
        A[row][i] =A[row][i]/key;
    b[row] /= key;
    for (i = 0; i < m;i++){
        if(i==row||A[i][col]==0)
            continue;
        double rate = A[i][col];
        for (j = 0; j < n + slack + aux;j++){
            A[i][j] = A[i][j] - rate * A[row][j];
        }
        b[i] = b[i] - rate * b[row];
    }
    Bi[row] = col;
    cout << "finish"
         << " pivot" << endl;
}
vector<int> initial(vector<vector<double> >&A,vector<double>& b,vector<int>&c){
    vector<int> ans;
    vector<double> sigma(n+slack+aux);
    double max_sigma;
    int max_sigma_no;
    double min_theta;
    int min_theta_no;
    int row, col;
    double epsilon = 0.000001;
    int i, j, k;
    if(aux==0){
        for (int i = 0; i < slack;i++)
            ans.push_back(n + i);
        return ans;
    }
    for (i = 0; i < aux;i++)
        ans.push_back(n + slack + i);
    while (1)
    {
        sigma.clear();
        sigma.resize(n + slack + aux);
        //calculate sigma
        for (i = 0; i < n + slack + aux; i++){
            for (j = 0; j < m;j++)
                sigma[i] -= A[j][i] * c[ans[j]];
            sigma[i] += c[i];
        }
        max_sigma = 0;
        max_sigma_no = -1;
        show(A, b, c, ans, sigma);
        //find the maximum of sigma
        for (i = 0; i < n + slack + aux;i++){
            if(sigma[i]>epsilon){
                for (j = 0; A[j][i] <= 0 && j < m;j++)
                    ;
                if(j==m){
                    ans.clear();
                    return ans;
                }
            }
            if(sigma[i]>max_sigma){
                max_sigma = sigma[i];
                max_sigma_no = i;
            }
        }
        if(max_sigma<epsilon){
            for (i = 0; i < m;i++){
                if(ans[i]>n+slack && b[i]>epsilon){
                    ans.clear();
                    return ans;
                }
            }
            return ans;
        }
        
        min_theta = -1;
        min_theta_no = 0;
        for (i = 0; i < m;i++){
            if(A[i][max_sigma_no]>epsilon && (min_theta==-1 || b[i]/A[i][max_sigma_no]<min_theta)){
                min_theta = b[i] / A[i][max_sigma_no];
                min_theta_no = i;
            }

        }
        row = min_theta_no;
        col = max_sigma_no;
        pivot(A, b, c, row, col,ans);
    }
}
void simplex(vector<vector<double> >& A,vector<double>& b,vector<int>& c){
    vector<int> Bi;
    vector<int> c1;
    double epsilon=0.00001;
    int i, j;
    // auxiliary 
    for (i = 0; i < n + slack;i++)
        c1.push_back(0);
    for (i = 0; i < aux;i++)
        c1.push_back(-1);
    Bi = initial(A, b, c1);
    if(Bi.size()==0){
        cout << "no solution" << endl;
        return;
    }

    vector<double> sigma(n+slack+aux);
    vector<double> sol(n + slack, 0);
    double max_sigma;
    int max_sigma_no;
    double min_theta;
    int min_theta_no;
    int row, col;
    while (1)
    {
        sigma.clear();
        sigma.resize(n + slack );
        //calculate sigma
        for (i = 0; i < n + slack ; i++){
            for (j = 0; j < m;j++)

                sigma[i] -= A[j][i] * c[Bi[j]];
            sigma[i] += c[i];
        }
        max_sigma = 0;
        max_sigma_no = -1;
        for (i = 0; i < n + slack ;i++){
            if(sigma[i]>0){
                for (j = 0; j<m && A[j][i] <= 0 ;j++)
                    ;
                if(j==m){
                    cout << "no bound"<<endl;
                    return;
                }
            }
            if(sigma[i]>max_sigma){
                max_sigma = sigma[i];
                max_sigma_no = i;
            }
        }
        show(A, b, c, Bi, sigma);
        if(max_sigma<epsilon){//change for precision
            set<int> s;
            for (i = 0; i < m;i++)
                s.insert(Bi[i]);
            int many = 0;
            for (i = 0; i < n + slack;i++){
                if(many)
                    break;
                if(sigma[i]<epsilon && sigma[i]>-epsilon && s.find(i)==s.end()){
                    for (j = 0; j < m;j++){
                        if(A[j][i]>0){
                            many = 1;
                            cout << "has many best solution" << endl;
                            break;
                        }
                    }
                }
            }
            for (i = 0; i < m; i++)
                sol[Bi[i]] = b[i];
            cout << "the best solution is:" << endl;
            for (i = 0; i < n+slack+aux;i++)
                cout << "x" << i+1 << "=" << sol[i] << endl;
            return;
        }
        min_theta = -1;
        min_theta_no = 0;
        for (i = 0; i < m;i++){
            if(A[i][max_sigma_no]>epsilon && (min_theta==-1 || b[i]/A[i][max_sigma_no]<min_theta)){
                min_theta = b[i] / A[i][max_sigma_no];
                min_theta_no = i;
            }

        }
        row = min_theta_no;
        col = max_sigma_no;
        pivot(A, b, c, row, col,Bi);
    }
}
int main(int argc,char* argv[]){
    fstream fin;
    fin.open(argv[1]);
    if(fin.is_open()==false){
        cout << "read in error";
        return 0;
    }
    int minmax;
    vector<int> c;
    fin >> minmax;
    string temp;
    getline(fin, temp);
    getline(fin,temp);
    stringstream sin(temp);
    int i,j,k;
    if(minmax==1){
        while(sin>>i)
            c.push_back(i);
    }else {
        minmax = 1;
        while(sin>>i)
            c.push_back(-i);
    }
    fin >> m >> n;
    vector<vector<double> > A(m);
    vector<int> sign(m);
    vector<double> b(m);
    for (i = 0; i < m;i++){
        for (j = 0; j < n;j++){
            fin >> k;
            A[i].push_back(k);
        }
        fin >> b[i];
        fin >> sign[i];
    }
    for (i = 0; i < m;i++){
        if(b[i]<0){
            for (j = 0; j < n;j++)
                A[i][j] =-A[i][j];
            b[i] *= -1;
            sign[i] = (3 - sign[i]) % 3;
        }
    }

    // add slack variable
    
    for (i = 0; i < m;i++){
        if(sign[i]==2){
            slack++;
            c.push_back(0);
            for (j = 0; j < m;j++){
                if(j==i)
                    A[j].push_back(-1);
                else
                    A[j].push_back(0);
            }
        }
        if(sign[i]==1){
            slack++;
            c.push_back(0);
            for (j = 0; j < m;j++){
                if(j==i)
                    A[j].push_back(1);
                else
                    A[j].push_back(0);
            }
        }
    }
    // add aux variable
    for (i = 0; i < m;i++)
        if(sign[i]!=1)
            break;
    if(i!=m){
        aux = m;
        for (j = 0; j < m;j++)
            for (k = 0; k < m;k++)
                if(k==j)
                    A[k].push_back(1);
                else
                    A[k].push_back(0);
    }
    simplex(A, b, c);
}