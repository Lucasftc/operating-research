#include<iostream>
#include<fstream>
#include <math.h>
#include <stdint.h>
#include<vector>
#include<sstream>
#include<set>
#include<queue>
#include<iomanip>
#include<cmath>
using namespace std;
void show(vector<vector<double> >& A,vector<double>&b,vector<double>& c,vector<int>& Bi,vector<double>& sigma){
    int i,j;
    cout << "c:";
    for (i = 0; i < c.size(); i++)
        cout<< setw(6)<<c[i] << "  ";
    cout << endl;
    cout <<setw(6)<<"Cb"<<setw(6)<<"Xb"<<"  "<<setw(6)<<"b"<<" "<< endl;
    for (i = 0; i < A.size();i++){
        cout <<setw(6)<< c[Bi[i]] <<setw(6)<< Bi[i] <<"  "<<setw(6) << b[i];
        for (j = 0; j < A[i].size();j++)
            cout <<setprecision(4)<<setw(6)<< A[i][j]<<"  ";
        cout << endl;
    }
    cout <<"  "<<setw(18)<< "sigma"<<" ";
    for (i = 0; i < sigma.size();i++)
        cout <<setw(6)<< sigma[i]<<" ";
    cout << endl;
}
void pivot(vector<vector<double> >&A,vector<double>&b,vector<double>&c,int row,int col,vector<int>& Bi,int m,int n,int slack,int aux){
    //cout << "row " << row << " col " << col << endl;
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
    //cout << "finish"
    //     << " pivot" << endl;
}
vector<int> initial(vector<vector<double> >&A,vector<double>& b,vector<double>&c,int m,int n,int slack,int aux){
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
        //show(A, b, c, ans, sigma);
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
        pivot(A, b, c, row, col,ans,m,n,slack,aux);
    }
}
vector<double> simplex(vector<vector<double> >& A,vector<double>& b,vector<double>& c,int m,int n,int slack,int aux){
    vector<int> Bi;
    vector<double> c1;
    vector<double> sol(n + slack, 0);
    double epsilon=0.00001;
    int i, j;
    // auxiliary 
    for (i = 0; i < n + slack;i++)
        c1.push_back(0);
    for (i = 0; i < aux;i++)
        c1.push_back(-1);
    Bi = initial(A, b, c1,m,n,slack,aux);
    if(Bi.size()==0){
        //cout << "no solution" << endl;
        sol.clear();
        return sol;
    }

    vector<double> sigma(n+slack+aux);
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
                    sol.clear();
                    return sol;
                }
            }
            if(sigma[i]>max_sigma){
                max_sigma = sigma[i];
                max_sigma_no = i;
            }
        }
        //show(A, b, c, Bi, sigma);
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
                            //cout << "has many best solution" << endl;
                            break;
                        }
                    }
                }
            }
            for (i = 0; i < m; i++)
                sol[Bi[i]] = b[i];
            //cout << "the best solution is:" << endl;
            //for (i = 0; i < n+slack+aux;i++)
            //    cout << "x" << i+1 << "=" << sol[i] << endl;
            return sol;
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
        pivot(A, b, c, row, col,Bi,m,n,slack,aux);
    }
}
vector<double> bound(vector<vector<double> >A,vector<double>b,vector<double>c,vector<int> sign, int m,int n){
    int i, j, k;
    int slack = 0, aux = 0;
    // add slack variable    
    for (i = 0; i < A.size();i++){
        if(sign[i]==2){
            slack++;
            c.push_back(0);
            for (j = 0; j < A.size();j++){
                if(j==i)
                    A[j].push_back(-1);
                else
                    A[j].push_back(0);
            }
        }
        if(sign[i]==1){
            slack++;
            c.push_back(0);
            for (j = 0; j < A.size();j++){
                if(j==i)
                    A[j].push_back(1);
                else
                    A[j].push_back(0);
            }
        }
    }
    // add aux variable
    for (i = 0; i < A.size();i++)
        if(sign[i]!=1)
            break;
    if(i!=A.size()){
        aux = A.size();
        for (j = 0; j < A.size();j++)
            for (k = 0; k < A.size();k++)
                if(k==j)
                    A[k].push_back(1);
                else
                    A[k].push_back(0);
    }
    return simplex(A, b, c, m, n, slack, aux);

} 
struct NODE{
    vector<vector<double> > A;
    vector<double> b;
    vector<int> sign;
    NODE(vector<vector<double> >_A,vector<double> _b,vector<int>_sign){
        A = _A;
        b=_b;
        sign = _sign;
    }
};
int main(int argc, char *argv[])
{
    fstream fin;
    fin.open("IPin.txt");
    if(fin.is_open()==false){
        fin.open(argv[1]);
        if(fin.is_open()==false){
            cout << "read in error";
            return 0;
        }
    }
    int minmax;
    vector<double> c;
    fin >> minmax;
    string temp;
    getline(fin, temp);
    getline(fin,temp);
    stringstream sin(temp);
    int m, n;
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
    queue<NODE> question;
    vector<double> ans;
    vector<double> res;
    double z=0;
    int deci_no=-1;
    double deci_value=0;
    double epsilon = 0.000001;
    double lb=-10000;
    double ub=10000;
    ans = bound(A, b,c, sign,m,n);
    if(ans.size()==0){
        cout << "no solution" << endl;
        return 0;
    }
    cout << "LP best solution" << endl;
    for (i = 0; i < c.size();i++)
        cout << ans[i] <<" ";
    cout << endl;
    z = 0;
    for (i = 0; i < c.size();i++)
        z += c[i] * ans[i];
    for (i = 0; i < c.size();i++)
        if(floor(ans[i])<ans[i]){
            deci_value = ans[i];
            deci_no = i;
            break;
        }
    if(deci_no==-1){
        cout << "res is:";
        for (i = 0; i < ans.size();i++)
            cout << ans[i]<<" ";
        return 0;
    }else{
        ub=z;
        NODE new1(A,b,sign);
        NODE new2(A,b,sign);
        vector<double> line;
        for (i = 0; i < A[0].size();i++){
            if(i!=deci_no)
                line.push_back(0);
            else
                line.push_back(1);
        }
        new1.A.push_back(line);
        new2.A.push_back(line);
        new1.b.push_back(floor(deci_value));
        new1.sign.push_back(1);
        new2.b.push_back(ceil(deci_value));
        new2.sign.push_back(2);
        question.push(new1);
        question.push(new2);
    }
    int time = 0;
    while(question.empty()==false && time<32){
        time++;
        NODE now = question.front();
        question.pop();
        deci_value = 0;
        deci_no = -1;
        ans = bound(now.A, now.b, c, now.sign, now.A.size(), now.A[0].size());
        cout << "node "<<time << " ";
        if(ans.size()==0){
            cout << "no solution" << endl;
            continue;
        }
        cout << "LP solution is";
        for (i = 0; i < c.size();i++)
            cout << ans[i] << " ";
        z = 0;
        for (i = 0; i < c.size();i++)
            z += c[i] * ans[i];
        if(ans.size()==0){
            continue;
        }
        for (i = 0; i < c.size();i++)
            if(floor(ans[i]+epsilon)+epsilon<ans[i]){
                deci_value = ans[i];
                deci_no = i;
                break;
            }
        if(deci_no==-1){
            cout << " integer " << endl;
            if(z>=lb){
                lb = z;
                res = ans;
            
            if(ub<lb+1){
                cout << "res is:";
                for (i = 0; i < res.size();i++)
                    cout << res[i]<<" ";
                cout << endl;
                //return 0;
            }
            }
        }else if (z>lb){
            cout << endl;
            //cout << " decino" << deci_no << "decivalue" << deci_value << endl;
            NODE new1(now.A,now.b,now.sign);
            NODE new2(now.A,now.b,now.sign);
            vector<double> line;
            for (i = 0; i < now.A[0].size();i++){
                if(i!=deci_no)
                    line.push_back(0);
                else
                    line.push_back(1);
            }
            new1.A.push_back(line);
            new2.A.push_back(line);
            new1.b.push_back(floor(deci_value));
            new1.sign.push_back(1);
            new2.b.push_back(ceil(deci_value));
            new2.sign.push_back(2);
            question.push(new1);
            question.push(new2);
        }else{
            cout << "abandon" << endl;
        }
    }
    cout << "one posible res is:" << endl;
    for (i = 0; i < res.size();i++)
        cout << res[i] << " ";

}
