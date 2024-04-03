#include <iRRAM.h>
#include "iRRAMx/polynomial.hpp"
using namespace iRRAM;
using std::vector;
struct Lpolynomial
    {
        int mindeg, maxdeg;
        vector<int> coeflist;
        // poly(e^(2*pi*(m/n)*i))
        COMPLEX subs(int m, int n){
            int deg = mindeg;
            COMPLEX ans = COMPLEX(0, 0);
            
            for(auto x: coeflist){
                // cout << x << " ";
                ans = ans + COMPLEX(x)*exp(COMPLEX(0, REAL(2*deg*m) * pi() / REAL(n)));
                deg++;
            }
            // cout << "\n";
            assert(deg == maxdeg + 1);
            return ans;
        }
                
    };
POLYNOMIAL det(vector<vector<POLYNOMIAL>> &MAT, vector<int> &perm, int row){
    int n = perm.size();
    if(n == 2){
        return MAT[row][perm[0]] * MAT[row+1][perm[1]] - MAT[row+1][perm[0]] * MAT[row][perm[1]];
    }
    POLYNOMIAL ans = POLYNOMIAL();
    for(int i = 0; i < n; i++){
        int pos = perm[i];
        perm.erase(perm.begin()+i);
        POLYNOMIAL minor;
        POLYNOMIAL minor = det(MAT, perm, row+1);
        perm.insert(perm.begin()+i, pos);
        if((row+i) % 2) ans = ans + MAT[row][perm[i]]*minor;
        else            ans = ans - MAT[row][perm[i]]*minor;
    }
    
}
    
void compute(){
    int MATDIM;
    vector<vector<Lpolynomial>> POLYNOMIALMAT;
    vector<vector<COMPLEX>> COMPLEXMAT;
    cout << "Input matrix dimension : ";
    cin >> MATDIM;
    for(int i = 0; i < MATDIM; i++){
        POLYNOMIALMAT.push_back({});
        for(int j = 0; j < MATDIM; j++){
            int mindeg, maxdeg;
            int coef;
            vector<int> poly;
            // a_0 + a_1x + ... + a_dx^d
            cin >> mindeg >> maxdeg;
            for(int k = mindeg; k <= maxdeg; k++){
                cin >> coef;
                poly.push_back(coef);
            }
            POLYNOMIALMAT[i].push_back({mindeg, maxdeg, poly});
        }
    }
    int m, n;
    cin >> m >> n;
    for(int i = 0; i < MATDIM; i++){
        COMPLEXMAT.push_back({});
        for(int j = 0; j < MATDIM; j++){
            COMPLEXMAT[i].push_back(POLYNOMIALMAT[i][j].subs(m, n));
            // cout << real(COMPLEXMAT[i][j]) << " " << imag(COMPLEXMAT[i][j]) << " ";
        }
        // cout << "\n";
    }
    // cout << "ggggg\n";
    vector<vector<POLYNOMIAL>> EIGENMAT;
    for(int i = 0; i < MATDIM; i++){
        EIGENMAT.push_back({});
        for(int j = 0; j < MATDIM; j++){
            vector<COMPLEX> temp;
            temp.push_back(COMPLEXMAT[i][j]);
            if(i == j){
                temp.push_back(COMPLEX(-1));
                EIGENMAT[i].push_back(POLYNOMIAL(1, temp));
            }
            else{
                EIGENMAT[i].push_back(POLYNOMIAL(0, temp));
            }
        }
        // cout << "\n";

    }
    vector<COMPLEX> root;
    if(MATDIM == 2){
        const POLYNOMIAL charpoly = (EIGENMAT[0][0])*(EIGENMAT[1][1]) - (EIGENMAT[0][1])*(EIGENMAT[1][0]);
        // const POLYNOMIAL charpoly = POLYNOMIAL(2, {COMPLEX(-30, 7), COMPLEX(-7, -1), COMPLEX(1)});
        cout << charpoly << "\n";
        root = roots(charpoly);
    }
    if(MATDIM == 3){
        const POLYNOMIAL charpoly = EIGENMAT[0][0]*EIGENMAT[1][1]*EIGENMAT[2][2]
                + EIGENMAT[0][1]*EIGENMAT[1][2]*EIGENMAT[2][0]
                + EIGENMAT[0][2]*EIGENMAT[1][0]*EIGENMAT[2][1]
                - EIGENMAT[0][0]*EIGENMAT[1][2]*EIGENMAT[2][1]
                - EIGENMAT[0][1]*EIGENMAT[1][0]*EIGENMAT[2][2]
                - EIGENMAT[0][2]*EIGENMAT[1][1]*EIGENMAT[2][0];
        root = roots(charpoly);
    }

    else{
        vector<int> perm;
        for(int i = 0; i < MATDIM; i++){
            perm.push_back(i);
        }
        const POLYNOMIAL charpoly = det(EIGENMAT, perm, 0);
        root = roots(charpoly);
    }
    for(auto x: root){
        cout << real(x) << " " << imag(x) << "\n";
    }
    
}