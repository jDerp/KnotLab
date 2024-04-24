#include <iRRAM.h>
#include "iRRAMx/polynomial.hpp"
using namespace iRRAM;
using std::vector;
struct Lpolynomial{
    int mindeg, maxdeg;
    vector<int> coeflist; // maxdeg -> mindeg

    // poly(e^(2*pi*(m/n)*i))
    COMPLEX subs(int m, int n){
        int deg = maxdeg;
        COMPLEX ans = COMPLEX(0, 0);
        
        for(auto x: coeflist){
            // cout << x << " ";
            ans = ans + COMPLEX(x)*exp(COMPLEX(0, REAL(2*deg*m) * pi() / REAL(n)));
            deg--;
        }
        // cout << "\n";
        assert(deg == mindeg - 1);
        return ans;
    }
    
            
};

Lpolynomial operator + (const Lpolynomial &m, const Lpolynomial &n){
    Lpolynomial sum = {min(m.mindeg, n.mindeg), max(m.maxdeg, n.maxdeg), vector<int>(max(m.maxdeg, n.maxdeg) - min(m.mindeg, n.mindeg) + 1, 0)};
    for(int i = 0; i <= m.maxdeg - m.mindeg; i++){
        sum.coeflist[sum.maxdeg - m.maxdeg + i] += m.coeflist[i];
    }
    for(int j = 0; j <= n.maxdeg - n.mindeg; j++){
        sum.coeflist[sum.maxdeg - n.maxdeg + j] += n.coeflist[j];
    }
    while(sum.mindeg < 0 && sum.coeflist.back() == 0){
        sum.mindeg++;
        sum.coeflist.pop_back();
    }
    while(sum.maxdeg > 0 && sum.coeflist.front() == 0){
        sum.maxdeg--;
        sum.coeflist.erase(sum.coeflist.begin());
    }
    
    return sum;
}

Lpolynomial operator * (const Lpolynomial &m, const Lpolynomial &n){
    Lpolynomial product = {m.mindeg + n.mindeg, m.maxdeg + n.maxdeg, vector<int>(m.maxdeg + n.maxdeg - m.mindeg - n.mindeg + 1, 0)};
    for(int i = 0; i <= m.maxdeg - m.mindeg; i++){
        for(int j = 0; j <= n.maxdeg - n.mindeg; j++){
            product.coeflist[i + j] += m.coeflist[i]*n.coeflist[j];
        }
    }
    while(product.mindeg < 0 && product.coeflist.back() == 0){
        product.mindeg++;
        product.coeflist.pop_back();
    }
    while(product.maxdeg > 0 && product.coeflist.front() == 0){
        product.maxdeg--;
        product.coeflist.erase(product.coeflist.begin());
    }
    
    return product;
}
    
void compute(){
    int MATDIM;
    
    cout << "Input n : ";
    cin >> MATDIM;
    MATDIM--;

    vector<vector<Lpolynomial>> POLYNOMIALMAT = vector<vector<Lpolynomial>>(MATDIM, vector<Lpolynomial>(MATDIM, Lpolynomial({0, 0, {0}})));
    // vector<vector<COMPLEX>> COMPLEXMAT;

    for(int i = 0; i < MATDIM; i++){
        POLYNOMIALMAT[i][i] = Lpolynomial({0, 0, {1}});
    }

    //Burau Representation of Braids

    int braid_num;
    cout << "Input number of braids : ";
    cin >> braid_num;
    assert(braid_num > 0);

    cout << "Input braids : ";
    while(braid_num--){
        int braid;
        cin >> braid;

        if(braid > MATDIM || braid < -MATDIM || braid == 0){
            cout << "ERROR: Braid out of bound, please input again.\n";
            return ;
        }

        vector<vector<Lpolynomial>> RESULTMAT = vector<vector<Lpolynomial>>(MATDIM, vector<Lpolynomial>(MATDIM, Lpolynomial({0, 0, {0}})));
        vector<vector<Lpolynomial>> BRAIDMAT = vector<vector<Lpolynomial>>(MATDIM, vector<Lpolynomial>(MATDIM, Lpolynomial({0, 0, {0}})));
        
        for(int i = 0; i < MATDIM; i++){
            BRAIDMAT[i][i] = Lpolynomial({0, 0, {1}});
        }

        if(braid > 0){
            BRAIDMAT[braid-1][braid-1] = Lpolynomial({0, 1, {-1, 0}}); // -t
            if(braid != 1)  BRAIDMAT[braid-1][braid-2] = Lpolynomial({0, 1, {1, 0}}); // t
            if(braid != MATDIM) BRAIDMAT[braid-1][braid] = Lpolynomial({0, 0, {1}}); // 1
        }

        else if(braid < 0){
            braid = -braid;
            BRAIDMAT[braid-1][braid-1] = Lpolynomial({-1, 0, {0, -1}}); // -t^-1
            if(braid != 1)  BRAIDMAT[braid-1][braid-2] = Lpolynomial({0, 0, {1}}); // 1
            if(braid != MATDIM) BRAIDMAT[braid-1][braid] = Lpolynomial({-1, 0, {0, 1}}); // t^-1
        }

        for(int i = 0; i < MATDIM; i++){
            for(int j = 0; j < MATDIM; j++){
                for(int k = 0; k < MATDIM; k++){
                    RESULTMAT[i][j] = RESULTMAT[i][j] + (POLYNOMIALMAT[i][k] * BRAIDMAT[k][j]);
                }
            }
        }

        POLYNOMIALMAT = RESULTMAT;

    }

    // for(auto x: POLYNOMIALMAT){
    //     for(auto y: x){
    //         cout << y.mindeg << ", " << y.maxdeg << ", ";
    //         for(auto z: y.coeflist){
    //             cout << z << " ";
    //         }
    //         cout << "| ";
    //     }
    //     cout << "\n";
    // }

    // Faddeev-Lavierrier Characteristic Polynomial

    vector<Lpolynomial> poly_charpoly;
    vector<vector<Lpolynomial>> CMAT = POLYNOMIALMAT;
    poly_charpoly.push_back({0, 0, {1}});
    for(int i = 0; i < MATDIM; i++){
        if(i > 0){
            // C = A*(C + c_prev*I)
            vector<vector<Lpolynomial>> RESULTMAT = vector<vector<Lpolynomial>>(MATDIM, vector<Lpolynomial>(MATDIM, Lpolynomial({0, 0, {0}})));

            for(int i = 0; i < MATDIM; i++){
                for(int j = 0; j < MATDIM; j++){
                    for(int k = 0; k < MATDIM; k++){
                        if(k == j)  RESULTMAT[i][j] = RESULTMAT[i][j] + (POLYNOMIALMAT[i][k] * (poly_charpoly.back() + CMAT[k][j]));
                        else        RESULTMAT[i][j] = RESULTMAT[i][j] + (POLYNOMIALMAT[i][k] * CMAT[k][j]);
                    }
                }
            }
            CMAT = RESULTMAT;
        }

        Lpolynomial charpolycoef = {0, 0, {0}};

        //c_i = -tr(C)/i
        for(int j = 0; j < MATDIM; j++){
            charpolycoef = charpolycoef + CMAT[j][j];
        }
        for(auto &x: charpolycoef.coeflist){
            x /= -(i+1);
        }
        poly_charpoly.push_back(charpolycoef);
    }

    // for(auto x: poly_charpoly){
    //     cout << x.mindeg << " " << x.maxdeg << " ";
    //     for(auto y: x.coeflist){
    //         cout << y << " ";
    //     }
    //     cout << "\n";
    // }

    // Complex Substitution
    
    vector<COMPLEX> complex_charpoly;
    int m, n;
    cout << "Input m n for e^(2*pi*i*m/n): ";
    cin >> m >> n;
    for(auto x: poly_charpoly){
        complex_charpoly.push_back(x.subs(m,n));
    }

    // Root Calculation 

    vector<COMPLEX> root;
    POLYNOMIAL charpoly = POLYNOMIAL(MATDIM, complex_charpoly);
    root = roots(charpoly);
    for(auto x: root){
        cout << real(x) << " + " << imag(x) << " i with abs = " << abs(x) << "\n";
    }
    
}