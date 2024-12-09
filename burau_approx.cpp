#include <iRRAM.h>
#include "iRRAMx/polynomial.hpp"
using namespace iRRAM;
using std::vector;

struct Lpolynomial{
    int mindeg, maxdeg;
    vector<int> coeflist; // maxdeg -> mindeg

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
    int MATDIM; // n-1
    
    cout << "Input n:\n";
    cin >> MATDIM;
    MATDIM--;

    vector<vector<Lpolynomial>> POLYNOMIALMAT = vector<vector<Lpolynomial>>(MATDIM, vector<Lpolynomial>(MATDIM, Lpolynomial({0, 0, {0}})));

    for(int i = 0; i < MATDIM; i++){
        POLYNOMIALMAT[i][i] = Lpolynomial({0, 0, {1}});
    }

    //Burau Representation of Braids

    int braid_num;
    cout << "Input number of braid generators:\n";
    cin >> braid_num;
    assert(braid_num > 0);

    cout << "Input braid generators:\n";
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

    // Faddeev-LaVerrier Characteristic Polynomial

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

    int prec = 20;

    int mode = 0;
    
    while(1){
        cout << "(1) calculate EVs\n(2) test reducibility condition\n(3) test exchangeability condition\n(4) set precision 2^-prec (current prec = " << prec << ")\n(5) exit\n";
        cin >> mode;
        if(mode == 4){
            cout << "Input prec:\n";
            cin >> prec;
        }
        if(mode == 5){
            return ;
        }
    
        if(mode == 1){
            int p, q;
            cout << "Input \"p q\" for e^(2*pi*i*p/q):\n";
            cin >> p >> q;
    
            // Complex Substitution
    
            vector<COMPLEX> complex_charpoly;
    
            for(auto x: poly_charpoly){
                complex_charpoly.push_back(x.subs(p,q));
            }
    
            // Root Calculation 
    
            vector<COMPLEX> root;
            POLYNOMIAL charpoly = POLYNOMIAL(MATDIM, complex_charpoly);
    
            root = roots(charpoly);
    
            for(auto x: root){
                cout << real(x) << " + " << imag(x) << " i with abs = " << abs(x) << "\n";
            }
        }
        if(mode == 2 || mode == 3){
            bool irr = 0, exch = 0;
            if(mode == 2) exch = 1;
            if(mode == 3) irr = 1;
            int MINQ, MAXQ;
            cout << "Input \"minq maxq\" for the range of q in e^(2*pi*i*p/q) to be tested:\n";
            cin >> MINQ >> MAXQ;
            MINQ = max(MINQ, MATDIM + 2);
            for(int q = MINQ; q <= MAXQ; q++){
                for(int p = q/(MATDIM+1) + 1; MATDIM*p < q; p++){
    
                    if(std::__gcd(p, q) != 1)    continue;
    
                    // Complex Substitution
    
                    vector<COMPLEX> complex_charpoly;
    
                    for(auto x: poly_charpoly){
                        complex_charpoly.push_back(x.subs(p,q));
                    }
    
                    cout << "testing " << p << " / " << q << "...\n";
    
                    // Root Calculation 
    
                    vector<COMPLEX> root;
                    POLYNOMIAL charpoly = POLYNOMIAL(MATDIM, complex_charpoly);
    
                    root = roots(charpoly);
    
                    // Root Testing
    
                    // Test Unit Circle
                    COMPLEX non_unit; bool flag = false;
                    for(auto x: root){
                        if(choose((abs(x) - REAL(1) < 0) || (0 < abs(x) - REAL(1)), (abs(x) - REAL(1) < power(2, -prec)) && (-power(2, -prec) < abs(x) - REAL(1))) == 1) {
                            non_unit = x;
                            flag = true;
                            break;
                        }
                    }
                    if(!flag){
                        cout << "no non unit norm for " << p << " / " << q << "\n";
                        continue;
                    }
                    
                    REAL gamma = sin(REAL(MATDIM*p*pi())/REAL(q))/sin(REAL(p*pi())/REAL(q));
    
                    // Reducibility Testing 
    
                    if(!irr){
                        COMPLEX center = -exp(COMPLEX(0, REAL((MATDIM+1) * p) * pi() / REAL(q))) / gamma;
                        REAL dist = abs(non_unit - center);
                        REAL disksize = sqrt((REAL(1)/(gamma*gamma)) - 1);
    
                        if(choose(dist - disksize <= 0, dist - disksize > -pow(2, -prec)) == 1){
                            cout << p << " / " << q << " failed irreducibility test\n";
                        }
                        else{
                            cout << p << " / " << q << " passed irreducibility test!\n";
                            irr = 1;
                        }
                    }
    
                    // Exchangablity Testing 
    
                    if(!exch){
                        REAL cond = abs(COMPLEX(1) - non_unit) - (1 + abs(non_unit))*sqrt(1 - gamma*gamma);
                        if(choose(cond <= 0, cond > -pow(2, -prec)) == 1){
                            cout << p << " / " << q << " failed unexchangeability test\n";
                        }
                        else{
                            cout << p << " / " << q << " passed unexchangeability test!\n";
                            exch = 1;
                        }
                    }
    
                    if(irr && exch) goto brk;
                }
            }
            brk:
            continue;
        }
    }
}
