#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <algorithm>
#include <cassert>
#include <string>
#include <vector>
#include <random>
using namespace std;

default_random_engine generator;
uniform_real_distribution<double> distribution_float(0.0, 1.0);
uniform_int_distribution<int> distribution_int(0,1);

class F2n_2{
public:
    int n;
    bool *v; // size = 2n

    F2n_2(){
        n = 0;
        v = NULL;
    }

    F2n_2(int _n){
        n = _n;
        v = (bool*)malloc((2 * n) * sizeof(bool));
        memset(v, 0, (2 * n) * sizeof(bool));
    }

    F2n_2(int _n, bool *_v){
        n = _n;
        v = _v;
    }

    F2n_2(F2n_2 &u, bool copy){
        n = u.n;

        if(copy){
            v = (bool*)malloc((2 * n) * sizeof(bool));
            for(int i = 0; i < 2 * n; i ++)
                v[i] = u.entry(i);
        }
        else{
            v = u.v;
        }
    }

    void free_class(){
        free(v);
    }

    bool operator ==(F2n_2 u){
        if(u.n != n) return false;
        for(int i = 0; i < 2 * n; i ++){
            if(u.v[i] != v[i]) return false;
        }
        return true;
    }

    void init(){
        memset(v, 0, (2 * n) * sizeof(bool));
    }

    void random(){
        for(int i = 0; i < 2 * n; i ++){
            v[i] = distribution_int(generator);
        }
    }

    void random_nonzero(){
        bool all_zero = true;
        while(all_zero){
            for(int i = 0; i < 2 * n; i ++){
                v[i] = distribution_int(generator);
                if(v[i] != 0)
                    all_zero = false;
            }
        }
    }

    bool& entry(int i){
        assert(i < 2 * n && i >= 0);
        return v[i];
    }

    void print(){
        printf("Printing an element in F2n_2:\n");
        for(int i = 0; i < 2 * n; i ++){
            printf("%d ", v[i]);
        }
        printf("\n");
    }
};

bool inner(F2n_2 v, F2n_2 w){ // symplectic inner product
    assert(v.n == w.n);
    int t = 0;
    for(int i = 0; i < v.n; i ++){
        t += v.entry(2 * i) * w.entry(2 * i + 1);
        t += w.entry(2 * i) * v.entry(2 * i + 1);
    }
    return t % 2;
}

F2n_2 transvection(F2n_2 k, F2n_2 v){
    F2n_2 ret(v, true);

    bool inn = inner(k, v);
    if(inn == true){
        for(int i = 0; i < 2 * ret.n; i++){
            ret.entry(i) = (ret.entry(i) != k.entry(i));
        }
    }
    return ret;
}

void findtransvection(F2n_2 x, F2n_2 y, F2n_2 &h1, F2n_2 &h2){ // finds h1,h2 such that y = Zh1 Zh2 x
    h1.init();
    h2.init();

    if(x == y){
        return;
    }
    if(inner(x, y)){
        for(int i = 0; i < 2 * h1.n; i ++)
            h1.entry(i) = (x.entry(i) != y.entry(i));
        return;
    }

    F2n_2 z(x.n);
    for(int i = 0; i < x.n; i ++){
        int ii = 2 * i;
        if((x.entry(ii) + x.entry(ii+1)) != 0 && (y.entry(ii) + y.entry(ii+1)) != 0){
            z.entry(ii) = (x.entry(ii) != y.entry(ii));
            z.entry(ii+1) = (x.entry(ii+1) != y.entry(ii+1));

            if(z.entry(ii) + z.entry(ii+1) == 0){
                z.entry(ii+1) = 1;
                if(x.entry(ii) != x.entry(ii+1)){
                    z.entry(ii) = 1;
                }
            }
            for(int i = 0; i < 2 * h1.n; i ++)
                h1.entry(i) = x.entry(i) != z.entry(i);
            for(int i = 0; i < 2 * h2.n; i ++)
                h2.entry(i) = y.entry(i) != z.entry(i);
            return;
        }
    }
    //didn’t find a pair

    // first y==00 and x doesn’t
    for(int i = 0; i < x.n; i ++){
        int ii = 2 * i;
        if(((x.entry(ii) + x.entry(ii+1)) != 0) && ((y.entry(ii) + y.entry(ii+1)) == 0)){
            if(x.entry(ii) == x.entry(ii+1)){
                z.entry(ii+1) = 1;
            }
            else{
                z.entry(ii+1) = x.entry(ii);
                z.entry(ii) = x.entry(ii+1);
            }
            break;
        }
    }
    // finally x==00 and y doesn’t
    for(int i = 0; i < x.n; i ++){
        int ii = 2 * i;
        if(((x.entry(ii) + x.entry(ii+1)) == 0) && ((y.entry(ii) + y.entry(ii+1)) != 0)){
            if(y.entry(ii) == y.entry(ii+1)){
                z.entry(ii+1) = 1;
            }
            else{
                z.entry(ii+1) = y.entry(ii);
                z.entry(ii) = y.entry(ii+1);
            }
            break;
        }
    }

    for(int i = 0; i < 2 * h1.n; i ++)
        h1.entry(i) = x.entry(i) != z.entry(i);
    for(int i = 0; i < 2 * h2.n; i ++)
        h2.entry(i) = y.entry(i) != z.entry(i);

    z.free_class();
    return;
}

class symplectic_matrix{
public:
    int n;
    bool *S; // size = 2n x 2n

    symplectic_matrix(){
        n = 0;
        S = NULL;
    }

    symplectic_matrix(int _n){
        n = _n;
        S = (bool*)malloc((2 * n) * (2 * n) * sizeof(bool));
        memset(S, 0, (2 * n) * (2 * n) * sizeof(bool));
    }

    void free_class(){
        free(S);
    }

    bool& entry(int i, int j){
        assert(i < 2 * n && i >= 0);
        assert(j < 2 * n && j >= 0);
        return S[i * (2 * n) + j];
    }

    F2n_2 row_slice(int i){
        F2n_2 r(n);
        for(int j = 0; j < 2 * n; j++){
            r.entry(j) = S[i * (2 * n) + j];
        }
        return r;
    }

    bool verify_symplectic(){
        for(int i = 0; i < 2 * n; i++){
            for(int j = 0; j < 2 * n; j++){
                int desired = 0;
                if(i / 2 == j / 2){
                    desired = (i % 2 != j % 2);
                }

                int sum = 0;
                for(int k = 0; k < n; k++){
                    sum += entry(i, 2*k) * entry(j, 2*k+1) + entry(i, 2*k+1) * entry(j, 2*k);
                }
                if(sum%2 != desired) return false;
            }
        }
        return true;
    }

    void print(){
        printf("Printing a symplectic matrix:\n");
        for(int i = 0; i < 2 * n; i++){
            for(int j = 0; j < 2 * n; j++){
                printf("%d ", entry(i, j));
            }
            printf("\n");
        }
    }
};

// Destroys m1 and m2 after returning the output
symplectic_matrix directsum(symplectic_matrix m1, symplectic_matrix m2){
    symplectic_matrix out(m1.n + m2.n);

    for(int i = 0; i < 2 * m1.n; i ++){
        for(int j = 0; j < 2 * m1.n; j ++){
            out.entry(i, j) = m1.entry(i, j);
        }
    }

    for(int i = 0; i < 2 * m2.n; i ++){
        for(int j = 0; j < 2 * m2.n; j ++){
            out.entry(i + 2 * m1.n, j + 2 * m1.n) = m2.entry(i, j);
        }
    }
    m1.free_class();
    m2.free_class();

    return out;
}

symplectic_matrix symplectic(int n){
    F2n_2 f1(n);
    f1.random_nonzero();

    F2n_2 e1(n);
    e1.entry(0) = 1;

    F2n_2 T0(n);
    F2n_2 T1(n);
    findtransvection(e1, f1, T0, T1);

    F2n_2 eprime(n);
    eprime.random();
    eprime.entry(0) = 1;
    eprime.entry(1) = 0;

    F2n_2 h0;
    h0 = transvection(T0, eprime);

    F2n_2 h0_new;
    h0_new = transvection(T1, h0);
    h0.free_class();

    if(distribution_int(generator) == 1)
        f1.init();

    symplectic_matrix id2(1);
    id2.entry(0, 0) = 1;
    id2.entry(1, 1) = 1;

    symplectic_matrix g;
    if(n != 1){
        g = directsum(id2, symplectic(n-1));
    }
    else{
        g = id2;
    }
    for(int i = 0; i < 2 * n; i ++){
        F2n_2 gi, gi_odd;
        gi = g.row_slice(i);
        gi_odd = transvection(T0, gi);
        gi.free_class();
        gi = transvection(T1, gi_odd);
        gi_odd.free_class();
        gi_odd = transvection(h0_new, gi);
        gi.free_class();
        gi = transvection(f1, gi_odd);
        gi_odd.free_class();
        for(int j = 0; j < 2 * n; j ++){
            g.entry(i, j) = gi.entry(j);
        }
        gi.free_class();
    }

    f1.free_class();
    e1.free_class();
    T0.free_class();
    T1.free_class();
    eprime.free_class();
    h0_new.free_class();

    return g;
}

class stabilizer_tableau{
public:
    int n;
    bool *x, *z, *r;

    stabilizer_tableau(int L){ // Toric code
        n = 2 * L * L;

        x = (bool*)malloc((2 * n + 1) * n * sizeof(bool));
        memset(x, 0, (2 * n + 1) * n * sizeof(bool));
        z = (bool*)malloc((2 * n + 1) * n * sizeof(bool));
        memset(z, 0, (2 * n + 1) * n * sizeof(bool));
        r = (bool*)malloc((2 * n + 1) * sizeof(bool));
        memset(r, 0, (2 * n + 1) * sizeof(bool));

        // First work out the stabilizers
        int cnt = n;

        // horizontal Z strings
        for(int i = 0; i < L; i++){
            int idx = L * (2 * L - 2) + i;
            z[cnt * n + (idx)] = 1;
        }
        cnt ++;

        // vertical Z strings
        for(int i = 1; i < 2 * L; i += 2){
            int idx = L * i + (L - 1);
            z[cnt * n + (idx)] = 1;
        }
        cnt ++;

        // Z plaquettes
        for(int _j = -1; _j < L-1; _j++){
            int j = (_j + L) % L;
            for(int _i = -1; _i < L-1; _i++){
                int i = (_i + L) % L;
                if(i == L-1 && j == L-1) continue;

                int idx1 = L * (2 * i) + j;
                z[cnt * n + (idx1)] = 1;
                int idx2 = L * (2 * i + 1) + j;
                z[cnt * n + (idx2)] = 1;
                int idx3 = L * (2 * i + 1) + ((j + 1) % L);
                z[cnt * n + (idx3)] = 1;
                int idx4 = L * ((2 * i + 2) % (2 * L)) + j;
                z[cnt * n + (idx4)] = 1;

                cnt ++;
            }
        }

        // X stars
        for(int j = L-1; j >= 0; j--){
            for(int i = L-1; i >= 0; i--){
                if(i == L-1 && j == L-1) continue;

                int idx1 = L * ((2 * i + 2 * L - 1) % (2 * L)) + j;
                x[cnt * n + (idx1)] = 1;
                int idx2 = L * (2 * i) + ((j + L - 1) % L);
                x[cnt * n + (idx2)] = 1;
                int idx3 = L * (2 * i) + j;
                x[cnt * n + (idx3)] = 1;
                int idx4 = L * ((2 * i + 1) % (2 * L)) + j;
                x[cnt * n + (idx4)] = 1;

                cnt ++;
            }
        }

        assert(cnt == 2 * n);

        // Now work out the destabilizers
        cnt = 0;

        // vertical X strings
        for(int i = 0; i < 2 * L; i += 2){
            int idx = L * i + (L - 1);
            x[cnt * n + (idx)] = 1;
        }
        cnt ++;

        // horizontal X strings
        for(int i = 0; i < L; i++){
            int idx = L * (2 * L - 1) + i;
            x[cnt * n + (idx)] = 1;
        }
        cnt ++;

        vector<int> one_string;

        // strings of X errors
        for(int y = 0; y < L; y++){
            one_string.clear();

            for(int i = 0; i < y; i++){
                int idx = L * (2 * L - 1) + i;
                one_string.push_back(idx);
            }
            for(int x = 0; x < 2 * L - 3; x += 2){
                int idx = L * x + ((y + L - 1) % L);
                one_string.push_back(idx);
            }

            for(int i = max(y-1, 0); i < (int)one_string.size(); i++){
                for(int j = 0; j <= i; j++){
                    x[cnt * n + one_string[j]] = 1;
                }
                cnt ++;
            }
        }

        // strings of Z errors
        for(int y = 0; y < L; y++){
            one_string.clear();

            for(int i = 1; i <= y; i++){
                int idx = L * (2 * L - 2) + (L - 1 - i);
                one_string.push_back(idx);
            }
            for(int x = 2 * L - 3; x >= 0; x -= 2){
                int idx = L * x + (L - 1 - y);
                one_string.push_back(idx);
            }

            for(int i = max(y-1, 0); i < (int)one_string.size(); i++){
                for(int j = 0; j <= i; j++){
                    z[cnt * n + one_string[j]] = 1;
                }
                cnt ++;
            }
        }

        assert(cnt == n);
    }

    stabilizer_tableau(int _n, bool random){
        n = _n;
        x = (bool*)malloc((2 * n + 1) * n * sizeof(bool));
        memset(x, 0, (2 * n + 1) * n * sizeof(bool));
        z = (bool*)malloc((2 * n + 1) * n * sizeof(bool));
        memset(z, 0, (2 * n + 1) * n * sizeof(bool));
        r = (bool*)malloc((2 * n + 1) * sizeof(bool));
        memset(r, 0, (2 * n + 1) * sizeof(bool));

        if(random == false){
            for(int i = 0; i < n; i ++){
                x[i * n + i] = 1;
                z[(n + i) * n + i] = 1;
            }
        }
        else{
            symplectic_matrix random_g = symplectic(n);
            assert(random_g.verify_symplectic());
            //printf("is symplectic? %d\n", random_g.verify_symplectic());
            //random_g.print();

            for(int i = 0; i < 2 * n; i ++){
                if(i >= n)
                    r[i] = distribution_int(generator);
                else
                    r[i] = 0;
                for(int j = 0; j < n; j ++){
                    if(i < n){
                        x[i * n + j] = random_g.entry(2 * i + 1, 2 * j + 1);
                        z[i * n + j] = random_g.entry(2 * i + 1, 2 * j);
                    }
                    else{
                        x[i * n + j] = random_g.entry(2 * (i - n), 2 * j + 1);
                        z[i * n + j] = random_g.entry(2 * (i - n), 2 * j);
                    }
                }
            }

            random_g.free_class();
        }
    }

    stabilizer_tableau(stabilizer_tableau const &ST, bool copy){
        if(copy){
            n = ST.n;
            x = (bool*)malloc((2 * n + 1) * n * sizeof(bool));
            memset(x, 0, (2 * n + 1) * n * sizeof(bool));
            z = (bool*)malloc((2 * n + 1) * n * sizeof(bool));
            memset(z, 0, (2 * n + 1) * n * sizeof(bool));
            r = (bool*)malloc((2 * n + 1) * sizeof(bool));
            memset(r, 0, (2 * n + 1) * sizeof(bool));

            for(int i = 0; i < (2 * n + 1) * n; i++){
                x[i] = ST.x[i];
                z[i] = ST.z[i];
            }
            for(int i = 0; i < 2 * n + 1; i++){
                r[i] = ST.r[i];
            }
        }
        else{
            n = ST.n;
            x = ST.x;
            z = ST.z;
            r = ST.r;
        }
    }

    void free_class(){
        free(x);
        free(z);
        free(r);
    }

    int g_func(int x1, int z1, int x2, int z2){
        if(x1 == 0 && z1 == 0){
            return 0;
        }
        if(x1 == 0 && z1 == 1){
            return x2 * (1 - 2 * z2);
        }
        if(x1 == 1 && z1 == 0){
            return z2 * (2 * x2 - 1);
        }
        if(x1 == 1 && z1 == 1){
            return z2 - x2;
        }
        assert(true);
        return 0;
    }

    void rowsum(int h, int i){
        if(h >= n){
            int val = 0;
            val = 2 * r[h] + 2 * r[i];
            for(int j = 0; j < n; j++){
                val += g_func(x[i * n + j], z[i * n + j], x[h * n + j], z[h * n + j]);
            }
            //printf("%d, %d val: %d\n", h, i, val);
            val = (val % 4 + 4) % 4;

            assert(val == 2 || val == 0);
            r[h] = (val == 2);
        }
        else{
            r[h] = 0;
        }

        for(int j = 0; j < n; j++){
            x[h * n + j] = (x[i * n + j] != x[h * n + j]);
            z[h * n + j] = (z[i * n + j] != z[h * n + j]);
        }
    }

    void CNOT(int a, int b){
        assert(a != b);
        for(int i = 0; i < 2 * n; i++){
            bool xia = x[i * n + a];
            bool xib = x[i * n + b];
            bool zia = z[i * n + a];
            bool zib = z[i * n + b];

            r[i] = (r[i] != (xia * zib * (xib != (zia != 1))));
            x[i * n + b] = (xib != xia);
            z[i * n + a] = (zia != zib);
        }
    }

    void Hadamard(int a){
        for(int i = 0; i < 2 * n; i++){
            bool xia = x[i * n + a];
            bool zia = z[i * n + a];

            r[i] = (r[i] != (xia * zia));
            z[i * n + a] = xia;
            x[i * n + a] = zia;
        }
    }

    void Phase(int a){
        for(int i = 0; i < 2 * n; i++){
            bool xia = x[i * n + a];
            bool zia = z[i * n + a];

            r[i] = (r[i] != (xia * zia));
            z[i * n + a] = (zia != xia);
        }
    }

    int measurement(int a, int forced_value = -999){
        int p = -1;
        for(int i = n; i < 2 * n; i++){
            if(x[i * n + a] == 1){
                p = i;
                break;
            }
        }
        //printf("%d %d\n", a, p);

        if(p == -1){
            for(int j = 0; j < n; j++){
                x[(2 * n) * n + j] = 0;
                z[(2 * n) * n + j] = 0;
            }
            r[2 * n] = 0;

            for(int i = 0; i < n; i++){
                if(x[i * n + a] == 1)
                    rowsum(2 * n, i + n);
            }
            return 2 * 0 + r[2 * n]; // deterministic, with outcome r[2 * n]
        }
        else{
            for(int i = 0; i < 2 * n; i++){
                if(i != p && x[i * n + a])
                    rowsum(i, p);
            }
            for(int j = 0; j < n; j++){
                x[(p - n) * n + j] = x[p * n + j];
                z[(p - n) * n + j] = z[p * n + j];
            }
            r[p - n] = r[p];
            for(int j = 0; j < n; j++){
                x[p * n + j] = 0;
                z[p * n + j] = 0;
            }
            z[p * n + a] = 1;

            if(forced_value != -999)
                r[p] = forced_value;
            else
                r[p] = distribution_int(generator);

            return 2 * 1 + r[p]; // random, with outcome r[p]
        }
    }

    void print(){
        printf("Printing a stabilizer state:\n");
        for(int i = 0; i < 2 * n; i ++){
            if(i == n){
                for(int j = 0; j < 2 * n + 3; j ++){
                    printf("--");
                }
                printf("\n");
            }

            for(int j = 0; j < n; j ++){
                printf("%d ", x[i * n + j]);
            }
            printf("| ");
            for(int j = 0; j < n; j ++){
                printf("%d ", z[i * n + j]);
            }
            printf("| ");
            printf("%d\n", r[i]);
        }
    }
};

int L;

vector<int> seq_singleQ_Clifford[24] = {{}, {0}, {1}, {0, 1}, {1, 0}, {1, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}, {0, 1, 0, 1}, {0, 1, 1, 0}, {0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1}, {0, 1, 0, 1, 1}, {0, 1, 1, 0, 1}, {1, 0, 1, 1, 0}, {1, 0, 1, 1, 1}, {1, 1, 0, 1, 1}, {0, 1, 0, 1, 1, 0}, {0, 1, 0, 1, 1, 1}, {0, 1, 1, 0, 1, 1}};

void apply_const_depth_Clifford_circuit(stabilizer_tableau& ST, int depth){
    for(int d = 0; d < depth; d++){
        for(int i = 0; i < ST.n; i++){
            int C = floor(distribution_float(generator) * 24);

            vector<int> seq = seq_singleQ_Clifford[C];
            reverse(seq.begin(), seq.end());
            for(int x: seq){
                if(x == 0)
                    ST.Hadamard(i);
                else
                    ST.Phase(i);
            }
        }

        if(d % 2 == 0){
            for(int i = 0; i < 2 * L; i++){
                for(int j = 0; j < L; j++){
                    if(j+1 >= L) continue;

                    double p = distribution_float(generator);
                    if(p < 1.0 / 3.0){
                        ST.CNOT(i*L+j, i*L+j+1);
                    }
                    else if(p < 2.0 / 3.0){
                        ST.CNOT(i*L+j+1, i*L+j);
                    }
                }
            }
        }
        else{
            for(int i = 0; i < 2 * L; i++){
                for(int j = 0; j < L; j++){
                    if(i+1 >= 2 * L) continue;

                    double p = distribution_float(generator);
                    if(p < 1.0 / 3.0){
                        ST.CNOT(i*L+j, (i+1)*L+j);
                    }
                    else if(p < 2.0 / 3.0){
                        ST.CNOT((i+1)*L+j, i*L+j);
                    }
                }
            }
        }
        // for(int i = (d % 2); i < ST.n; i+=2){
        //     if(i+1 >= ST.n) continue;
        //
        //     double p = distribution_float(generator);
        //     if(p < 1.0 / 3.0){
        //         ST.CNOT(i, i+1);
        //     }
        //     else if(p < 2.0 / 3.0){
        //         ST.CNOT(i+1, i);
        //     }
        // }
    }
}

vector<int> classical_shadow(stabilizer_tableau& ST){
    stabilizer_tableau ST_m(ST, true);
    vector<int> shadow;

    for(int i = 0; i < ST.n; i++){
        int outcome = 0;
        double p = distribution_float(generator);
        if(p < 1.0 / 3.0){ // X basis
            ST_m.Hadamard(i);
            outcome += 2;
        }
        else if(p < 2.0 / 3.0){ // Y basis
            ST_m.Phase(i);
            ST_m.Hadamard(i);
            outcome += 4;
        }
        outcome += (ST_m.measurement(i) % 2);

        shadow.push_back(outcome);
    }
    return shadow;
}

vector<vector<vector<int> > > training_data;

int main(int argc, char *argv[]){
    assert(argc == 7);

    L = stoi(argv[1]);
    int Nsample = stoi(argv[2]);
    int Ntrain = stoi(argv[3]);
    int depth = stoi(argv[4]);

    generator = default_random_engine(stoi(argv[5]));

    float gamma = stof(argv[6]);

    int n = 2 * L * L;

    printf("[");
    for(int dep = 0; dep < depth; dep++){
        training_data.clear();

        stabilizer_tableau ToricCode(L);
        for(int t = 0; t < Ntrain / 2; t++){
            stabilizer_tableau ToricPhase(ToricCode, true);
            apply_const_depth_Clifford_circuit(ToricPhase, dep);

            vector<vector<int> > one_training_data;
            for(int s = 0; s < Nsample; s++){
                one_training_data.push_back(classical_shadow(ToricPhase));
            }
            training_data.push_back(one_training_data);
        }

        for(int t = Ntrain / 2; t < Ntrain; t++){
            stabilizer_tableau TrivialPhase(n, false);
            apply_const_depth_Clifford_circuit(TrivialPhase, dep);

            vector<vector<int> > one_training_data;
            for(int s = 0; s < Nsample; s++){
                one_training_data.push_back(classical_shadow(TrivialPhase));
            }
            training_data.push_back(one_training_data);
        }

        // True kernel
        printf("[");
        for(int t1 = 0; t1 < Ntrain; t1++){
            printf("[");
            for(int t2 = 0; t2 < Ntrain; t2++){
                float kernel = 0;
                for(int r1 = 0; r1 < Nsample; r1++){
                    for(int r2 = 0; r2 < Nsample; r2++){
                        if(r1 == r2) continue;

                        vector<float> list_aff(n+1);
                        list_aff[0] = 0;
                        for(int i = 0; i < n; i++){
                            list_aff[i+1] = list_aff[i] + (training_data[t1][r1][i]/2 == training_data[t2][r2][i]/2 ? ((training_data[t1][r1][i]%2 == training_data[t2][r2][i]%2)? 5: -4) : 0.5);
                        }

                        float affinity = 0;
                        affinity = exp(gamma * list_aff[n] / n);

                        kernel += affinity / Nsample / (Nsample - 1);
                    }
                }
                printf("%.3f", kernel);
                if(t2 < Ntrain - 1) printf(", ");
            }
            printf("],\n");
        }
        printf("],\n");
    }
    printf("]\n");
}
