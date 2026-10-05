// Same scalar rounding order as the notes. Compile without FMA or fast-math.
#include <cmath>
#include <cstdint>
#include <algorithm>
template<class T> void initialize(T* a, int n, int q, int p) {
    int k=1<<(2*q), b=7*k;
    T t=std::ldexp(T(3),-(p+q+2)/2);
    T eta=std::ldexp(T(12),q-p), delta=T(k)*eta;
    std::fill(a,a+int64_t(n)*n,T(0)); // Touch every page, including padding.
    auto put=[&](int i,int j,T x){a[int64_t(j)*n+i]=x; a[int64_t(i)*n+j]=x;};
    for(int i=0;i<n;++i) put(i,i,T(4));
    for(int j=0;j<2*k;++j) {
        for(int i=0;i<6*k;++i) put(b+j,i,T(2)*t);
        for(int i=0;i<k;++i) {
            T x=j<k ? T(i==j) : std::ldexp(T((__builtin_popcount(unsigned(i&(j-k)))&1)?-1:1),-q);
            put(b+j,6*k+i,T(3)*x);
        }
    }
    for(int j=0;j<2*k;++j) for(int i=j;i<2*k;++i) {
        T base=i==j ? T(2.25) : T(0);
        if(i>=k && j<k) base=std::ldexp(T((__builtin_popcount(unsigned((i-k)&j))&1)?-2.25:2.25),-q);
        T x=base+eta;
        if(i==j) x+=delta;
        put(b+i,b+j,x);
    }
}
template<class T> int factor(T* a,int n,int q,int p,int64_t* mismatches) {
    int k=1<<(2*q), b=7*k;
    T eta=std::ldexp(T(12),q-p);
    *mismatches=-1;
    for(int r=0;r<n;++r) {
        if(r==b) {
            *mismatches=0;
            for(int j=0;j<2*k;++j) for(int i=j;i<2*k;++i) {
                T expected=i==j ? T(k+1)*eta : ((i>=k && j<k)?-eta:-eta/T(8));
                *mismatches+=a[int64_t(b+j)*n+b+i]!=expected;
            }
        }
        T* col=a+int64_t(r)*n;
        if(!(col[r]>T(0)) || !std::isfinite(col[r])) return r+1;
        col[r]=std::sqrt(col[r]);
        // Skip zero contributions without regrouping effective operations.
        int end=n;
        while(end>r+1 && col[end-1]==T(0)) --end;
        for(int i=r+1;i<end;++i) col[i]/=col[r];
        for(int j=r+1;j<end;++j) {
            T x=col[j];
            if(x==T(0)) continue;
            T* dest=a+int64_t(j)*n;
            for(int i=j;i<end;++i) {
                T product=col[i]*x;
                dest[i]=dest[i]-product;
            }
        }
    }
    return 0;
}
extern "C" {
void init64(double* a,int n,int q){initialize(a,n,q,53);}
void init32(float* a,int n,int q){initialize(a,n,q,24);}
int chol64(double* a,int n,int q,int64_t* m){return factor(a,n,q,53,m);}
int chol32(float* a,int n,int q,int64_t* m){return factor(a,n,q,24,m);}
}
