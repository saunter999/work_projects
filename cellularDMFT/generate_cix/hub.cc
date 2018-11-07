#define NO_SYMMETRY // superconducting state
//#define MOMENTUM_SYMMETRY

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <algorithm>
#include "sutil.h"
#include "complex.h"
#include "sblas.h"
#include "sfunction.h"

using namespace std;

// Just old good !
int Binomial(int n, int m)
{
  double Mf = 1;
  for (int i=2; i<=m; i++) Mf *= i;
  double r = 1;
  for (int i=n; i>=n-m+1; i--) r*=i;
  return static_cast<int>(r/Mf);
}
double dFactorial(int j)
{
  if (j<0) {cerr<<"Factorial defined only for non-negative numbers!"<<endl;}
  if (j<=1) return 1;
  double x=1;
  for (int i=2; i<=j; i++) x*=i;
  return x;
}
double dFactorial(double x)
{
  double r=1;
  double y=x;
  while(y>1.0){
    r *= y;
    y--;
  }
  if (fabs(y-0.5)<1e-10) r/=M_2_SQRTPI;
  return r;
}

class PrimitiveShift{
public:
  // shift and inverse shift
  vector<int> shx, shy, ishx, ishy;
  vector<pair<double,double> > K;
  vector<pair<int,int> > shxy;
  vector<vector<dcomplex> > expiK;
public:
  PrimitiveShift(){};
  void Init(int N)
  {
    shx.resize(N);
    shy.resize(N);
    ishx.resize(N);
    ishy.resize(N);
    K.resize(N);
    shxy.resize(N);
    expiK.resize(N);
    if (N==4){
      //  2 3
      //  0 1
      shxy[0] = make_pair(0,0);
      shxy[1] = make_pair(1,0);
      shxy[2] = make_pair(0,1);
      shxy[3] = make_pair(1,1);
      shx[0] = 1;
      shx[1] = 0;
      shx[2] = 3;
      shx[3] = 2;
      shy[0] = 2;
      shy[1] = 3;
      shy[2] = 0;
      shy[3] = 1;
      K[0] = make_pair(0,   0);
      K[1] = make_pair(M_PI,0);
      K[2] = make_pair(0,   M_PI);
      K[3] = make_pair(M_PI,M_PI);
      for (int k=0; k<N; k++){
	expiK[k].resize(N);
	for (int s=0; s<N; s++){
	  double phi = shxy[s].first*K[k].first + shxy[s].second*K[k].second;
	  expiK[k][s] = dcomplex(cos(phi),sin(phi));
	}
      }
    } else if (N==8){
      //      7
      //    4 6 0
      //  7 2 3 5 7
      //  6 0 1 4 6
      //    5 7 2
      //      6
      shxy[0] = make_pair(0,0);
      shxy[1] = make_pair(1,0);
      shxy[2] = make_pair(0,1);
      shxy[3] = make_pair(1,1);
      shxy[4] = make_pair(2,0);
      shxy[5] = make_pair(2,1);
      shxy[6] = make_pair(1,2);
      shxy[7] = make_pair(1,3);// (1,-1)
      shx[0] = 1;
      shx[1] = 4;
      shx[2] = 3;
      shx[3] = 5;
      shx[4] = 6;
      shx[5] = 7;
      shx[6] = 0;
      shx[7] = 2;
      shy[0] = 2;
      shy[1] = 3;
      shy[2] = 4;
      shy[3] = 6;
      shy[4] = 5;
      shy[5] = 0;
      shy[6] = 7;
      shy[7] = 1;
      K[0] = make_pair(0,   0);
      K[1] = make_pair(M_PI,0);
      K[2] = make_pair(0,   M_PI);
      K[3] = make_pair(M_PI,M_PI);
      K[4] = make_pair( 0.5*M_PI, 0.5*M_PI);
      K[5] = make_pair(-0.5*M_PI, 0.5*M_PI);
      K[6] = make_pair( 0.5*M_PI,-0.5*M_PI);
      K[7] = make_pair(-0.5*M_PI,-0.5*M_PI);
      for (int k=0; k<N; k++){
	expiK[k].resize(N);
	for (int s=0; s<N; s++){
	  double phi = shxy[s].first*K[k].first + shxy[s].second*K[k].second;
	  expiK[k][s] = dcomplex(cos(phi),sin(phi));
	}
      }
    }
    for (int i=0; i<shx.size(); i++) ishx[shx[i]]=i;
    for (int i=0; i<shy.size(); i++) ishy[shy[i]]=i;
  }
  PrimitiveShift(int N)
  { Init(N);}
  int Plus(int k1, int k2){
    double kx = K[k1].first + K[k2].first;
    double ky = K[k1].second + K[k2].second;
    pair<double,double> Kn = make_pair(fmod(kx,2*M_PI),fmod(ky,2*M_PI));
    if (Kn.first>M_PI) Kn.first -= 2*M_PI;
    if (Kn.second>M_PI) Kn.second -= 2*M_PI;
    if (Kn.first<=-M_PI) Kn.first += 2*M_PI;
    if (Kn.second<=-M_PI) Kn.second += 2*M_PI;
    for (int i=0; i<K.size(); i++){
      if (sqr(Kn.first-K[i].first)+sqr(Kn.second-K[i].second)<1e-6) return i;
    }
    cerr<<"Should't happen! "<<kx/M_PI<<" "<<ky/M_PI<<" "<<k1<<" "<<k2<<endl;
    return 0;
  }
  int eqK(int k){
#ifdef NO_SYMMETRY
    return k;
#else
    if (K.size()==4){
      switch(k){
      case 0: return 0;
      case 1: return 1;
      case 2: return 1;
      case 3: return 2;
      default:  cerr<<"Shuldn't happen!"<<endl; return 0;
      }
    }
    if (K.size()==8){
      switch(k){
      case 0: return 0;
      case 1: return 1;
      case 2: return 1;
      case 3: return 2;
      case 4: return 3;
      case 5: return 3;
      case 6: return 3;
      case 7: return 3;
      default:  cerr<<"Shuldn't happen!"<<endl; return 0;
      }
    }
#endif      
  }
  int nK(){
#ifdef NO_SYMMETRY
    return K.size();
#else	
    if (K.size()==4) return 3;
    if (K.size()==8) return 4;
#endif    
  }
  int nK(int k){
#ifndef NO_SYMMETRY    
    if (K.size()==4 && k==1) return 2;
    if (K.size()==8 && k==1) return 2;
    if (K.size()==8 && k==3) return 4;
#endif    
    return 1;
  }
  bool equivalent(int K1, int K2)
  {
    if (K1==K2) return true;
#ifdef NO_SYMMETRY
    return false;
#else
    if (K.size()==4) if (K1==1 && K2==2 || K1==2 && K2==1) return true; else return false;
    if (K.size()==8) if (K1==1 && K2==2 || K1==2 && K2==1 || K1>3 && K2>3) return true; else return false;
#endif    
  }
  int minus(int k){
    if (K.size()==4) return k;
    if (K.size()==8)
      switch(k){
      case 4: return 7;
      case 5: return 6;
      case 6: return 5;
      case 7: return 4;
      default: return k;
      }
    return k;
  }
};


class operateLS{
  int N, baths;
  vector<int> mask, mask_u, mask_d;
  vector<int> sz;
  void Init()
  {
    for (int i=0; i<baths; i++) mask[i] = 1<<i;
    for (int i=0; i<N; i++){mask_d[i] = mask[i];}
    for (int i=N; i<baths; i++){mask_u[i-N] = mask[i];}
    for (int i=0; i<N; i++){ sz[i]=-1; sz[i+N]=1;}
  }
public:
  operateLS() : N(0){};
  operateLS(int N_) : N(N_), baths(2*N), mask(2*N), mask_u(N), mask_d(N), sz(baths) {Init();}
  void Init(int N_)
  {
    N=N_;
    baths = 2*N;
    mask.resize(2*N);
    mask_u.resize(N);
    mask_d.resize(N);
    sz.resize(baths);
    Init();
  }
  int estate(int state, int i) const{
    if (N==0) {cerr<<"Did not initialize operateLS"<<endl; exit(1);}
    bool lup = state & mask_u[i];
    bool ldo = state & mask_d[i];
    if (lup && ldo) return 2;
    if (lup) return 1;
    if (ldo) return -1;
    return 0;
  }
  int Nel(int state) const
  {
    int n=0;
    for (int k=0; k<baths; k++) if (mask[k]&state) n++;
    return n;
  }
  double Sz(int state) const
  {
    int nu=0, nd=0;
    for (int i=0; i<N; i++){
      if (mask_u[i]&state) nu++;
      if (mask_d[i]&state) nd++;
    }
    return 0.5*(nu-nd);
  }
  int Doubles(int state) const
  {
    int d=0;
    for (int i=0; i<N; i++)
      if (mask_u[i]&state && mask_d[i]&state) d++;
    return d;
  }
  int sign(int state, int mask_min, int mask_max) const
  {
    int mask = mask_min<<1;
    int n=0;
    while (mask<mask_max){
      if (mask&state) n++;
      mask = mask<<1;
    }
    return 1-2*(n%2);
  }
  int sign(int state, int mask) const
  {
    int m = 1;
    int n=0;
    while (m<mask){
      if (m&state) n++;
      m = m<<1;
    }
    return 1-2*(n%2);
  }
  void S2(int state, list<pair<int,double> >& sts) const
  {// S^2=S_z^2+1/2\sum_{i}(S_i^+ S_i^- + S_i^- S_i^+) + \sum_{i!=j}(S_i^- S_j^+)
    sts.clear();
    // diagonal part
    double dd=0;
    for (int i=0; i<N; i++){
      int up = (mask_u[i]&state) ? 1 : 0;
      int dn = (mask_d[i]&state) ? 1 : 0;
      // if only up or only down at certain side
      if (up+dn==1) dd += 0.5;
    }
    // Sz^2
    double fct = sqr(Sz(state)) + dd;
    // store diagonal
    sts.push_back(make_pair(state,fct));
    // off diagonal
    for (int i=0; i<N; i++){
      int im1 = mask_u[i];
      int im2 = mask_d[i];
      bool ib1 = state & im1;
      bool ib2 = state & im2;
      if (ib1 && !ib2){// S^-_i gives nonzero
	int isig = sign(state,im2,im1);
	int istate = state^im1^im2;
	for (int j=0; j<N; j++){
	  if (i==j) continue;
	  int jm1 = mask_d[j];
	  int jm2 = mask_u[j];
	  bool jb1 = state & jm1;
	  bool jb2 = state & jm2;
	  if (jb1 && !jb2){//S^+_j gives nonzero
	    int jsig = sign(istate,jm1,jm2);
	    int jstate = istate^jm1^jm2;
	    sts.push_back(make_pair(jstate, isig*jsig));
	  }
	}
      }
    }
  }
  int sign(vector<int>& a){
    int n=a.size();
    int s=0;
    for (int i=0; i<n-1; i++) {
      for (int j=0; j<n-1-i; j++)
	if (a[j+1] < a[j]) {  /* compare the two neighbors */
	  int tmp = a[j];         /* swap a[j] and a[j+1]      */
	  a[j] = a[j+1];
	  a[j+1] = tmp;
	  s++;
	}
    }
    return 1-2*(s%2);
  }
  void Shifts(int state, int nshx, const vector<int>& shx, int nshy, const vector<int>& shy, int& new_state, int& sgn)
  {
    vector<int> orderd, orderu;
    int newstate=0;
    for (int i=0; i<N; i++){
      if (state&mask_d[i]){
	int j=i;
	for (int l=0; l<nshx; l++) j = shx[j];
	for (int l=0; l<nshy; l++) j = shy[j];
	newstate = newstate^mask_d[j];
	orderd.push_back(j);
      }
      if (state&mask_u[i]){
	int j=i;
	for (int l=0; l<nshx; l++) j = shx[j];
	for (int l=0; l<nshy; l++) j = shy[j];
	newstate = newstate^mask_u[j];
	orderu.push_back(j);
      }
    }
    sgn = sign(orderd)*sign(orderu);
    new_state = newstate;
  }
  void Hop(int state, const PrimitiveShift& s, map<int,double>& sts, map<int,double>& stsp) const
  {
    sts.clear(), stsp.clear();
    vector<int> sosed(8);
    for (int i=0; i<N; i++){
      sosed[0] = s.shx[i];
      sosed[1] = s.ishx[i];
      sosed[2] = s.shy[i];
      sosed[3] = s.ishy[i];
      sosed[4] = s.shy[s.shx[i]];
      sosed[5] = s.shy[s.ishx[i]];
      sosed[6] = s.ishy[s.shx[i]];
      sosed[7] = s.ishy[s.ishx[i]];
      if (state&mask_d[i]){// hops from i
	for (int l=0; l<sosed.size(); l++){
	  int j = sosed[l]; // hops to j
	  if (!(state&mask_d[j])){// Pauli principle
	    int newstate = state^mask_d[i]^mask_d[j];
	    int sgn = sign(state,min(mask_d[i],mask_d[j]),max(mask_d[i],mask_d[j]));
	    if (l<4) sts[newstate] += sgn; // hopping with t
	    else stsp[newstate] += sgn;    // hopping with t'
	  }
	}
      }
      if (state&mask_u[i]){// hops from i
	for (int l=0; l<sosed.size(); l++){
	  int j = sosed[l]; // hops to j
	  if (!(state&mask_u[j])){// Pauli principle
	    int newstate = state^mask_u[i]^mask_u[j];
	    int sgn = sign(state,min(mask_u[i],mask_u[j]),max(mask_u[i],mask_u[j]));
	    if (l<4) sts[newstate] += sgn; // hopping with t
	    else stsp[newstate] += sgn;    // hopping with t'
	  }
	}
      }
    }
  }
  void Flip(int state, const PrimitiveShift& s, map<int,double>& sts, map<int,double>& stsp) const
  {// \sum_{ij} Si.Sj = \sum_{ij} (Si^z Sj^z + Si^- Sj^+)
    sts.clear(); stsp.clear();
    vector<int> sosed(8);
    for (int i=0; i<N; i++){
      sosed[0] = s.shx[i];
      sosed[1] = s.ishx[i];
      sosed[2] = s.shy[i];
      sosed[3] = s.ishy[i];
      sosed[4] = s.shy[s.shx[i]];
      sosed[5] = s.shy[s.ishx[i]];
      sosed[6] = s.ishy[s.shx[i]];
      sosed[7] = s.ishy[s.ishx[i]];
      int mid = state&mask_d[i];
      int miu = state&mask_u[i];
      if ((mid && !miu) || (miu && !mid)){// something at i
	int isz = (miu) ? 1 : -1;
	for (int l=0; l<sosed.size(); l++){
	  int j = sosed[l]; // flips with j
	  int mjd = state&mask_d[j];
	  int mju = state&mask_u[j];
	  if ((mjd && !mju) || (mju && !mjd)){
	    int jsz = (mju) ? 1 : -1;
	    if (l<4) sts[state] += 0.25*isz*jsz;
	    else stsp[state] += 0.25*isz*jsz;
	    if (isz==1 && jsz==-1){// flip possible
	      int newstate = state^mask_u[i]^mask_d[j]^mask_u[j]^mask_d[i];
	      int sgn1 = sign(state,mask_d[i],mask_u[i]);
	      int tstate = state^mask_u[i]^mask_d[i];
	      int sgn2 = sign(tstate,mask_d[j],mask_u[j]);
	      if (l<4) sts[newstate] += sgn1*sgn2;
	      else stsp[newstate] += sgn1*sgn2;
	    }
	  }
	}
      }
    }
  }
  void Sp(int dS, int state, map<int,double>& sts0) const
  {
    sts0.clear();
    sts0[state] = 1.0;
    for (int n=0; n<dS; n++){
      map<int,double> sts;
      for (map<int,double>::const_iterator it=sts0.begin(); it!=sts0.end(); it++){
	int state = it->first;
	double fact = it->second;
	for (int i=0; i<N; i++){
	  if (mask_d[i]&state && !(mask_u[i]&state)){
	    int newstate = state^mask_d[i]^mask_u[i];
	    int sgn = sign(state,mask_d[i],mask_u[i]);
	    sts[newstate] += fact*sgn;
	  }
	}
      }
      sts0=sts;
    }
    for (int n=0; n>dS; n--){
      map<int,double> sts;
      for (map<int,double>::const_iterator it=sts0.begin(); it!=sts0.end(); it++){
	int state = it->first;
	double fact = it->second;
	for (int i=0; i<N; i++){
	  if (mask_u[i]&state && !(mask_d[i]&state)){
	    int newstate = state^mask_u[i]^mask_d[i];
	    int sgn = sign(state,mask_d[i],mask_u[i]);
	    sts[newstate] += fact*sgn;
	  }
	}
      }
      sts0=sts;
    }
  }
  void Fp(int state, list<pair<int,dcomplex> >& sts, const PrimitiveShift& sh, int K, int ud) const
  {// K=0...N, ud=+1,-1
    sts.clear();
    if (ud>0){
      for (int i=0; i<N; i++){
	if (!(state&mask_u[i])){
	  dcomplex fct = sign(state,mask_u[i])*sh.expiK[K][i]/sqrt(static_cast<double>(N));
	  sts.push_back(make_pair(state^mask_u[i], fct));
	}
      }
    } else{
      for (int i=0; i<N; i++){
	if (!(state&mask_d[i])){
	  dcomplex fct = sign(state,mask_d[i])*sh.expiK[K][i]/sqrt(static_cast<double>(N));
	  sts.push_back(make_pair(state^mask_d[i], fct));
	}
      }
    }
  }
  void Sz(int state, dcomplex& fct, const PrimitiveShift& sh, int K) const
  {
    fct = 0;
    for (int i=0; i<N; i++){
      if (state&mask_u[i])
	fct += 0.5*sh.expiK[K][i]/sqrt(static_cast<double>(N));
      if (state&mask_d[i])
	fct -= 0.5*sh.expiK[K][i]/sqrt(static_cast<double>(N));
    }
  }
  void Nks(int state, int updown, dcomplex& fct, const PrimitiveShift& sh, int K) const
  {
    fct = 0;
    for (int i=0; i<N; i++){
      if (updown==1){
	if (state&mask_u[i]) fct += sh.expiK[K][i]/sqrt(static_cast<double>(N));
      }else{
	if (state&mask_d[i]) fct += sh.expiK[K][i]/sqrt(static_cast<double>(N));
      }
    }
  }

  void Hop_down(int i, int j, int state, int& newstate, int& sgn) const
  {
    newstate=-1; sgn=0;
    if (state&mask_d[i]){// hops from i
      if (!(state&mask_d[j])){// Pauli principle
	newstate = state^mask_d[i]^mask_d[j];
	sgn = sign(state,min(mask_d[i],mask_d[j]),max(mask_d[i],mask_d[j]));
      }
    }
  }
  void Hop_up(int i, int j, int state, int& newstate, int& sgn) const
  {
    newstate=-1; sgn=0;
    if (state&mask_u[i]){// hops from i
      if (!(state&mask_u[j])){// Pauli principle
	newstate = state^mask_u[i]^mask_u[j];
	sgn = sign(state,min(mask_u[i],mask_u[j]),max(mask_u[i],mask_u[j]));
      }
    }
  }
  void Hop(int state, const PrimitiveShift& s, map<int,dcomplex>& sts, int K) const
  {
    sts.clear();
    for (int i=0; i<N; i++){
      for (int l=0; l<N; l++){
	int sx = s.shxy[l].first;
	int sy = s.shxy[l].second;
	int j=i;
	for (int m=0; m<sx; m++) j = s.shx[j];
	for (int m=0; m<sy; m++) j = s.shy[j];
	dcomplex pf = s.expiK[K][l];
	int newstate=-1, sgn=0;
	Hop_up(i,j,state,newstate,sgn);
	if (newstate!=-1 && sgn!=0) sts[newstate] += pf*sgn/N;
	Hop_down(i,j,state,newstate,sgn);
	if (newstate!=-1 && sgn!=0) sts[newstate] += pf*sgn/N;
      }
    }
  }
  void Sp(int state, map<int,dcomplex>& sts, const PrimitiveShift& sh, int K) const
  {
    sts.clear();
    dcomplex fct = 0;
    for (int i=0; i<N; i++){
      if ((state&mask_d[i]) && !(state&mask_u[i])){
	int newstate = state^mask_d[i]^mask_u[i];
	int sgn = sign(state,mask_d[i],mask_u[i]);
	sts[newstate] += sgn*sh.expiK[K][i]/sqrt(static_cast<double>(N));
      }
    }
  }
  void Sm(int state, map<int,dcomplex>& sts, const PrimitiveShift& sh, int K) const
  {
    dcomplex fct = 0;
    for (int i=0; i<N; i++){
      if (!(state&mask_d[i]) && (state&mask_u[i])){
	int newstate = state^mask_d[i]^mask_u[i];
	int sgn = sign(state,mask_d[i],mask_u[i]);
	sts[newstate] += sgn*sh.expiK[K][i]/sqrt(static_cast<double>(N));
      }
    }
  }
  string print(int state) const
  {
    stringstream stream;
    for (int i=0; i<N; i++)
      stream<<setw(3)<<estate(state,i)<<" ";
    return stream.str();
  }
  string printn(int state) const
  {
    stringstream stream;
    for (int i=0; i<N; i++){
      int e = estate(state,i);
      char ch;
      switch(e){
      case  0: ch='0'; break;
      case  2: ch='2'; break;
      case  1: ch='u'; break;
      case -1: ch='d'; break;
      }
      stream<<ch;
    }
    return stream.str();
  }
};

void CreateBase_n_Sz(int Nm, int n, double Sz, const operateLS& op, vector<int>& base)
{
  int imax = 1<<(2*Nm);
  list<int> tbase;
  for (int i=0; i<imax; i++){
    if (/*op.Doubles(i)==0 &&*/ op.Nel(i)==n && op.Sz(i)==Sz){
      cout<<setw(3)<<i<<"  "<<op.printn(i)<<" "<<setw(3)<<op.Nel(i)<<"  "<<setw(5)<<op.Sz(i)<<endl;
      tbase.push_back(i);
    }
  }
  base.resize(tbase.size());
  int l=0;
  for (list<int>::const_iterator i=tbase.begin(); i!=tbase.end(); i++,l++) base[l] = *i;
}


void ComputeHamiltonian(double t0, double tp0, double U, const vector<int>& base, map<int,int>& ind, const operateLS& op, const PrimitiveShift& sh, function2D<dcomplex>& Ham)
{
  Ham=0;
  // Hopping
  for (int i=0; i<base.size(); i++){
    map<int,double> sts, stsp;
    op.Hop(base[i], sh, sts,stsp);

    for (map<int,double>::const_iterator l=sts.begin(); l!=sts.end(); l++){
      if (ind.find(l->first)!=ind.end()) // otherwise is not part of the base
	Ham(i,ind[l->first]) += -t0*l->second;
    }
    for (map<int,double>::const_iterator l=stsp.begin(); l!=stsp.end(); l++){
      if (ind.find(l->first)!=ind.end()) // otherwise is not part of the base
	Ham(i,ind[l->first]) += -tp0*l->second;
    }
  }
  // Hubbard U
  for (int i=0; i<base.size(); i++){
    int dbl = op.Doubles(base[i]);
    Ham(i,i) += dbl*U;
  }
  
//   // Flipping
//   for (int i=0; i<base.size(); i++){
//     map<int,double> sts, stsp;
//     op.Flip(base[i], sh, sts,stsp);

//     for (map<int,double>::const_iterator l=sts.begin(); l!=sts.end(); l++){
//       if (ind.find(l->first)!=ind.end()) // otherwise is not part of the base
// 	Ham(i,ind[l->first]) += 0.5*J0*l->second;
//     }
//     for (map<int,double>::const_iterator l=stsp.begin(); l!=stsp.end(); l++){
//       if (ind.find(l->first)!=ind.end()) // otherwise is not part of the base
// 	Ham(i,ind[l->first]) += 0.5*Jp0*l->second;
//     }
//   }
}

void ComputeS2(const vector<int>& base, map<int,int>& ind, const operateLS& op, function2D<dcomplex>& S2)
{
  S2=0;
  for (int i=0; i<base.size(); i++){
    list<pair<int,double> > sts;
    op.S2(base[i], sts);
    for (list<pair<int,double> >::const_iterator l=sts.begin(); l!=sts.end(); l++){
      if (ind.find(l->first)!=ind.end()) // otherwise is not part of the base
	S2(i,ind[l->first]) += l->second;
    }
  }
}

void ComputeSp(int dSz, const vector<int>& base, map<int,int>& ind, const vector<int>& basen, map<int,int>& indn,
	       const operateLS& op, function2D<dcomplex>& Sp)
{
  Sp.resize(base.size(),basen.size());
  Sp=0;
  for (int i=0; i<base.size(); i++){
    map<int,double> sts;
    op.Sp(dSz, base[i], sts);
    for (map<int,double>::const_iterator l=sts.begin(); l!=sts.end(); l++){
      if (indn.find(l->first)!=indn.end()) // otherwise is not part of the base
	Sp(i,indn[l->first]) += l->second;
    }
  }
}

void Test(int Nm, const function2D<dcomplex>& Ham, const function1D<int>& NumPsis)
{
  for (int i=0; i<Ham.size_N(); i++){
    for (int j=0; j<Ham.size_Nd(); j++){
      dcomplex t = Ham(i,j);
      if (norm(t)<1e-10) t=0;
      if (fabs(t.imag())<1e-10) t = dcomplex(t.real(),0);
      if (fabs(t.real())<1e-10) t = dcomplex(0,t.imag());
      cout<<setw(10)<<t<<" ";
    }
    cout<<endl;
  }
  
  for (int i=0; i<Ham.size_N(); i++)
    for (int j=0; j<Ham.size_Nd(); j++)
      if (norm(Ham(i,j)-Ham(j,i).conj())>1e-10) cerr<<"Not Hermitian "<<i<<" "<<j<<" "<<Ham(i,j)<<" "<<Ham(j,i)<<endl;

  int l=0;
  for (int k=0; k<Nm; k++){
    for (int i=l; i<l+NumPsis[k]; i++){
      for (int j=0; j<Ham.size_Nd(); j++){
	if (j>=l && j<l+NumPsis[k]) continue;
	if (norm(Ham(i,j))>1e-10) cerr<<"Off-diagonal elements in Hamiltonian! ERROR! "<<i<<" "<<j<<" "<<Ham(i,j)<<endl;
      }
    }
    for (int j=l; j<l+NumPsis[k]; j++){
      for (int i=0; i<Ham.size_N(); i++){
	if (i>=l && i<l+NumPsis[k]) continue;
	if (norm(Ham(i,j))>1e-10) cerr<<"Off-diagonal elements in Hamiltonian! ERROR! "<<i<<" "<<j<<" "<<Ham(i,j)<<endl;
      }
    }
    l += NumPsis[k];
  }
}
void Test(int Nm, const function2D<dcomplex>& Ham, const list<pair<int,int> >& start_end)
{
  
  for (int i=0; i<Ham.size_N(); i++){
    for (int j=0; j<Ham.size_Nd(); j++){
      dcomplex t = Ham(i,j);
      if (norm(t)<1e-10) t=0;
      if (fabs(t.imag())<1e-10) t = dcomplex(t.real(),0);
      if (fabs(t.real())<1e-10) t = dcomplex(0,t.imag());
      cout<<setw(10)<<t.real()<<" ";
    }
    cout<<endl;
  }
  
  
  for (int i=0; i<Ham.size_N(); i++)
    for (int j=0; j<Ham.size_Nd(); j++)
      if (norm(Ham(i,j)-Ham(j,i).conj())>1e-10) cerr<<"Not Hermitian "<<i<<" "<<j<<" "<<Ham(i,j)<<" "<<Ham(j,i)<<endl;

  for (list<pair<int,int> >::const_iterator l=start_end.begin(); l!=start_end.end(); l++){
    for (int i=l->first; i<l->second; i++){
      for (int j=0; j<Ham.size_Nd(); j++){
	if (j>=l->first && j<l->second) continue;
	if (norm(Ham(i,j))>1e-10) cerr<<"Off-diagonal elements in Hamiltonian! ERROR! "<<i<<" "<<j<<" "<<Ham(i,j)<<endl;
      }
    }
  }
}

string readableForm(const dcomplex& z)
{
  if (norm(z)<1e-6) return "0";
  stringstream ss;
  if (fabs(z.imag())<1e-6) {ss<<z.real(); return ss.str();}
  if (fabs(z.real())<1e-6) {ss<<z.imag()<<"*i"; return ss.str();}
  ss<<z.real()<<"+"<<z.imag()<<"*i";
  return ss.str();
}

void ComplexEigenvectors(function2D<dcomplex>& Teig, const vector<vector<double> >& KSE)
{
  deque<int> unresolved;
  for (int r=0; r<Teig.size_N(); r++){
    double isum=0, rsum=0;
    for (int p=0; p<Teig.size_Nd(); p++){
      isum += fabs(Teig(r,p).imag());
      rsum += fabs(Teig(r,p).real());
    }
    if (isum>1e-6){
      cerr<<"Eigenvector "<<r<<" is not purely real! Trying to resolve the issue.... "<<KSE[r][2]<<endl;
      if (rsum<1e-6){
 	for (int p=0; p<Teig.size_Nd(); p++) Teig(r,p) *= dcomplex(0,-1);
 	cerr<<"Problem solved!"<<endl;
 	continue;
      }
      cerr<<endl;
      unresolved.push_back(r);
    }
  }
  if (unresolved.size()>0){
    deque<deque<int> > unres;
    int ir=0;
    while (ir<unresolved.size()){
      double Ene = KSE[unresolved[ir]][2];
      deque<int> segment;
      while (ir<unresolved.size() && fabs(KSE[unresolved[ir]][2]-Ene)<1e-7){ segment.push_back(unresolved[ir]); ir++;}
      unres.push_back(segment);
    }
    for (int iq=0; iq<unres.size(); iq++){
      
      if (unres[iq].size()==2){// take difference and average of both eigenvectors
	vector<int> ri(2);
	ri[0] = unres[iq][0];
	ri[1] = unres[iq][1];
	function2D<dcomplex> V(2,Teig.size_Nd());
	function1D<double> sumr(2), sumi(2);
	sumr=0; sumi=0;
	for (int p=0; p<Teig.size_Nd(); p++){
	  V(0,p) = (Teig(ri[0],p)-Teig(ri[1],p))/sqrt(2);// difference
	  V(1,p) = (Teig(ri[0],p)+Teig(ri[1],p))/sqrt(2);// average
	  for (int j=0; j<2; j++){
	    sumr[j] += fabs(V(j,p).real());// check if they are now purely real or purely imaginary
	    sumi[j] += fabs(V(j,p).imag());
	  }
	}
	for (int j=0; j<2; j++){
	  if (fabs(sumr[j])<1e-6){// if vector is purely imaginary, make it real
	    for (int p=0; p<V.size_Nd(); p++) V(j,p) *= dcomplex(0,-1);
	    sumr[j] = sumi[j]; sumi[j]=0;
	  }
	}
	if (sumi[0]<1e-6 && sumi[1]<1e-6){// worked out, everything is fine
	  dcomplex prod=0;// need to orthogonalize
	  for (int p=0; p<V.size_Nd(); p++) prod += V(0,p)*V(1,p);
	  if (norm(prod)>1e-6){
	    for (int p=0; p<V.size_Nd(); p++) V(1,p) -= V(0,p)*prod;
	    double nrm=0;
	    for (int p=0; p<V.size_Nd(); p++) nrm += norm(V(1,p));
	    nrm = 1/sqrt(nrm);
	    for (int p=0; p<V.size_Nd(); p++) V(1,p) *= nrm;
	  }
	  cerr<<"Problem solved for "<<ri[0]<<" and "<<ri[1]<<endl;
	  for (int p=0; p<V.size_Nd(); p++){
	    Teig(ri[0],p) = V(0,p);
	    Teig(ri[1],p) = V(1,p);
	  }
	  for (int p=0; p<V.size_Nd(); p++)  cerr<<setw(12)<<readableForm(Teig(ri[0],p))<<", ";
	  cerr<<endl;
	  for (int p=0; p<V.size_Nd(); p++)  cerr<<setw(12)<<readableForm(Teig(ri[1],p))<<", ";
	  cerr<<endl;
	  continue;
	}
      }
      
      cerr<<"ERROR: Eigenvectors are not purely real nor purely imaginary! You need to take care of this case .... not yet supported! "<<endl;
      cerr<<"segment"<<iq<<" "<<endl;
      for (int ir=0; ir<unres[iq].size(); ir++){
	int r = unres[iq][ir];
	for (int p=0; p<Teig.size_Nd(); p++) cerr<<setw(12)<<readableForm(Teig(r,p))<<", ";
	cerr<<endl;
      }
    }
  }
}

class Ecmp{
  const vector<vector<double> >& KSE;
public:
  Ecmp(const vector<vector<double> >& KSE_) : KSE(KSE_){};
  bool operator()(int a, int b)
  {
    if (KSE[a][0]!=KSE[b][0]) return KSE[a][0]<KSE[b][0];
    return KSE[a][2]<KSE[b][2];
  }
};

class Eigenbase{
  int Nm, n;
  double Sz;
public:
  double t0, tp0, U;
private:
  operateLS op;
  PrimitiveShift sh;
  vector<int> base;
  map<int,int> ind;
  double maxSz_;
  vector<int> StateNumber;
  vector<int> basen;
  map<int,int> indn;
  // for Sz=Sz_min
  function2D<dcomplex> Teig;
  function2D<dcomplex> Teigp;
public:
  vector<int> parent;
  vector<vector<double> > KSE;
  vector<vector<int> > KSi;
public:
  // for the particular Sz
  function2D<dcomplex> Tn;
  function2D<dcomplex> Tnp;
  vector<vector<double> > KSEn;
  vector<vector<int> > KSin;
  int max_StateNumber;
  int maxStNum;
  vector<int> degeneracy;
  vector<double> energy;
  vector<double> spin;
  vector<int> momentum;
  vector<int> gamma;
  int start_dg;
public:
  Eigenbase(int Nm_, int n_, double t0_, double tp0_, double U_)
    : Nm(Nm_), n(n_), Sz(0.5*(n%2)), t0(t0_), tp0(tp0_), op(Nm), sh(Nm), U(U_) {}
  Eigenbase() : Nm(0) {};
  void Init(int Nm_, int n_, double t0_, double tp0_, double U_)
  {
    Nm = Nm_;
    n = n_;
    Sz = 0.5*(n%2);
    t0 = t0_; tp0 = tp0_; U=U_;
    op.Init(Nm);
    sh.Init(Nm);
  }
  void CreateDirectBase();
  void ConstructEigenbase();
  bool ConstructEigenbase_Sz(double Szn);
  int size0() const {return KSE.size();}
  int size() const {return KSEn.size();}
  double maxSz() const {return maxSz_;}
  int N() const {return Nm;}
  int n_() const {return n;}
  void NumberStates(int start);
  int StNumber(int i){return StateNumber[parent[i]];}
  int Degeneracy(int i){return degeneracy[StateNumber[parent[i]]-start_dg];}
  friend void ComputeFp(int K, int ud, Eigenbase& Eb1, Eigenbase& Eb2, function2D<dcomplex>& Fp0);
  friend void ComputeFp_debug(int K, int ud, Eigenbase& Eb1, Eigenbase& Eb2, function2D<dcomplex>& Fp0, function2D<dcomplex>& FpFm);
  friend void ComputeGz(int K, Eigenbase& Eb, function2D<dcomplex>& Gz0);
  friend void ComputeNk(int K, int updown, Eigenbase& Eb, function2D<dcomplex>& Nk0);
  friend void ComputeGp(int K, Eigenbase& Eb1, Eigenbase& Eb2, function2D<dcomplex>& Gp0);
  friend void ComputeHop(int K, Eigenbase& Eb, function2D<dcomplex>& hop);
};

void Eigenbase::CreateDirectBase()
{
  if (Nm==0){cerr<<"Did not initialize Eigenbase!"<<endl; exit(1);}
  CreateBase_n_Sz(Nm,n,Sz,op,base);
  if (base.size()==0){ cerr<<"Base is empty!"<<endl; return;}
  for (int i=0; i<base.size(); i++) ind[base[i]]=i;
  parent.resize(base.size());
  for (int i=0; i<base.size(); i++) parent[i]=i;
}
void Eigenbase::ConstructEigenbase()
{
  // First, find base (from direct base) with good quantum number K
  function2D<dcomplex> Tk(base.size(),base.size());
  Tk=0;
  function1D<int> NumPsis(Nm); NumPsis=0;
  function1D<dcomplex> psi(base.size());
  int l=0;
  for (int k=0; k<Nm; k++){// over all cluster K points
    map<int,int> database;
    int newstate, sgn, sgn1;
    for (int i=0; i<base.size(); i++){
      if (database.find(i)==database.end()){// this state did not yet occur
	if (l>base.size()) cerr<<"Something wrong! Too many eigenvectors! "<<endl;
	psi=0;
	for (int s=0; s<Nm; s++){           // creates |i>+exp(ikx)shx|i>+exp(iky)shy|i>+...
	  op.Shifts(base[i], sh.shxy[s].first, sh.shx, sh.shxy[s].second, sh.shy, newstate, sgn);
	  int j = ind[newstate];
	  psi[j] += sgn*sh.expiK[k][s];
	  database[j]++;
	}
	// normalization
	double nrm=0;
	for (int j=0; j<base.size(); j++) nrm += norm(psi[j]);
	// This new vector psi is nonzero (if the state is more symmetric, it can be zero) and needs to be stored
	if (nrm>1e-6){
	  nrm = 1/sqrt(nrm);
	  for (int j=0; j<base.size(); j++) Tk(l,j) = psi[j]*nrm;
	  l++;
	  NumPsis[k]++;
	}
      }
    }
  }
  if (l<base.size()) cerr<<"Number of transformation vectors is too small!"<<endl;

  // Conjugate Transform
  function2D<dcomplex> Tkp(base.size(),base.size());
  for (int i=0; i<base.size(); i++)
    for (int j=0; j<base.size(); j++)
      Tkp(i,j) = Tk(j,i).conj();

  
  //Next, find base with good quantum number K and S^2
  function2D<dcomplex> S2(base.size(), base.size());
  ComputeS2(base, ind, op, S2);
    
  function2D<dcomplex> temp(base.size(),base.size());
  
  temp.MProduct(Tk, S2);
  S2.MProduct(temp,Tkp);

  function2D<dcomplex> TKS(base.size(),base.size());
  KSE.resize(base.size());
  KSi.resize(base.size());
    
  function1D<dcomplex> work(5*base.size());
  function1D<double> rwork(5*base.size());
  int il=0;
  for (int k=0; k<Nm; k++){
    int size = NumPsis[k];
    if (size==0) continue;
    int offset = il;
    il += NumPsis[k];
    // Takes out one block of S2
    function2D<dcomplex> Snew(size,size);
    for (int i=0; i<size; i++)
      for (int j=0; j<size; j++) Snew(i,j) = S2(offset+i,offset+j);
    // Solves for eigensystem of a block of S2
    function1D<double> Seig(size);
    xsyev(size, Snew.MemPt(), Snew.fullsize_Nd(), Seig.MemPt(), work.MemPt(), work.size(), rwork.MemPt());
    // Constructs eigenvectors in terms of direct-base functions.
    // Also saves values of K and S for each eigenvector
    int kappa=0, S_last=-1;// index that distinguishes states with the same n,Sz,K,S
    for (int r=0; r<size; r++){
      for (int p=0; p<Tk.size_Nd(); p++){
	dcomplex sum=0;
	for (int j=0; j<size; j++) sum += Snew(r,j)*Tk(j+offset,p);
	TKS(r+offset,p) = sum;
      }
      // S for this eigenvector
      double St = 0.5*round(sqrt(1+4*Seig[r])-1);
      vector<double> tt(3);
      tt[0] = k;
      tt[1] = St;
      KSE[r+offset] = tt;
      vector<int> lsk(3);
      lsk[0] = k;
      lsk[1] = static_cast<int>(round(St-0.25));
      if (static_cast<int>(round(St-0.25))==S_last) kappa++;
      else{kappa=0; S_last=static_cast<int>(round(St-0.25));}
      lsk[2] = kappa;
      KSi[r+offset] = lsk;
    }
  }
  // Finds which blocks of H are nonzero
  list<pair<int,int> > start_end;
  l=0;
  for (int i=1; i<KSi.size(); i++){
    if (KSi[i][2]>KSi[i-1][2]) l++;
    else{
      start_end.push_back(make_pair(i-l-1,i));
      l=0;
    }
  }
  start_end.push_back(make_pair(KSi.size()-l-1,KSi.size()));

  function2D<dcomplex> TKSp(base.size(),base.size());
  for (int i=0; i<TKS.size_N(); i++)
    for (int j=0; j<TKS.size_Nd(); j++)
      TKSp(i,j) = TKS(j,i).conj();


  // Here, we calculate Hamiltonian in direct base
  function2D<dcomplex> Ham(base.size(), base.size());
  ComputeHamiltonian(t0, tp0, U, base, ind, op, sh, Ham);

  temp.MProduct(TKS, Ham);
  Ham.MProduct(temp,TKSp);
  Test(Nm, Ham, start_end);
    
  // Finally, calculating eigenstates within subblock of good quantum number K and S^2
  Teig.resize(base.size(),base.size());
  for (list<pair<int,int> >::const_iterator l=start_end.begin(); l!=start_end.end(); l++){
    int size = l->second-l->first;
    if (size==0) continue;
    int offset = l->first;
    // Takes out one block of H
    function2D<dcomplex> Hnew(size,size);
    for (int i=0; i<size; i++)
      for (int j=0; j<size; j++) Hnew(i,j) = Ham(offset+i,offset+j);
    // Solves for eigensystem of a block of H
    function1D<double> Heig(size);
    xsyev(size, Hnew.MemPt(), Hnew.fullsize_Nd(), Heig.MemPt(), work.MemPt(), work.size(), rwork.MemPt());

    // Constructs eigenvectors in terms of direct-base functions.
    // Also saves values of K and S for each eigenvector
    for (int r=0; r<size; r++){
      for (int p=0; p<TKS.size_Nd(); p++){
	dcomplex sum=0;
	for (int j=0; j<size; j++) sum += Hnew(r,j)*TKS(j+offset,p);
	Teig(r+offset,p) = sum;
      }
      if (fabs(Heig[r])<1e-10) Heig[r]=0;
      KSE[r+offset][2] = Heig[r];
    }
  }

  // What if eigenvectors are complex?
  ComplexEigenvectors(Teig,KSE);
  
  // Want to sort eigenvectors differently?
  if (true){
    vector<int> ind(KSE.size());
    for (int i=0; i<ind.size(); i++) ind[i]=i;
    Ecmp cmp(KSE);
    sort(ind.begin(),ind.end(),cmp);
    cout<<"New energies are:"<<endl;
    for (int i=0; i<KSE.size(); i++)
      cout<<setw(3)<<i<<" "<<setw(12)<<KSE[ind[i]][2]<<endl;
    vector<vector<double> > tKSE(KSE.size());
    function2D<dcomplex> tTeig(Teig.size_N(),Teig.size_Nd());
    for (int i=0; i<KSE.size(); i++){
      tKSE[i] = KSE[ind[i]];
      for (int j=0; j<Teig.size_Nd(); j++) tTeig[i][j] = Teig[ind[i]][j];
    }
    for (int i=0; i<KSE.size(); i++){
      KSE[i] = tKSE[i];
      for (int j=0; j<Teig.size_Nd(); j++) Teig[i][j] = tTeig[i][j];
    }
  }
  
  Teigp.resize(base.size(),base.size());
  for (int i=0; i<Teig.size_N(); i++)
    for (int j=0; j<Teig.size_Nd(); j++)
      Teigp(i,j) = Teig(j,i).conj();

 
  maxSz_=0;
  for (int i=0; i<KSE.size(); i++) if (KSE[i][1]>maxSz_) maxSz_=KSE[i][1];
}

bool Eigenbase::ConstructEigenbase_Sz(double Szn)
{
  int dS = static_cast<int>(Szn-Sz);
  if (dS==0){
    KSEn = KSE;
    KSin = KSi;
    Tn = Teig;
    Tnp = Teigp;
    basen = base;
    indn = ind;
    parent.resize(base.size());
    for (int i=0; i<basen.size(); i++) parent[i]=i;
    return true;
  }
  CreateBase_n_Sz(Nm,n,Szn,op,basen);
  if (basen.size()==0){ /*cerr<<"Base is empty!"<<endl;*/ return false;}
  for (int i=0; i<basen.size(); i++) indn[basen[i]]=i;
  // operatir S^+ or S^- is constructed
  function2D<dcomplex> Sp(base.size(),basen.size());
  ComputeSp(dS, base, ind, basen, indn, op, Sp);
  // temporary transformation matrix
  function2D<dcomplex> tTn(base.size(),basen.size());
  tTn.MProduct(Teig,Sp);
  
  KSEn.resize(basen.size());
  KSin.resize(basen.size());
  Tn.resize(basen.size(),basen.size());
  parent.resize(basen.size());
  
  int il=0;
  for (int i=0; i<base.size(); i++){
    if (il>basen.size()) {cerr<<"There should be no more than basen.size() eigenvectors!"<<endl; break;}
    double nrm=0;
    for (int j=0; j<basen.size(); j++) nrm += norm(tTn(i,j));
    if (nrm>1e-6) {
      nrm = 1/sqrt(nrm);
      for (int j=0; j<basen.size(); j++) Tn(il,j) = tTn(i,j)*nrm;
      KSEn[il] = KSE[i];
      KSin[il] = KSi[i];
      parent[il] = i;
      il++;
    }
  }
  if (il<basen.size()) cerr<<"There should be exactly basen.size() eigenvectors!"<<endl;
  // Conjugate of the transformation matrix
  Tnp.resize(basen.size(),basen.size());
  for (int i=0; i<Tn.size_N(); i++)
    for (int j=0; j<Tn.size_Nd(); j++)
      Tnp(i,j) = Tn(j,i).conj();

  return true;
}

void Eigenbase::NumberStates(int start)
{
  StateNumber.resize(base.size());
  for (int i=0; i<base.size(); i++) StateNumber[i]=-1;

  vector<double> St(base.size());
  
  int m = -1;
  int l = 0;
  StateNumber[l] = (++m)+start;
  St[m]=KSE[l][1];
  while(1){
#ifdef NO_SYMMETRY
#elif defined MOMENTUM_SYMMETRY
    for (int i=l+1; i<base.size(); i++){
      if (fabs(KSE[i][2]-KSE[l][2])<1e-8 && StateNumber[i]==-1 && KSi[i][1]==KSi[l][1] && KSi[i][2]==KSi[l][2]){
	if (sh.equivalent(KSi[i][0],KSi[l][0])) StateNumber[i] = m+start;// same S and E, equivalent K and same gamma
      }
    }
#else
    for (int i=l+1; i<base.size(); i++){
      if (fabs(KSE[i][2]-KSE[l][2])<1e-8 && StateNumber[i]==-1 && KSi[i][1]==KSi[l][1]){
	if (sh.equivalent(KSi[i][0],KSi[l][0])) StateNumber[i] = m+start;// same S and E, equivalent K and different gamma
      }
    }
#endif
    while(l<base.size() && StateNumber[++l]>0) ;
    if (l>=base.size()) break;
    StateNumber[l] = (++m)+start;
    St[m]=KSE[l][1];
  }
  max_StateNumber = m+start;
  maxStNum = m+1;
  
  start_dg = start;
  degeneracy.resize(m+1);
  energy.resize(m+1);
  spin.resize(m+1);
  momentum.resize(m+1);
  gamma.resize(m+1);
  for (int i=0; i<m+1; i++) degeneracy[i]=0;
  for (int i=0; i<base.size(); i++) degeneracy[StateNumber[i]-start]++;
  for (int i=0; i<m+1; i++) degeneracy[i] *= static_cast<int>(2*St[i]+1);
  for (int i=0; i<base.size(); i++) {
    energy[StateNumber[i]-start] = KSE[i][2];
    spin[StateNumber[i]-start] = KSE[i][1];
    momentum[StateNumber[i]-start] = sh.eqK(KSi[i][0]);
    gamma[StateNumber[i]-start] = KSi[i][2];
  }
  for (int i=base.size()-1; i>=0; i--)
    gamma[StateNumber[i]-start] = KSi[i][2];
}

void ComputeFp(int K, int ud, Eigenbase& Eb1, Eigenbase& Eb2, function2D<dcomplex>& Fp0)
{
  function2D<dcomplex> Fp(Eb1.size(),Eb2.size());
  Fp=0;
  list<pair<int,dcomplex> > sts;
  if (Eb1.N()!=Eb2.N()) cerr<<"Not equal N. Troubles!"<<endl;
  PrimitiveShift sh(Eb1.N());
  operateLS op(Eb1.N());
  for (int i=0; i<Eb1.size(); i++){
    int state = Eb1.basen[i];
    op.Fp(state, sts, sh, K, ud);
    for (list<pair<int,dcomplex> >::const_iterator l=sts.begin(); l!=sts.end(); l++){
      if (Eb2.indn.find(l->first)!=Eb2.indn.end()){// this state is part of the base
	int j = Eb2.indn[l->first];
	Fp(i,j) += l->second;
      }else{
	cerr<<"In Hubbard model should work!"<<endl;
      }
    }
  }
    
  // Transforming to eigenbases
  function2D<dcomplex> temp(Eb1.size(),Eb2.size());
  temp.MProduct(Eb1.Tn, Fp);
  Fp.MProduct(temp,Eb2.Tnp);
  for (int i=0; i<Eb1.size(); i++){
    for (int j=0; j<Eb2.size(); j++){
      //      if (fabs(Fp(i,j).imag())>1e-6) cerr<<"Fp is imaginary! Did not expect that! Was there a magnetic field? "<<Fp(i,j)<<endl;
      Fp0(j,i) = Fp(i,j);
    }
  }
}
void ComputeFp_debug(int K, int ud, Eigenbase& Eb1, Eigenbase& Eb2, function2D<dcomplex>& Fp0,
		     function2D<dcomplex>& FpFm)
{
  function2D<dcomplex> Fp(Eb1.size(),Eb2.size());
  Fp=0;
  list<pair<int,dcomplex> > sts;
  if (Eb1.N()!=Eb2.N()) cerr<<"Not equal N. Troubles!"<<endl;
  PrimitiveShift sh(Eb1.N());
  operateLS op(Eb1.N());
  for (int i=0; i<Eb1.size(); i++){
    int state = Eb1.basen[i];
    op.Fp(state, sts, sh, K, ud);
    for (list<pair<int,dcomplex> >::const_iterator l=sts.begin(); l!=sts.end(); l++){
      if (Eb2.indn.find(l->first)!=Eb2.indn.end()){// this state is part of the base
	int j = Eb2.indn[l->first];
	Fp(i,j) += l->second;
      }else{
	cerr<<"In Hubbard model should work!"<<endl;
      }
    }
  }


  // Start Transforming to eigenbases
  function2D<dcomplex> temp(Eb1.size(),Eb2.size());
  temp.MProduct(Eb1.Tn, Fp);
  Fp.MProduct(temp,Eb2.Tnp);
  // End Transforming to eigenbases
  
  
  // start debugging
  function2D<dcomplex> Fm(Eb2.size(),Eb1.size());
  for (int i=0; i<Fp.size_N(); i++)
    for (int j=0; j<Fp.size_Nd(); j++)
      Fm(j,i) = Fp(i,j).conj();
  
  function2D<dcomplex> tmp(Eb2.size(),Eb2.size());
  tmp.MProduct(Fm,Fp);
  cout<<"OFP IS"<<endl;
  for (int i=0; i<tmp.size_N(); i++){
    for (int j=0; j<tmp.size_Nd(); j++){
      if (norm(tmp[i][j])<1e-10) tmp[i][j]=0;
      cout<<setw(12)<<tmp[i][j].real()<<" ";
    }
    cout<<endl;
  }
  FpFm=tmp;
  // end debugging
    
  
  
  for (int i=0; i<Eb1.size(); i++){
    for (int j=0; j<Eb2.size(); j++){
      //      if (fabs(Fp(i,j).imag())>1e-6) cerr<<"Fp is imaginary! Did not expect that! Was there a magnetic field? "<<Fp(i,j)<<endl;
      Fp0(j,i) = Fp(i,j).conj();
    }
  }
}

void ComputeGz(int K, Eigenbase& Eb, function2D<dcomplex>& Gz0)
{
  function2D<dcomplex> Gz(Eb.size(),Eb.size());
  Gz=0;
  dcomplex fct;
  PrimitiveShift sh(Eb.N());
  operateLS op(Eb.N());
  for (int i=0; i<Eb.size(); i++){
    int state = Eb.basen[i];
    op.Sz(state, fct, sh, K);
    Gz(i,i) += fct;
  }
  // Transforming to eigenbases
  function2D<dcomplex> temp(Eb.size(),Eb.size());
  temp.MProduct(Eb.Tn, Gz);
  Gz.MProduct(temp,Eb.Tnp);
  for (int i=0; i<Eb.size(); i++){
    for (int j=0; j<Eb.size(); j++){
      Gz0(j,i) = Gz(i,j);
    }
  }
}

void ComputeNk(int K, int updown, Eigenbase& Eb, function2D<dcomplex>& Nk0)
{
  function2D<dcomplex> Nk(Eb.size(),Eb.size());
  Nk=0;
  dcomplex fct;
  PrimitiveShift sh(Eb.N());
  operateLS op(Eb.N());
  for (int i=0; i<Eb.size(); i++){
    int state = Eb.basen[i];
    op.Nks(state, updown, fct, sh, K);
    Nk(i,i) += fct;
  }
  // Transforming to eigenbases
  function2D<dcomplex> temp(Eb.size(),Eb.size());
  temp.MProduct(Eb.Tn, Nk);
  Nk.MProduct(temp,Eb.Tnp);
  for (int i=0; i<Eb.size(); i++){
    for (int j=0; j<Eb.size(); j++){
      Nk0(j,i) = Nk(i,j);
    }
  }
}

void ComputeHop(int K, Eigenbase& Eb, function2D<dcomplex>& hop)
{
  hop=0;
  operateLS op(Eb.N());
  PrimitiveShift sh(Eb.N());
  for (int i=0; i<Eb.size(); i++){
    map<int,dcomplex> sts;
    int state = Eb.basen[i];
    op.Hop(state, sh, sts, K);
    for (map<int,dcomplex>::const_iterator l=sts.begin(); l!=sts.end(); l++){
      if (Eb.indn.find(l->first)!=Eb.indn.end()) // otherwise is not part of the base
	hop(i,Eb.indn[l->first]) += l->second;
    }
  }
  // Transforming to eigenbases
  function2D<dcomplex> temp(Eb.size(),Eb.size());
  temp.MProduct(Eb.Tn, hop);
  hop.MProduct(temp,Eb.Tnp);
}

void ComputeGp(int K, Eigenbase& Eb1, Eigenbase& Eb2, function2D<dcomplex>& Gp0)
{
  function2D<dcomplex> Gp(Eb1.size(),Eb2.size());
  Gp=0;
  map<int,dcomplex> sts;
  if (Eb1.N()!=Eb2.N()) cerr<<"Not equal N. Troubles!"<<endl;
  PrimitiveShift sh(Eb1.N());
  operateLS op(Eb1.N());
  for (int i=0; i<Eb1.size(); i++){
    int state = Eb1.basen[i];
    op.Sp(state, sts, sh, K);
    for (map<int,dcomplex>::const_iterator l=sts.begin(); l!=sts.end(); l++){
      if (Eb2.indn.find(l->first)!=Eb2.indn.end()){// this state is part of the base
	int j = Eb2.indn[l->first];
	Gp(i,j) += l->second;
      }
    }
  }
  // Transforming to eigenbases
  function2D<dcomplex> temp(Eb1.size(),Eb2.size());
  temp.MProduct(Eb1.Tn, Gp);
  Gp.MProduct(temp,Eb2.Tnp);
  for (int i=0; i<Eb1.size(); i++){
    for (int j=0; j<Eb2.size(); j++){
      Gp0(j,i) = Gp(i,j);
    }
  }
}


void TestFp(int K, const Eigenbase& Eb1, const Eigenbase& Eb2, const function2D<dcomplex>& Fp)
{
  PrimitiveShift sh(Eb1.N());
  for (int i=0; i<Eb1.size(); i++){
    int K1 = static_cast<int>(Eb1.KSEn[i][0]);
    double S1 = Eb1.KSEn[i][1];
    for (int j=0; j<Eb2.size(); j++){
      int K2 = static_cast<int>(Eb2.KSEn[j][0]);
      double S2 = Eb2.KSEn[j][1];
      bool bS = (S2==S1+0.5) || (S2==S1-0.5);
      bool bK = sh.Plus(K1,K)==K2;
      if (bK && bS) continue;
      if (norm(Fp(j,i))>1e-10) cerr<<"TEZAVE!: "<<setw(2)<<i<<" "<<setw(2)<<j<<" "<<setw(11)<<Fp(j,i)<<endl;
    }
  }
}
void PrintFp(const Eigenbase& Eb1, const Eigenbase& Eb2, const function2D<dcomplex>& Fp)
{
  for (int i=0; i<Eb1.size(); i++){
    cout<<setw(2)<<i<<" "<<setw(3)<<Eb1.KSEn[i][0]<<" "<<setw(3)<<Eb1.KSEn[i][1]<<" "<<setw(3)<<Eb1.KSin[i][2]<<" "<<setw(11)<<Eb1.KSEn[i][2]<<endl;
  }
  cout<<endl;
  for (int i=0; i<Eb2.size(); i++){
    cout<<setw(2)<<i<<" "<<setw(3)<<Eb2.KSEn[i][0]<<" "<<setw(3)<<Eb2.KSEn[i][1]<<" "<<setw(3)<<Eb2.KSin[i][2]<<" "<<setw(11)<<Eb2.KSEn[i][2]<<endl;
  }
  cout<<endl;
  
  for (int j=0; j<Eb2.size(); j++){
    for (int i=0; i<Eb1.size(); i++){
      double r = Fp(j,i).real();
      if (fabs(r)<1e-6) r=0;
      if (fabs(Fp(j,i).imag())<1e-6)  cout<<setw(11)<<r<<" ";
      else {
	stringstream ss; ss<<r<<"+"<<Fp(j,i).imag()<<"i";
	cout<<setw(11)<<ss.str()<<" ";
      }
    }
    cout<<endl;
  }
}


class Ind_N_K_Sz{
  int Nm;
public:
  Ind_N_K_Sz(int Nm_) : Nm(Nm_){};
  int ind(int N, int K, double Sz)
  {
    return static_cast<int>(2*Sz+Nm) + (2*Nm+1)*K + (2*Nm+1)*Nm*N;
  }
  void inv(int ind, int& N, int& K, double& Sz)
  {
    int iSz = ind%(2*Nm+1);
    Sz = 0.5*(iSz-Nm);
    ind = ind/(2*Nm+1);
    K = ind%Nm;
    ind = ind/Nm;
    N = ind;
  }
  int max_size(){return (Nm+1)*Nm*(1+2*Nm);}
};
class State{
public:
  bool empty;
  int ii;
  int istart, iend;
  int N, K;
  double Sz;
  deque<double> Ene, Spin;
  State() : empty(true){};
  State(int ii_, int istart_, int iend_, deque<double>& Ene_, deque<double>& Spin_, int N_, int K_, double Sz_) :
    empty(false), ii(ii_), istart(istart_), iend(iend_), N(N_), K(K_), Sz(Sz_), Ene(Ene_), Spin(Spin_)  {};
  void Set(int istart_, int iend_, deque<double>& Ene_, deque<double>& Spin_, int N_, int K_, double Sz_)
  {
    empty=false;
    istart=istart_;
    iend=iend_;
    Ene = Ene_;
    Spin = Spin_;
    N = N_;
    K = K_;
    Sz = Sz_;
  }
  void Set_num(int ii_){ii=ii_;}
};

void Generate_normal(vector<Eigenbase>& Eb, const function2D<double>& eps_k_mq, const string& filename)
{
  int Nm = Eb[0].N();
  PrimitiveShift sh(Nm);
  Ind_N_K_Sz inde(2*Nm);
  
  Eb[0].ConstructEigenbase_Sz(0);

  vector<State> state(inde.max_size());
  {
    int ind = inde.ind(0,0,0.0);
    deque<double> Ene, Spin;
    Ene.push_back(0);
    Spin.push_back(0);
    state[ind].Set(0,1,Ene,Spin, 0, 0, 0.0);
  }
  
  int max_size=0;
  for (int n=1; n<=2*Nm; n++){
    cout<<"***************** n="<<n<<" ***********************"<<endl;
    for (double Sz=-Eb[n].maxSz(); Sz<=Eb[n].maxSz(); Sz++){
      Eb[n].ConstructEigenbase_Sz(Sz);
      cout<<"Sz="<<Sz<<endl;
      int i=0;
      while (i<Eb[n].size()){
	int kt = static_cast<int>(Eb[n].KSEn[i][0]);
	int ind = inde.ind(n,kt,Sz);
	int istart=i;
	deque<double> Ene, Spin;
	while (i<Eb[n].size() && Eb[n].KSEn[i][0]==kt){
	  Ene.push_back(Eb[n].KSEn[i][2]);
	  Spin.push_back(Eb[n].KSEn[i][1]);
	  i++;
	}
	int iend=i;
	state[ind].Set(istart,iend,Ene,Spin,n,kt,Sz);
	if (Ene.size()>max_size) max_size = Ene.size();
      }
    }
  }
  int ii=0;
  int ss=0;
  for (int i=0; i<state.size(); i++){
    if (!state[i].empty){
      state[i].Set_num(ii++);
      ss += state[i].Ene.size();
    }
  }
  int state_size=ii;
  vector<vector<vector<function2D<double> > > > Fp_matrix(state_size);
  vector<vector<function2D<double> > > M1_matrix(state_size);
  vector<vector<vector<int> > > Fp_index(state_size);
  for (int i=0; i<state_size; i++){
    Fp_matrix[i].resize(Nm);
    Fp_index[i].resize(Nm);
    M1_matrix[i].resize(Nm);
    for (int j=0; j<Nm; j++){
      Fp_matrix[i][j].resize(2);
      Fp_index[i][j].resize(2);
      for (int in=0; in<2; in++){
	Fp_index[i][j][in]=-1;
      }
    }
  }
  vector<int> equiv_k(Nm), deg_k(3);
  equiv_k[0]=0;
  equiv_k[1]=1;
  equiv_k[2]=1;
  equiv_k[3]=2;
  deg_k[0]=1;
  deg_k[1]=2;
  deg_k[2]=1;
  
  for (int n=1; n<=2*Nm; n++){
    cout<<"***************** n="<<n<<" ***********************"<<endl;
    for (double Sz=-Eb[n].maxSz(); Sz<=Eb[n].maxSz(); Sz++){
      Eb[n].ConstructEigenbase_Sz(Sz);
      
      for (int ud=-1; ud<=1; ud+=2){
	cout<<"************* n="<<n<<" Sz[n]="<<Sz<<" Sz[n-1]="<<Sz-0.5*ud<<"*****************"<<endl;
	int iud = (ud+1)/2;
	if (!Eb[n-1].ConstructEigenbase_Sz(Sz-0.5*ud)) continue;
	for (int K=0; K<Nm; K++){
	  function2D<dcomplex> Fp(Eb[n].size(),Eb[n-1].size());
	  ComputeFp(K, ud, Eb[n-1], Eb[n], Fp);
	  
	  for (int ik=0; ik<Nm; ik++){
	    int ii = inde.ind(n,ik,Sz);
	    if (state[ii].empty) continue;
	    int is = state[ii].istart;
	    int ie = state[ii].iend;
	    int ikq = sh.Plus(ik,sh.minus(K));
	    int iin = state[ii].ii;
	    int jj = inde.ind(n-1,ikq,Sz-0.5*ud);
	    if (state[jj].empty) continue;
	    int js = state[jj].istart;
	    int je = state[jj].iend;
	    int jjn = state[jj].ii;
	    Fp_index[jjn][K][iud]=iin;
	    Fp_matrix[jjn][K][iud].resize(je-js,ie-is);
	    for (int i=is; i<ie; i++){
	      for (int j=js; j<je; j++){
		double f = Fp[i][j].real();
		if (fabs(f)<1e-10) f=0;
		Fp_matrix[jjn][K][iud](j-js,i-is) = f;
	      }
	    }
	  }
	}
      }
    }
  }

  ofstream gout(filename.c_str()); gout.precision(12);
  gout<<"# Cix file for cluster DMFT with CTQMC"<<endl;
  gout<<"# cluster_size, number of states, number of baths, maximum_matrix_size "<<endl;
  gout<<Nm<<" "<<state_size<<" "<<2*Nm<<" "<<max_size<<endl;

  gout<<"# baths, dimension, global flip"<<endl;
  for (int ik=0; ik<equiv_k.size(); ik++){
    gout<<left<<setw(3)<<2*ik<<" "<<1<<" "<<setw(4)<<equiv_k[ik]<<" "<<setw(4)<<ik<<endl;
    gout<<left<<setw(3)<<2*ik+1<<" "<<1<<" "<<setw(4)<<equiv_k[ik]<<" "<<setw(4)<<ik<<endl;
  }
  gout<<right;
  gout<<"# cluster energies for non-equivalent baths, eps[k]"<<endl;
  gout<<eps_k_mq[0][0]<<" "<<eps_k_mq[0][1]<<" "<<eps_k_mq[0][3]<<endl;
 
  gout<<setw(2)<<"#"<<" "<<setw(3)<<"N"<<" "<<setw(3)<<"K"<<" "<<setw(4)<<"Sz"<<" "<<setw(3)<<"size"<<endl;
  for (int i=0; i<state.size(); i++){
    if (!state[i].empty){
      gout<<setw(2)<<state[i].ii+1<<" "<<setw(3)<<state[i].N<<" "<<setw(3)<<state[i].K<<" "<<setw(4)<<state[i].Sz<<" "<<setw(3)<<state[i].Ene.size()<<"    ";
      for (int ik=0; ik<Nm; ik++)
	for (int iud=0; iud<2; iud++)
	  gout<<setw(3)<<Fp_index[state[i].ii][ik][iud]+1<<" ";
      for (int l=0; l<state[i].Ene.size(); l++) gout<<setw(18)<<state[i].Ene[l]<<" ";
      for (int l=0; l<state[i].Ene.size(); l++) gout<<setw(18)<<state[i].Spin[l]<<" ";
      gout<<endl;
    }
  }
  gout<<"# matrix elements"<<endl;
  for (int i=0; i<state.size(); i++){
    if (!state[i].empty){
      for (int ik=0; ik<Nm; ik++)
	for (int iud=0; iud<2; iud++){
	  gout<<setw(2)<<state[i].ii+1<<" ";
	  gout<<setw(3)<<Fp_index[state[i].ii][ik][iud]+1<<"  ";
	  function2D<double>& fp = Fp_matrix[state[i].ii][ik][iud];
	  gout<<setw(3)<<fp.size_N()<<" "<<setw(3)<<fp.size_Nd()<<"  ";
	  for (int it=0; it<fp.size_N(); it++)
	    for (int jt=0; jt<fp.size_Nd(); jt++)
	      gout<<setw(20)<<fp(it,jt)<<" ";
	    
	  gout<<endl;


	  if (state[i].ii+1==63 || Fp_index[state[i].ii][ik][iud]+1==63){
	    cout<<"Fp: "<<state[i].ii+1<<" "<<Fp_index[state[i].ii][ik][iud]+1<<endl;
	    for (int jt=0; jt<fp.size_Nd(); jt++){
	      for (int it=0; it<fp.size_N(); it++)
		cout<<setw(20)<<fp(it,jt)<<" ";
	      cout<<endl;
	    }
	    cout<<endl;
	  }
	}
    }
  }
  gout<<"HB1"<<endl;
  gout<<"# number of operators needed"<<endl;
  gout<<0<<endl;
  /*
  gout<<"# high frequency expansion (Real and Imaginary part of Sigma)"<<endl;
  gout<<"$ReSigma =  "<<Eb[0].U<<"*$nf/2;"<<endl;
  gout<<"$ImSigma =  -"<<Eb[0].U<<"**2*$nf/2*(1-$nf/2)/$om[$im];"<<endl;
  gout<<"# number of operators needed"<<endl;
  gout<<"0"<<endl;
  */
}


void Generate_superc(vector<Eigenbase>& Eb, const function2D<double>& eps_k_mq, const string& filename)
{
  int Nm = Eb[0].N();
  PrimitiveShift sh(Nm);
  Ind_N_K_Sz inde(2*Nm);
  
  Eb[0].ConstructEigenbase_Sz(0);

  vector<State> state(inde.max_size());
  {
    int ind = inde.ind(0,0,0.0);
    deque<double> Ene, Spin;
    Ene.push_back(0);
    Spin.push_back(0);
    state[ind].Set(0,1,Ene,Spin, 0, 0, 0.0);
  }
  
  int max_size=0;
  for (int n=1; n<=2*Nm; n++){
    cout<<"***************** n="<<n<<" ***********************"<<endl;
    for (double Sz=-Eb[n].maxSz(); Sz<=Eb[n].maxSz(); Sz++){
      Eb[n].ConstructEigenbase_Sz(Sz);
      cout<<"Sz="<<Sz<<endl;
      int i=0;
      while (i<Eb[n].size()){
	int kt = static_cast<int>(Eb[n].KSEn[i][0]);
	int ind = inde.ind(n,kt,Sz);
	int istart=i;
	deque<double> Ene, Spin;
	while (i<Eb[n].size() && Eb[n].KSEn[i][0]==kt){
	  Ene.push_back(Eb[n].KSEn[i][2]);
	  Spin.push_back(Eb[n].KSEn[i][1]);
	  i++;
	}
	int iend=i;
	state[ind].Set(istart,iend,Ene,Spin,n,kt,Sz);
	if (Ene.size()>max_size) max_size = Ene.size();
      }
    }
  }
  int ii=0;
  int ss=0;
  for (int i=0; i<state.size(); i++){
    if (!state[i].empty){
      state[i].Set_num(ii++);
      ss += state[i].Ene.size();
    }
  }
  int state_size=ii;
  vector<vector<vector<function2D<double> > > > Fp_matrix(state_size);
  vector<vector<function2D<double> > > M1_matrix(state_size);
  vector<vector<vector<int> > > Fp_index(state_size);
  for (int i=0; i<state_size; i++){
    Fp_matrix[i].resize(Nm);
    Fp_index[i].resize(Nm);
    M1_matrix[i].resize(Nm);
    for (int j=0; j<Nm; j++){
      Fp_matrix[i][j].resize(2);
      Fp_index[i][j].resize(2);
      for (int in=0; in<2; in++){
	Fp_index[i][j][in]=-1;
      }
    }
  }
  vector<int> equiv_k(Nm), deg_k(3);
  equiv_k[0]=0;
  equiv_k[1]=1;
  equiv_k[2]=1;
  equiv_k[3]=2;
  deg_k[0]=1;
  deg_k[1]=2;
  deg_k[2]=1;
  
  for (int n=1; n<=2*Nm; n++){
    cout<<"***************** n="<<n<<" ***********************"<<endl;
    for (double Sz=-Eb[n].maxSz(); Sz<=Eb[n].maxSz(); Sz++){
      Eb[n].ConstructEigenbase_Sz(Sz);
      for (int ud=-1; ud<=1; ud+=2){
	cout<<"************* n="<<n<<" Sz[n]="<<Sz<<" Sz[n-1]="<<Sz-0.5*ud<<"*****************"<<endl;
	int iud = (ud+1)/2;
	if (!Eb[n-1].ConstructEigenbase_Sz(Sz-0.5*ud)) continue;
	for (int K=0; K<Nm; K++){
	  function2D<dcomplex> Fp(Eb[n].size(),Eb[n-1].size());
	  ComputeFp(K, ud, Eb[n-1], Eb[n], Fp);
	  for (int ik=0; ik<Nm; ik++){
	    int ii = inde.ind(n,ik,Sz);
	    if (state[ii].empty) continue;
	    int is = state[ii].istart;
	    int ie = state[ii].iend;
	    int ikq = sh.Plus(ik,sh.minus(K));
	    int iin = state[ii].ii;
	    int jj = inde.ind(n-1,ikq,Sz-0.5*ud);
	    if (state[jj].empty) continue;
	    int js = state[jj].istart;
	    int je = state[jj].iend;
	    int jjn = state[jj].ii;
	    Fp_index[jjn][K][iud]=iin;
	    Fp_matrix[jjn][K][iud].resize(je-js,ie-is);
	    for (int i=is; i<ie; i++){
	      for (int j=js; j<je; j++){
		double f = Fp[i][j].real();
		if (fabs(f)<1e-10) f=0;
		Fp_matrix[jjn][K][iud](j-js,i-is) = f;
	      }
	    }
	  }
	}
      }
      for (int Q=0; Q<Nm; Q++){
	int ii = inde.ind(n,Q,Sz);
	int is=0, ie=0, iin=-1;
	if (!state[ii].empty){
	  is = state[ii].istart;
	  ie = state[ii].iend;
	  iin = state[ii].ii;
	}
	cout<<"states Q="<<Q<<" ("<<is<<" "<<ie<<") "<<iin<<endl;
      }
    }
  }

   vector<vector<vector<function2D<double> > > > Fm_matrix(state_size);
   vector<vector<vector<int> > > Fm_index(state_size);
   for (int i=0; i<state_size; i++){
     Fm_matrix[i].resize(Nm);
     Fm_index[i].resize(Nm);
     for (int ik=0; ik<Nm; ik++){
       Fm_matrix[i][ik].resize(2);
       Fm_index[i][ik].resize(2);
       for (int iud=0; iud<2; iud++){
	 Fm_index[i][ik][iud]=-1;
       }
     }
   }
   for (int i=0; i<state_size; i++){
     for (int ik=0; ik<Nm; ik++){
       for (int iud=0; iud<2; iud++){
	 int j = Fp_index[i][ik][iud];
	 if (j<0) continue;
	 Fm_index[j][ik][iud]=i;
	 function2D<double>& fp = Fp_matrix[i][ik][iud];
	 Fm_matrix[j][ik][iud].resize(fp.size_Nd(),fp.size_N());
	 for (int l=0; l<fp.size_Nd(); l++)
	   for (int m=0; m<fp.size_N(); m++)
	     Fm_matrix[j][ik][iud](l,m) = fp(m,l);
       }
     }
   }
   
  vector<vector<string> > baths(2*Nm-2);
  baths[0].push_back("0");
  baths[1].push_back("-0*");
  baths[2].push_back("1");
  baths[2].push_back("3");
  baths[2].push_back("3");
  baths[2].push_back("-1*");
  baths[3].push_back("1");
  baths[3].push_back("-3");
  baths[3].push_back("-3");
  baths[3].push_back("-1*");
  baths[4].push_back("2");
  baths[5].push_back("-2*");
  vector<int> bath_sym(2*Nm);
  bath_sym[0]=0;
  bath_sym[1]=0;
  bath_sym[2]=1;
  bath_sym[3]=1;
  bath_sym[4]=2;
  bath_sym[5]=2;
  
  
  ofstream gout(filename.c_str()); gout.precision(12);
  gout<<"# Cix file for cluster DMFT with CTQMC"<<endl;
  gout<<"# cluster_size, number of states, number of baths, maximum_matrix_size "<<endl;
  gout<<Nm<<" "<<state_size<<" "<<2*Nm-2<<" "<<max_size<<endl;
  gout<<"# baths, dimension, symmetry, global flip"<<endl;
  for (int ik=0; ik<baths.size(); ik++){
    gout<<setw(3)<<ik<<" "<<setw(3)<<sqrt(baths[ik].size())<<"    ";
    for (int j=0; j<baths[ik].size(); j++) gout<<setw(4)<<baths[ik][j]<<" ";
    for (int j=baths[ik].size(); j<4; j++) gout<<setw(4)<<" "<<" ";
    gout<<setw(4)<<ik<<endl;
  }

  gout<<"# cluster energies for non-equivalent baths, eps[k]"<<endl;
  gout<<eps_k_mq[0][0]<<" "<<eps_k_mq[1][0]<<" "<<eps_k_mq[3][0]<<"  0"<<endl;  
  
  gout<<setw(2)<<"#"<<" "<<setw(3)<<"N"<<" "<<setw(3)<<"K"<<" "<<setw(4)<<"Sz"<<" "<<setw(3)<<"size"<<endl;
  for (int i=0; i<state.size(); i++){
    if (!state[i].empty){
      gout<<setw(2)<<state[i].ii+1<<" "<<setw(3)<<state[i].N<<" "<<setw(3)<<state[i].K<<" "<<setw(4)<<state[i].Sz<<" "<<setw(3)<<state[i].Ene.size()<<"    ";
      for (int ik=0; ik<Nm; ik++){
	gout<<setw(3)<<Fp_index[state[i].ii][ik][1]+1<<" ";
	gout<<setw(3)<<Fm_index[state[i].ii][ik][0]+1<<" ";
      }
      for (int l=0; l<state[i].Ene.size(); l++) gout<<setw(18)<<state[i].Ene[l]<<" ";
      for (int l=0; l<state[i].Ene.size(); l++) gout<<setw(18)<<state[i].Spin[l]<<" ";
      gout<<endl;
    }
  }
  gout<<"# matrix elements"<<endl;
  for (int i=0; i<state.size(); i++){
    if (!state[i].empty){
      for (int ik=0; ik<Nm; ik++)
	{
	  int iud=1;
	  gout<<setw(2)<<state[i].ii+1<<" ";
	  gout<<setw(3)<<Fp_index[state[i].ii][ik][iud]+1<<"  ";
	  function2D<double>& fp = Fp_matrix[state[i].ii][ik][iud];
	  gout<<setw(3)<<fp.size_N()<<" "<<setw(3)<<fp.size_Nd()<<"  ";
	  for (int it=0; it<fp.size_N(); it++)
	    for (int jt=0; jt<fp.size_Nd(); jt++)
	      gout<<setw(20)<<fp(it,jt)<<" ";
	  gout<<endl;
	  iud=0;
	  gout<<setw(2)<<state[i].ii+1<<" ";
	  gout<<setw(3)<<Fm_index[state[i].ii][ik][iud]+1<<"  ";
	  function2D<double>& fm = Fm_matrix[state[i].ii][ik][iud];
	  gout<<setw(3)<<fm.size_N()<<" "<<setw(3)<<fm.size_Nd()<<"  ";
	  for (int it=0; it<fm.size_N(); it++)
	    for (int jt=0; jt<fm.size_Nd(); jt++)
	      gout<<setw(20)<<fm(it,jt)<<" ";
	  gout<<endl;
	}
    }
  }
  gout<<"# high frequency expansion (Real and Imaginary part of Sigma)"<<endl;
  gout<<"$ReSigma =  ($k<3) ? "<<Eb[0].U<<"*$nf/2 : 0;"<<endl;
  gout<<"$ImSigma =  ($k<3) ? -"<<Eb[0].U<<"**2*$nf/2*(1-$nf/2)/$om[$im] : 0;"<<endl;
  gout<<"# number of operators needed"<<endl;
  gout<<"0"<<endl;

}


void resize(vector<vector<vector<vector<dcomplex> > > >& V, int spn)
{
  V.resize(spn);
  for (int i=0; i<spn; i++){
    V[i].resize(spn);
    for (int j=0; j<spn; j++){
      V[i][j].resize(spn);
      for (int l=0; l<spn; l++){
	V[i][j][l].resize(spn);
      }
    }
  }
}

void RenumberPseudoParticles(const function2D<bool>& sgm_exist, map<int,int>& off, vector<pair<int,int> >& pp0, list<list<int> >& pp)
{
  map<int,int> ioff;
  int spn = sgm_exist.size_N();
  for (int i=0; i<spn; i++){
    off[i*spn+i] = i;
    ioff[i] = i*spn+i;
  }
  int l=spn-1;
  for (int i=0; i<spn; i++){
    for (int j=i+1; j<spn; j++){
      if (!sgm_exist[i][j]) continue;
      off[i*spn+j] = ++l;
      ioff[l] = i*spn+j;
      off[j*spn+i] = l;
    }
  }
  list<int> pi;
  vector<int> rr(spn);
  for (int i=0; i<spn; i++) rr[i]=-1;
  int m=-1;
  l=0;
  rr[l] = ++m;
  pi.push_back(0);
  while(1){
    for (int i=l+1; i<spn; i++)
      if (sgm_exist[l][i] && rr[i]<0) {rr[i]=m; pi.push_back(i);}
    pp.push_back(pi); pi.clear();
    while(l<spn && rr[++l]>=0);
    if (l>=spn) break;
    rr[l] = ++m;
    pi.push_back(l);
  }

  pp0.resize(ioff.size());
  for (map<int,int>::const_iterator l=ioff.begin(); l!=ioff.end(); l++) pp0[l->first] = make_pair((l->second)/spn,(l->second)%spn);
  
  cout<<" rr = "<<endl;
  for (int i=0; i<spn; i++) cout<<setw(3)<<i<<setw(5)<<rr[i]<<endl;
  
  cout<<"Pseudoparticles:"<<endl;
  for (int i=0; i<pp0.size(); i++)
    cout<<setw(3)<<i<<"      "<<setw(3)<<pp0[i].first<<"  "<<setw(3)<<pp0[i].second<<endl;
  
  cout<<"pairs : "<<endl;
  for (list<list<int> >::const_iterator i=pp.begin(); i!=pp.end(); i++){
    for (list<int>::const_iterator j=i->begin(); j!=i->end(); j++){
      for (list<int>::const_iterator k=j; k!=i->end(); k++){
	cout<<setw(3)<<off[(*j)*spn+(*k)]<<"  ";
      }
    }
    cout<<endl;
  }
}


int main(int argc, char *argv[])
{// Usgae: hub -U 6.0 -s 1
  double tp=-0.;
  double U = 6.;
  int Nm=4;
  // N=4 -> spn=26
  // N=8 -> spn=673
  int superc=0;
  
  stringstream inp;
  for (int i=0; i<argc/2; i++){
    inp<<argv[2*i+1]<<" "<<argv[2*i+2]<<endl;
  }
  while(inp){
    string str;
    inp>>str;
    if (str[0]=='#') inp.ignore(2000,'\n');
    if (str=="-U") inp>>U;
    if (str=="-s") inp>>superc;
    inp.ignore(1000,'\n');
  }
  
  double t0, tp0;
  if (Nm==4) {t0=0.5; tp0=0.25*tp;}//{t0=2/M_PI; tp0=4*tp/(M_PI*M_PI);}
  if (Nm==8) {t0=8/(M_PI*M_PI);}
  
  vector<Eigenbase> Eb(2*Nm+1);

 
  int maxStNum=0;
  for (int n=0; n<=2*Nm; n++){
    cout<<"********** n="<<n<<endl;
    Eb[n].Init(Nm, n, t0, tp0, U);
    Eb[n].CreateDirectBase();
    Eb[n].ConstructEigenbase();
    Eb[n].NumberStates(maxStNum);
    maxStNum += Eb[n].maxStNum;
  }
  int spn = maxStNum;

  vector<int> degeneracy(spn);
  degeneracy[0]=1;
  for (int n=1; n<=Nm; n++)
    for (int i=0; i<Eb[n].max_StateNumber-Eb[n-1].max_StateNumber; i++)
      degeneracy[i+Eb[n-1].max_StateNumber+1] = Eb[n].degeneracy[i];

  vector<pair<int,int> > state_index(Eb[Nm].max_StateNumber+1);
  for (int n=0; n<=Nm; n++)
    for (int i=0; i<Eb[n].maxStNum; i++){
      state_index[i+Eb[n].start_dg] = make_pair(n,i);
    }

  cout<<"Starting (original) states: "<<endl;
  cout<<setw(2)<<"#"<<" "<<setw(3)<<"n"<<" "<<setw(3)<<"deg"<<" "<<setw(11)<<"energy"<<" "<<setw(5)<<"S"<<" "<<setw(3)<<"K"<<" "<<setw(3)<<"g"<<endl;
  for (int j=0; j<state_index.size(); j++){
    int n = state_index[j].first;
    int i = state_index[j].second;
    cout<<setw(2)<<j<<" "<<setw(3)<<n<<" "<<setw(3)<<Eb[n].degeneracy[i]<<" "<<setw(11)<<Eb[n].energy[i]<<" "<<setw(5)<<Eb[n].spin[i]<<" "<<setw(3)<<Eb[n].momentum[i]<<" "<<setw(3)<<Eb[n].gamma[i]<<endl;
  }

  function2D<double> eps_k_mq(Nm,Nm);
  PrimitiveShift sh(Nm);
  for (int i=0; i<Nm; i++){
    for (int j=0; j<Nm; j++){
      eps_k_mq[i][j] = -2*t0*(cos(sh.K[i].first-sh.K[j].first)+cos(sh.K[i].second-sh.K[j].second))
	-4*tp0*cos(sh.K[i].first-sh.K[j].first)*cos(sh.K[i].second-sh.K[j].second);
    }
  }

  stringstream fname;
  if (!superc){
    fname<<"hubbard_U_"<<U<<"_normal.cix";
    Generate_normal(Eb, eps_k_mq, fname.str());
  }else{
    fname<<"hubbard_U_"<<U<<"_superc.cix";
    Generate_superc(Eb, eps_k_mq, fname.str());
  }
  
  return 0;
}
  
