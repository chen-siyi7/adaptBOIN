#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>
#include <utility>

using namespace Rcpp;

// ------------------------------------------------------------
// SPECIAL FUNCTIONS
// Regularised incomplete beta via Lentz continued fraction.
// ------------------------------------------------------------

static double betacf(double a, double b, double x) {
  const int    MAXIT = 200;
  const double EPS = 3e-10, FPMIN = 1e-30;
  double qab = a+b, qap = a+1.0, qam = a-1.0;
  double c = 1.0, d = 1.0 - qab*x/qap;
  if (std::abs(d) < FPMIN) d = FPMIN;
  d = 1.0/d; double h = d;
  for (int m = 1; m <= MAXIT; ++m) {
    double m2 = 2*m;
    double aa = m*(b-m)*x/((qam+m2)*(a+m2));
    d = 1.0+aa*d; if (std::abs(d)<FPMIN) d=FPMIN;
    c = 1.0+aa/c; if (std::abs(c)<FPMIN) c=FPMIN;
    d = 1.0/d; h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d = 1.0+aa*d; if (std::abs(d)<FPMIN) d=FPMIN;
    c = 1.0+aa/c; if (std::abs(c)<FPMIN) c=FPMIN;
    d = 1.0/d; double del = d*c; h *= del;
    if (std::abs(del-1.0) < EPS) break;
  }
  return h;
}

static double pbeta_cpp(double x, double a, double b) {
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0;
  double lbeta = lgamma(a)+lgamma(b)-lgamma(a+b);
  double bt = std::exp(std::log(x)*a + std::log(1.0-x)*b - lbeta);
  if (x < (a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
  else                         return 1.0 - bt*betacf(b,a,1.0-x)/b;
}

static double log_choose(int n, int k) {
  if (k<0||k>n) return -1e300;
  return lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1);
}

static double binom_logpmf(int k, int n, double p) {
  if (p<=0.0) return (k==0)?0.0:-1e300;
  if (p>=1.0) return (k==n)?0.0:-1e300;
  return log_choose(n,k)+k*std::log(p)+(n-k)*std::log(1.0-p);
}

// ------------------------------------------------------------
// RANDOM NUMBER GENERATION  (per-trial mt19937)
// ------------------------------------------------------------

static double rgamma_mt(double shape, std::mt19937& rng) {
  if (shape < 1.0) {
    std::uniform_real_distribution<double> u(0.0,1.0);
    return rgamma_mt(1.0+shape,rng)*std::pow(u(rng),1.0/shape);
  }
  double d=shape-1.0/3.0, c=1.0/std::sqrt(9.0*d);
  std::normal_distribution<double>       norm(0.0,1.0);
  std::uniform_real_distribution<double> u(0.0,1.0);
  for (;;) {
    double x, v;
    do { x=norm(rng); v=1.0+c*x; } while(v<=0.0);
    v=v*v*v; double uu=u(rng);
    if (uu < 1.0-0.0331*(x*x)*(x*x)) return d*v;
    if (std::log(uu) < 0.5*x*x+d*(1.0-v+std::log(v))) return d*v;
  }
}

// Dirichlet(alpha,...,alpha) with K+2 components -> cumulative weight vector
// w[k] = sum(v[0..k]), length K+1, monotone increasing in [0,1]
static std::vector<double> rdirichlet_w(int K, double alpha, std::mt19937& rng) {
  int n = K+2;
  std::vector<double> g(n); double s=0.0;
  for (int i=0;i<n;++i){ g[i]=rgamma_mt(alpha,rng); s+=g[i]; }
  std::vector<double> w(K+1); double cs=0.0;
  for (int k=0;k<=K;++k){ cs+=g[k]/s; w[k]=cs; }
  return w;
}

static int rbinom_cpp(int n, double p, std::mt19937& rng) {
  if (n==0||p<=0.0) return 0;
  if (p>=1.0) return n;
  return std::binomial_distribution<int>(n,p)(rng);
}

// ------------------------------------------------------------
// BERNSTEIN BASIS & EVALUATION
// ------------------------------------------------------------

static double bern_basis(double d, int k, int K) {
  if (d<=0.0) return (k==0)?1.0:0.0;
  if (d>=1.0) return (k==K)?1.0:0.0;
  return std::exp(log_choose(K,k)+k*std::log(d)+(K-k)*std::log(1.0-d));
}

static double bern_eval(const std::vector<double>& w, double d, int K) {
  double v=0.0;
  for (int k=0;k<=K;++k) v+=w[k]*bern_basis(d,k,K);
  return v;
}

// ------------------------------------------------------------
// PAVA (equal-weight isotonic regression, in-place)
// ------------------------------------------------------------

static void pava(std::vector<double>& y) {
  int n=(int)y.size();
  std::vector<double> mu(y), wt(n,1.0);
  std::vector<int>    sz(n,1);
  int nb=n, i=0;
  while (i<nb-1) {
    if (mu[i]>mu[i+1]) {
      double nw=wt[i]+wt[i+1];
      mu[i]=(wt[i]*mu[i]+wt[i+1]*mu[i+1])/nw;
      wt[i]=nw; sz[i]+=sz[i+1];
      mu.erase(mu.begin()+i+1);
      wt.erase(wt.begin()+i+1);
      sz.erase(sz.begin()+i+1);
      --nb; if(i>0)--i;
    } else ++i;
  }
  int pos=0;
  for (int b=0;b<nb;++b)
    for (int j=0;j<sz[b];++j) y[pos++]=mu[b];
}

// ------------------------------------------------------------
// POSTERIOR STATE PROBABILITIES
// p1=P(pi<phi1), p2=P(phi1<=pi<=phi2), p3=P(pi>phi2)
// ------------------------------------------------------------

struct PP { double p1,p2,p3; };

static PP post_probs(int s, int n,
                     double phi1, double phi2,
                     double a0,   double b0) {
  double a=a0+s, b=b0+n-s;
  double p1=pbeta_cpp(phi1,a,b);
  double p3=1.0-pbeta_cpp(phi2,a,b);
  return {p1, std::max(0.0,1.0-p1-p3), p3};
}

// ------------------------------------------------------------
// LOSS MATRIX & BAYES ACTION
// L[action*3 + state], actions 0=E,1=S,2=D; states 0=H1(under),1=H2(target),
// 2=H3(over). The matrix is now a parameter (see adapt_params$loss). The
// symmetric default reproduces the original paper. An overdose-aversion
// multiplier k>1 on the H3-column penalties (escalate-when-toxic, stay-when-
// toxic) pushes both boundaries down; single-crossing (Assumption 3.1) still
// holds for k > 0.5 since ell_E3 > ell_E2 + ell_S3 remains satisfied.
// ------------------------------------------------------------

static const double LOSS_SYM[9] = { 0.0,0.5,1.5,  0.5,0.0,0.5,  1.5,0.5,0.0 };

static int bayes_action(int s, int n,
                         double phi1, double phi2,
                         double a0,   double b0,
                         const double* L) {
  auto pp = post_probs(s,n,phi1,phi2,a0,b0);
  double p[3]={pp.p1,pp.p2,pp.p3}, risk[3];
  for (int a=0;a<3;++a){
    risk[a]=0.0;
    for (int h=0;h<3;++h) risk[a]+=L[a*3+h]*p[h];
  }
  return (int)(std::min_element(risk,risk+3)-risk);
}

// ------------------------------------------------------------
// SAFETY RULE  P(pi > phi2 | s,n) > rho  -> eliminate
// ------------------------------------------------------------

// phi_elim: toxicity level used in the overdose test (standard BOIN uses
//   phi_tgt=0.25; the original submission used phi2=0.35).
// ea0,eb0: prior for the elimination posterior (standard BOIN uses Beta(1,1);
//   the original submission reused the Jeffreys 0.5,0.5).
static bool is_elim(int s, int n,
                    double phi_elim, double rho,
                    double ea0,      double eb0) {
  if (n==0) return false;
  return (1.0-pbeta_cpp(phi_elim,ea0+s,eb0+n-s)) > rho;
}

// ------------------------------------------------------------
// ISOTONIC MTD SELECTION
// Apply PAVA to empirical rates at visited non-eliminated doses;
// select dose whose smoothed rate is closest to phi_tgt.
// ------------------------------------------------------------

// tie_high: when true, exact ties in |smoothed rate - phi_tgt| are broken
//   toward the higher dose (manuscript convention); when false, toward the
//   lower dose (Referee 1 conservative convention). Applied consistently with
//   the true-MTD definition so PCS is computed against a matching convention.
static int iso_mtd(const std::vector<int>& nv,
                   const std::vector<int>& sv,
                   const std::vector<bool>& elim,
                   double phi_tgt, int J, bool tie_high) {
  std::vector<int> adm;
  for (int j=0;j<J;++j) if(nv[j]>0&&!elim[j]) adm.push_back(j);
  if (adm.empty()) return -1;
  std::vector<double> r(adm.size());
  for (int i=0;i<(int)adm.size();++i) r[i]=(double)sv[adm[i]]/nv[adm[i]];
  pava(r);
  int best=adm[0]; double bd=std::abs(r[0]-phi_tgt);
  for (int i=1;i<(int)adm.size();++i){
    double d=std::abs(r[i]-phi_tgt);
    if(tie_high ? (d<=bd) : (d<bd)){bd=d;best=adm[i];}
  }
  return best;
}

// ------------------------------------------------------------
// BERNSTEIN MTD  (IS with ESS fallback)
// ------------------------------------------------------------

struct BernRes { int mtd; double ess; bool fell_back; };

static BernRes bern_mtd(const std::vector<int>& nv,
                         const std::vector<int>& sv,
                         const std::vector<bool>& elim,
                         double phi_tgt,
                         int K, double alpha_d,
                         int M, double ESS_thr,
                         int J, bool tie_high, std::mt19937& rng) {
  std::vector<int> adm;
  for (int j=0;j<J;++j) if(nv[j]>0&&!elim[j]) adm.push_back(j);
  if (adm.empty()) return {-1,0.0,true};

  std::vector<double> ds(J);
  for (int j=0;j<J;++j) ds[j]=(double)j/(J-1);

  std::vector<double> lw(M);
  std::vector<std::vector<double>> Wmat(M);
  for (int i=0;i<M;++i){
    Wmat[i]=rdirichlet_w(K,alpha_d,rng);
    double ll=0.0;
    for (int j=0;j<J;++j) if(nv[j]>0){
      double ph=std::max(1e-10,std::min(1.0-1e-10,bern_eval(Wmat[i],ds[j],K)));
      ll+=binom_logpmf(sv[j],nv[j],ph);
    }
    lw[i]=ll;
  }

  double lmax=*std::max_element(lw.begin(),lw.end());
  std::vector<double> w(M); double ws=0.0;
  for (int i=0;i<M;++i){ w[i]=std::exp(lw[i]-lmax); ws+=w[i]; }
  for (int i=0;i<M;++i) w[i]/=ws;
  double ssq=0.0;
  for (int i=0;i<M;++i) ssq+=w[i]*w[i];
  double ess=1.0/ssq;

  if (ess<ESS_thr){
    return {iso_mtd(nv,sv,elim,phi_tgt,J,tie_high), ess, true};
  }

  std::vector<double> ppi(J,0.0);
  for (int j=0;j<J;++j)
    for (int i=0;i<M;++i)
      ppi[j]+=w[i]*bern_eval(Wmat[i],ds[j],K);

  int best=adm[0]; double bd=std::abs(ppi[adm[0]]-phi_tgt);
  for (int k=1;k<(int)adm.size();++k){
    double d=std::abs(ppi[adm[k]]-phi_tgt);
    if(tie_high ? (d<=bd) : (d<bd)){bd=d;best=adm[k];}
  }
  return {best,ess,false};
}

// ------------------------------------------------------------
// CRM  (one-parameter power model, grid posterior)
// pi_j(a) = skel_j ^ exp(a),  a ~ N(0, sd^2)
// ------------------------------------------------------------

struct CrmRes {
  int next_dose;
  std::vector<bool>   elim;
  std::vector<double> pi_mean;
};

static CrmRes crm_step(const std::vector<int>&    nv,
                        const std::vector<int>&    sv,
                        std::vector<bool>          elim,
                        const std::vector<double>& skel,
                        double phi_tgt, double elim_thr,
                        double sd, int ng, int J) {
  std::vector<double> ag(ng), lp(ng);
  double lo=-3.0*sd, hi=3.0*sd;
  for (int k=0;k<ng;++k) ag[k]=lo+(hi-lo)*k/(ng-1);

  for (int k=0;k<ng;++k){
    double ea=std::exp(ag[k]);
    lp[k]=-0.5*(ag[k]*ag[k])/(sd*sd);
    for (int j=0;j<J;++j) if(nv[j]>0){
      double p=std::max(1e-10,std::min(1.0-1e-10,std::pow(skel[j],ea)));
      lp[k]+=binom_logpmf(sv[j],nv[j],p);
    }
  }
  double lmax=*std::max_element(lp.begin(),lp.end());
  std::vector<double> post(ng); double ps=0.0;
  for (int k=0;k<ng;++k){ post[k]=std::exp(lp[k]-lmax); ps+=post[k]; }
  for (int k=0;k<ng;++k) post[k]/=ps;

  std::vector<double> pm(J,0.0);
  for (int j=0;j<J;++j)
    for (int k=0;k<ng;++k)
      pm[j]+=post[k]*std::pow(skel[j],std::exp(ag[k]));

  for (int j=0;j<J;++j)
    if(nv[j]>0&&pm[j]>elim_thr) elim[j]=true;

  int best=-1; double bd=1e9;
  for (int j=0;j<J;++j) if(!elim[j]){
    double d=std::abs(pm[j]-phi_tgt);
    if(d<bd){bd=d;best=j;}
  }
  return {best,elim,pm};
}

// ------------------------------------------------------------
// mTPI-2  (UPM criterion)
// ------------------------------------------------------------

static int mtpi2_action(int s, int n,
                         double phi1, double phi2,
                         double a0,   double b0) {
  double a=a0+s, b=b0+n-s;
  double m1=pbeta_cpp(phi1,a,b);
  double m2=pbeta_cpp(phi2,a,b)-m1;
  double m3=1.0-pbeta_cpp(phi2,a,b);
  double u1=m1/phi1, u2=m2/(phi2-phi1), u3=m3/(1.0-phi2);
  if(u1>=u2&&u1>=u3) return 0;
  if(u3>=u2&&u3>=u1) return 2;
  return 1;
}

// ------------------------------------------------------------
// gBOINS shrinkage boundaries (Mu, Hu, Xu, Pan 2021, BMC Med Res Methodol),
// binary endpoint. Sample-size-dependent (shrinking) escalation/de-escalation
// rate boundaries. Lead-in: fixed gBOIN (= BOIN) boundaries for n <= N0, then
// shrinkage boundaries. Validated against the paper's Table 1 for phi0=0.2,0.3.
// gamma_k = exp(c_k * n^eps); phi1* = argmax_{mu<phi0} f(mu,g1),
// phi2* = argmin_{mu>phi0} f(mu,g2); boundaries then follow the BOIN form.
// ------------------------------------------------------------

static std::pair<double,double> gboins_bounds(double phi0, int n,
                                              double c1, double c2,
                                              double eps, int N0) {
  double phi1 = 0.6*phi0, phi2 = 1.4*phi0;
  double le, ld;
  if (n <= N0) {
    le = std::log((1.0-phi1)/(1.0-phi0)) /
         std::log(phi0*(1.0-phi1)/((1.0-phi0)*phi1));
    ld = std::log((1.0-phi0)/(1.0-phi2)) /
         std::log(phi2*(1.0-phi0)/((1.0-phi2)*phi0));
    return {le, ld};
  }
  double g1 = std::exp(c1*std::pow((double)n,eps));
  double g2 = std::exp(c2*std::pow((double)n,eps));
  auto fval = [&](double mu, double g){
    double num = std::log(g) - n*(std::log(1.0-mu) - std::log(1.0-phi0));
    double den = std::log(mu/(1.0-mu)) - std::log(phi0/(1.0-phi0));
    return num/den;
  };
  const double step = 1e-4;
  double p1 = 0.6*phi0, bestv = -1e300;
  for (double mu=step; mu<phi0; mu+=step){ double v=fval(mu,g1); if(v>bestv){bestv=v;p1=mu;} }
  double p2 = 1.4*phi0, bestw = 1e300;
  for (double mu=phi0+step; mu<1.0; mu+=step){ double v=fval(mu,g2); if(v<bestw){bestw=v;p2=mu;} }
  le = std::log((1.0-p1)/(1.0-phi0)) /
       std::log(phi0*(1.0-p1)/((1.0-phi0)*p1));
  ld = std::log((1.0-phi0)/(1.0-p2)) /
       std::log(p2*(1.0-phi0)/((1.0-p2)*phi0));
  return {le, ld};
}

// ============================================================
// TRIAL SIMULATION
// design: 0=adaptive_iso 1=adaptive_bern 2=boin_bern
//         3=boin 4=crm 5=mtpi2
// ============================================================

struct TrialRes {
  int mtd_sel, true_mtd, pcs;
  double pod; int dlts;
  int over_sel, under_sel, any_above, n_above;
  double ess; int fell_back;
  std::vector<int> nv, sv;
  int mtd_elim;            // 1 if the TRUE MTD was eliminated during the trial
  std::vector<int> elimv;  // per-dose elimination flags at end of trial
};

static TrialRes sim_trial(
    const std::vector<double>& pi,
    int design,
    const List& atbl,          // adaptive boundary table (R list)
    double phi_tgt, double phi1, double phi2,
    double rho, double a0, double b0,
    int N, int m, int J,
    double lam1, double lam2,  // BOIN fixed boundaries
    int K, double alpha_d, int M_IS, double ESS_thr,
    const std::vector<double>& skel,
    double crm_sd, int crm_ng, double crm_thr,
    bool tie_high,                            // higher-dose tie-break if true
    double phi_elim, double elim_a0, double elim_b0, // elimination rule
    double gb_c1, double gb_c2, double gb_eps, int gb_N0, // gBOINS calibration
    uint32_t seed)
{
  std::mt19937 rng(seed);

  std::vector<int>  nv(J,0), sv(J,0);
  std::vector<bool> elim(J,false);
  int cur=0, enrolled=0;

  while (enrolled < N) {
    int cn = std::min(m, N-enrolled);
    int nd = rbinom_cpp(cn, pi[cur], rng);
    nv[cur]+=cn; sv[cur]+=nd; enrolled+=cn;
    if (enrolled>=N) break;

    // safety check on current dose
    if (is_elim(sv[cur],nv[cur],phi_elim,rho,elim_a0,elim_b0))
      for(int j=cur;j<J;++j) elim[j]=true;

    // check all doses at/below current eliminated
    bool all_e=true;
    for(int j=0;j<=cur;++j) if(!elim[j]){all_e=false;break;}
    if(all_e) break;

    // escalation decision
    int action=1;  // default stay
    int next=cur;

    if (design==4) {
      // CRM: recommend globally
      CrmRes cr=crm_step(nv,sv,elim,skel,phi_tgt,crm_thr,crm_sd,crm_ng,J);
      elim=cr.elim;
      if(cr.next_dose<0) break;
      next=cr.next_dose;
    } else {
      // derive action from the relevant rule
      if (design==0||design==1) {
        // adaptive table lookup
        std::string key=std::to_string(nv[cur]);
        if (atbl.containsElementNamed(key.c_str())){
          IntegerVector av=atbl[key];
          action=av[sv[cur]];
        } else {
          // Unreachable safety net: the adaptive table always covers every
          // cohort size n in {m,...,N}, so this branch is never taken in
          // practice. Uses the symmetric default only as a fallback.
          action=bayes_action(sv[cur],nv[cur],phi1,phi2,a0,b0,LOSS_SYM);
        }
      } else if (design==2||design==3) {
        double rate=(double)sv[cur]/nv[cur];
        if(rate<=lam1) action=0;
        else if(rate>=lam2) action=2;
        else action=1;
      } else if (design==6) {
        // gBOINS: sample-size-dependent shrinkage boundaries
        double rate=(double)sv[cur]/nv[cur];
        std::pair<double,double> gb=gboins_bounds(phi_tgt,nv[cur],gb_c1,gb_c2,gb_eps,gb_N0);
        if(rate<=gb.first) action=0;
        else if(rate>=gb.second) action=2;
        else action=1;
      } else { // mtpi2 (design==5)
        action=mtpi2_action(sv[cur],nv[cur],phi1,phi2,a0,b0);
      }

      // translate action to next dose
      if(action==0){           // escalate
        next=cur+1;
        if(next>=J||elim[next]) next=cur;
        // safety pre-check on next
        if(next>cur&&nv[next]>0&&
           is_elim(sv[next],nv[next],phi_elim,rho,elim_a0,elim_b0)){
          for(int j=next;j<J;++j) elim[j]=true;
          next=cur;
        }
      } else if(action==2) {   // de-escalate
        next=cur-1;
        if(next<0) next=0;
        while(next>0&&elim[next]) --next;
        if(elim[next]) break;
      }
    }
    cur=next;
  }

  // ---- end-of-trial MTD selection ----
  double ess=-1.0; int fell_back=0, mtd_sel=-1;

  if (design==0||design==3||design==5||design==6) {
    mtd_sel=iso_mtd(nv,sv,elim,phi_tgt,J,tie_high);
  } else if (design==1||design==2) {
    BernRes br=bern_mtd(nv,sv,elim,phi_tgt,K,alpha_d,M_IS,ESS_thr,J,tie_high,rng);
    mtd_sel=br.mtd; ess=br.ess; fell_back=br.fell_back?1:0;
  } else { // crm
    std::vector<int> adm;
    for(int j=0;j<J;++j) if(nv[j]>0&&!elim[j]) adm.push_back(j);
    if(!adm.empty()){
      CrmRes cr2=crm_step(nv,sv,elim,skel,phi_tgt,crm_thr,crm_sd,crm_ng,J);
      int best=-1; double bd=1e9;
      for(int j:adm) if(!cr2.elim[j]){
        double d=std::abs(cr2.pi_mean[j]-phi_tgt);
        if(tie_high ? (d<=bd) : (d<bd)){bd=d;best=j;}
      }
      mtd_sel=best;
    }
  }

  // ---- operating characteristics ----
  // True MTD tie-break MUST match the end-of-trial selection convention above,
  // otherwise PCS is scored against a different dose than the estimators target.
  // tie_high=true reproduces the manuscript's higher-dose convention.
  int true_mtd=0; double bd=1e9;
  for(int j=0;j<J;++j){
    double d=std::abs(pi[j]-phi_tgt);
    if(tie_high ? (d<=bd) : (d<bd)){bd=d;true_mtd=j;}
  }

  int nt=std::accumulate(nv.begin(),nv.end(),0);
  int na=0;
  for(int j=true_mtd+1;j<J;++j) na+=nv[j];

  std::vector<int> elimv(J,0);
  for(int j=0;j<J;++j) elimv[j] = elim[j] ? 1 : 0;

  return {mtd_sel, true_mtd,
          (mtd_sel>=0&&mtd_sel==true_mtd)?1:0,
          (nt>0)?100.0*na/nt:0.0,
          std::accumulate(sv.begin(),sv.end(),0),
          (mtd_sel>=0)?(int)(mtd_sel>true_mtd):-1,
          (mtd_sel>=0)?(int)(mtd_sel<true_mtd):-1,
          (na>0)?1:0, na,
          ess, fell_back, nv, sv,
          elim[true_mtd] ? 1 : 0, elimv};
}

// ============================================================
// EXPORTED: build adaptive table
// [[Rcpp::export]]
List build_adap_tbl(double phi1, double phi2,
                    double a0,   double b0,
                    int N,       int m,
                    Nullable<NumericMatrix> loss = R_NilValue) {
  double L[9];
  if (loss.isNotNull()) {
    NumericMatrix Lm(loss.get());
    for (int a=0;a<3;++a) for (int h=0;h<3;++h) L[a*3+h]=Lm(a,h);
  } else {
    for (int i=0;i<9;++i) L[i]=LOSS_SYM[i];
  }
  List tbl; CharacterVector nms;
  for (int n=m; n<=N; n+=m) {
    IntegerVector acts(n+1);
    for (int s=0;s<=n;++s)
      acts[s]=bayes_action(s,n,phi1,phi2,a0,b0,L);
    tbl.push_back(acts);
    nms.push_back(std::to_string(n));
  }
  tbl.attr("names")=nms;
  return tbl;
}

// ============================================================
// EXPORTED: run n_sim trials, return aggregated results
// [[Rcpp::export]]
List run_scenario_cpp(NumericVector pi_r,
                  int design,
                  List atbl,
                  List p,          // params list
                  int n_sim,
                  int base_seed) {

  int J=pi_r.size();
  std::vector<double> pi(pi_r.begin(),pi_r.end());

  double phi_tgt=p["phi_tgt"],phi1=p["phi1"],phi2=p["phi2"],
         rho=p["rho"],a0=p["a0"],b0=p["b0"];
  int    N=p["N"], m=p["m"];
  double lam1=p["lam1"],lam2=p["lam2"];
  int    K=p["K"]; double alpha_d=p["alpha_d"],ESS_thr=p["ESS_thr"];
  int    M_IS=p["M_IS"];
  NumericVector skr=p["crm_skel"];
  std::vector<double> skel(skr.begin(),skr.end());
  double crm_sd=p["crm_sd"],crm_thr=p["crm_thr"];
  int    crm_ng=p["crm_ng"];

  // New (revision) parameters, with defaults that preserve behavior if absent.
  // tie_high defaults to TRUE (manuscript higher-dose convention).
  // Elimination defaults to STANDARD BOIN: overdose test at phi_tgt under Beta(1,1).
  bool   tie_high = p.containsElementNamed("tie_high") ? (bool)(int)p["tie_high"] : true;
  double phi_elim = p.containsElementNamed("phi_elim") ? (double)p["phi_elim"] : phi_tgt;
  double elim_a0  = p.containsElementNamed("elim_a0")  ? (double)p["elim_a0"]  : 1.0;
  double elim_b0  = p.containsElementNamed("elim_b0")  ? (double)p["elim_b0"]  : 1.0;

  // gBOINS calibration (Mu et al. 2021); defaults interpolate the paper's
  // phi0=0.2 and phi0=0.3 calibrations to phi_tgt=0.25.
  double gb_c1  = p.containsElementNamed("gb_c1")  ? (double)p["gb_c1"]  : std::log(1.075);
  double gb_c2  = p.containsElementNamed("gb_c2")  ? (double)p["gb_c2"]  : std::log(1.075)/3.0;
  double gb_eps = p.containsElementNamed("gb_eps") ? (double)p["gb_eps"] : 0.5;
  int    gb_N0  = p.containsElementNamed("gb_N0")  ? (int)p["gb_N0"]     : 6;

  // accumulators
  double spcs=0,spcs2=0,spod=0,sdlts=0;
  double sosel=0,susel=0,sany=0,snab=0,sfb=0,selim=0;
  std::vector<double> ess_v;
  std::vector<int> mtd_dist(J+1,0);

  for (int i=0;i<n_sim;++i){
    auto r=sim_trial(pi,design,atbl,
                     phi_tgt,phi1,phi2,rho,a0,b0,N,m,J,
                     lam1,lam2,K,alpha_d,M_IS,ESS_thr,
                     skel,crm_sd,crm_ng,crm_thr,
                     tie_high,phi_elim,elim_a0,elim_b0,
                     gb_c1,gb_c2,gb_eps,gb_N0,
                     (uint32_t)(base_seed+i));
    spcs+=r.pcs; spcs2+=(double)r.pcs*r.pcs; spod+=r.pod;
    sdlts+=r.dlts;
    if(r.over_sel>=0) sosel+=r.over_sel;
    if(r.under_sel>=0) susel+=r.under_sel;
    sany+=r.any_above; snab+=r.n_above; sfb+=r.fell_back;
    selim+=r.mtd_elim;
    if(r.ess>=0) ess_v.push_back(r.ess);
    int sel=(r.mtd_sel<0)?0:(r.mtd_sel+1);
    mtd_dist[sel]++;
  }

  double pcs=spcs/n_sim*100.0;
  double se=std::sqrt((spcs2/n_sim-(spcs/n_sim)*(spcs/n_sim))/n_sim)*100.0;

  double e_med=-1,e_p10=-1;
  if(!ess_v.empty()){
    std::sort(ess_v.begin(),ess_v.end());
    int ne=(int)ess_v.size();
    e_med=ess_v[ne/2]; e_p10=ess_v[(int)(ne*0.10)];
  }

  return List::create(
    Named("pcs")          =pcs,
    Named("pcs_se")       =se,
    Named("pod")          =spod/n_sim,
    Named("mean_dlts")    =sdlts/n_sim,
    Named("pct_over_sel") =sosel/n_sim*100.0,
    Named("pct_under_sel")=susel/n_sim*100.0,
    Named("pct_any_above")=sany/n_sim*100.0,
    Named("mean_n_above") =snab/n_sim,
    Named("pct_fallback") =sfb/n_sim*100.0,
    Named("pct_mtd_elim") =selim/n_sim*100.0,
    Named("ess_median")   =e_med,
    Named("ess_p10")      =e_p10,
    Named("mtd_dist")     =IntegerVector(mtd_dist.begin(),mtd_dist.end()),
    Named("n_sim")        =n_sim
  );
}

// ============================================================
// EXPORTED: prior mean curve (Section 4.3 validation)
// [[Rcpp::export]]
List prior_mean_curve(int K, int J) {
  std::vector<double> ds(J),Ewk(K+1),ppi(J),dlt(J);
  for(int j=0;j<J;++j) ds[j]=(double)j/(J-1);
  for(int k=0;k<=K;++k) Ewk[k]=(double)(k+1)/(K+2);
  for(int j=0;j<J;++j){
    double v=0.0;
    for(int k=0;k<=K;++k) v+=Ewk[k]*bern_basis(ds[j],k,K);
    ppi[j]=v; dlt[j]=std::abs(v-0.25);
  }
  return List::create(Named("prior_pi")=NumericVector(ppi.begin(),ppi.end()),
                      Named("delta")   =NumericVector(dlt.begin(),dlt.end()),
                      Named("E_wk")    =NumericVector(Ewk.begin(),Ewk.end()));
}

// ============================================================
// EXPORTED: single trial (testing/illustration)
// [[Rcpp::export]]
List one_trial(NumericVector pi_r, int design, List atbl, List p, int seed) {
  int J=pi_r.size();
  std::vector<double> pi(pi_r.begin(),pi_r.end());
  double phi_tgt=p["phi_tgt"],phi1=p["phi1"],phi2=p["phi2"],
         rho=p["rho"],a0=p["a0"],b0=p["b0"];
  int N=p["N"],m=p["m"];
  double lam1=p["lam1"],lam2=p["lam2"];
  int K=p["K"]; double alpha_d=p["alpha_d"],ESS_thr=p["ESS_thr"];
  int M_IS=p["M_IS"];
  NumericVector skr=p["crm_skel"];
  std::vector<double> skel(skr.begin(),skr.end());
  double crm_sd=p["crm_sd"],crm_thr=p["crm_thr"];
  int crm_ng=p["crm_ng"];

  bool   tie_high = p.containsElementNamed("tie_high") ? (bool)(int)p["tie_high"] : true;
  double phi_elim = p.containsElementNamed("phi_elim") ? (double)p["phi_elim"] : phi_tgt;
  double elim_a0  = p.containsElementNamed("elim_a0")  ? (double)p["elim_a0"]  : 1.0;
  double elim_b0  = p.containsElementNamed("elim_b0")  ? (double)p["elim_b0"]  : 1.0;
  double gb_c1  = p.containsElementNamed("gb_c1")  ? (double)p["gb_c1"]  : std::log(1.075);
  double gb_c2  = p.containsElementNamed("gb_c2")  ? (double)p["gb_c2"]  : std::log(1.075)/3.0;
  double gb_eps = p.containsElementNamed("gb_eps") ? (double)p["gb_eps"] : 0.5;
  int    gb_N0  = p.containsElementNamed("gb_N0")  ? (int)p["gb_N0"]     : 6;

  auto r=sim_trial(pi,design,atbl,phi_tgt,phi1,phi2,rho,a0,b0,N,m,J,
                   lam1,lam2,K,alpha_d,M_IS,ESS_thr,
                   skel,crm_sd,crm_ng,crm_thr,
                   tie_high,phi_elim,elim_a0,elim_b0,
                   gb_c1,gb_c2,gb_eps,gb_N0,(uint32_t)seed);

  return List::create(
    Named("mtd_sel") =r.mtd_sel+1,
    Named("true_mtd")=r.true_mtd+1,
    Named("pcs")     =r.pcs,
    Named("pod")     =r.pod,
    Named("dlts")    =r.dlts,
    Named("ess")     =r.ess,
    Named("fell_back")=r.fell_back,
    Named("n_vec")   =IntegerVector(r.nv.begin(),r.nv.end()),
    Named("s_vec")   =IntegerVector(r.sv.begin(),r.sv.end()),
    Named("mtd_elim")=r.mtd_elim,
    Named("elim_vec")=IntegerVector(r.elimv.begin(),r.elimv.end())
  );
}
