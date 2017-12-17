#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

const double err = 1e-10;

arma::mat projection_matrix(const arma::mat X, arma::mat projection,const int new_idx){
    int n = X.n_rows;
		arma::mat I(n,n,arma::fill::eye);
		arma::mat new_projection;
		new_projection = projection;
	  arma::mat xs=X.col(new_idx);
    new_projection = new_projection * (I - ((xs * xs.t() * new_projection) / arma::as_scalar(xs.t() * new_projection * xs))) ;

		return new_projection;
}

double logcomb(int n, int k){
    double ans=0;
    for(int i=0;i<k;i++){
        ans = ans + log(n-i) - log(i+1);
    }
    return ans;
}

double EBIC(const arma::colvec res_y, int s, int p, double r){
    double ebic;
    int n = res_y.size();
    ebic =  n*log(arma::as_scalar(res_y.t()*res_y)/n);
    ebic = ebic + s*log(n) + 2*(1- log(n)/(r*log(p)))*logcomb(p,s);
    return ebic;
}

arma::vec Normalization(const arma::vec x){
    int n = x.size();
    arma::vec nx(n);
    double m = arma::mean(x);
    double sd = arma::stddev(x);
    for (int i=0; i<n; i++){
        nx(i) = (x(i) - m) / sd;
    }
    return nx;
}

// [[Rcpp::export]]
List SLasso(const arma::colvec y, const arma::mat X, double r = 2){
    int n = X.n_rows;
    int p = X.n_cols;
    
    //normalization
    arma::vec norm_y;
    arma::mat norm_X(n,p);
    

    norm_y = Normalization(y);
    for (int i =0; i<p; i++){
        norm_X.col(i) = Normalization(X.col(i));
    }
    arma::uvec candidate = arma::linspace<arma::uvec>(0,p-1,p);
    arma::uvec selected;
    arma::uword temp_select;
    arma::mat proj(n,n,arma::fill::eye);
    double ebic;
    arma::colvec res_y=y;   
    int s = 0;
    ebic = arma::datum::inf;
    

    while(true){
        if (candidate.is_empty()){
            cout<<"All selected"<<endl;
            break;
        } else{
            arma::mat tmp = norm_X.cols(candidate);
            tmp = tmp.t() * res_y;
            tmp = arma::abs(tmp);
            temp_select = tmp.index_max();
            
            arma::uword tmp_s = candidate(temp_select);
            candidate.shed_row(temp_select);
            proj = projection_matrix(norm_X,proj,temp_select);
            res_y = proj * norm_y;
            double new_ebic = EBIC(res_y, s+1, p, r);
            cout<<"Iter:"<< s<< "  EBIC: " << new_ebic <<endl;//
            if(isnan( new_ebic )){
               break;
            }
            if ( (s>0) && (new_ebic>ebic) ){
                break;
            } else{
                arma::uvec a={tmp_s};
                selected = arma::join_cols(selected,a);
                ebic =new_ebic;
                s++;
            }
        }
    }
    arma::mat select_X;
    select_X = arma::join_rows(arma::ones<arma::vec>(n),X.cols(selected));
    arma::colvec coef;
    coef = arma::inv_sympd(select_X.t() * select_X) * select_X.t() * y;
    selected = selected + arma::ones<arma::uvec>(selected.size());
    return List::create(Named("Select.Features") = selected,
                        Named("Coef") = coef);
    
}
