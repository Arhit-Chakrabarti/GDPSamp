#include <RcppDist.h> 
#include <RcppArmadilloExtensions/sample.h> 
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace Rcpp;
using namespace RcppArmadillo;

// Function to calculate which Z sample = h for a fixed h
Rcpp::IntegerVector which_my(const arma::vec &x, const int &h){
  int n = x.size(); // Define the size of vector
  
  Rcpp::IntegerVector z(n); // To return the vector of labels
  for(int i = 0; i < n; i++){
    if(x(i) == h){  // If x_i = h, then return z_i = 1
      z(i) = 1;
    }else{          // If x_i != h, then return z_i = 0
      z(i) = 0;
    }
  }
  return z;         // Return the vector of labels
}
// Function to calculate n_h values. i.e. Rcpp equivalent of sum(data[which(z == h), ])
int count_my(const arma::vec &x, const int &h){
  int n = x.size();  // Define the size of vector
  
  Rcpp::IntegerVector z(n);  // To return the vector of labels
  for(int i = 0; i < n; i++){
    if(x(i) == h){           // If x_i = h, then return z_i = 1
      z(i) = 1;
    }else{
      z(i) = 0;              // If x_i != h, then return z_i = 0
    }
  }
  int count = sum(z);        // Calculate sum(z == h)
  return count;              // Return the count
}

// Function to subset the data. Rcpp equivalent of data[which(z == h), ]
arma::mat data_subset(const arma::mat &x,  const arma::vec &z, const int &h){
  int p = x.n_cols;   // Number of columns correspond to number of variables
  arma::uvec index = find(z == h); // Find which elements of z == h
  
  int n = index.size(); // Calculate the size of the index vector, which stores how many z == h
  arma::mat sub_mat(n, p, arma::fill::zeros); // Define a sub-matrix to store the subsetted matrix data[which(z == h), ]
  
  sub_mat = x.rows(index);  // Calculates data[which(z == h), ]
  return sub_mat; // Return the sub matrix
}

// Function to calculate cluster specific means. Rcpp equivalent of mean_my_vec(data[ , z == h, drop = FALSE]) 
// [[Rcpp::export]]
arma::rowvec mean_my_vec_c(const arma::mat& x) {
  int n = x.n_rows; // Defining the number of observations
  int p = x.n_cols; // Defining the number of variables
  // Initiate a column vector y, filled with 0 s
  arma::rowvec y(p, arma::fill::zeros);
  if(n > 1){
    y = mean(x, 0); // Here x is a data matrix of dimension nxp. Take column-means if number of observations >1
  }
  return(y);        // Return the mean vector in the form of a row-vector
}

// Function to calculate cluster specific SS matrix. Rcpp equivalent of ss_my(data[ , z == h, drop = FALSE])
// [[Rcpp::export]]
arma::mat ss_my_c(const arma::mat& x) {
  int n = x.n_rows; // Defining the number of observations
  int p = x.n_cols; // Defining the number of variables
  // Initiate a column vector y, filled with 0 s
  arma::mat y(p, p, arma::fill::zeros);
  if(n > 1){
    y = (n - 1)*cov(x); // Here x is a data matrix of dimension nxp. Take variance of each column if n >1 and multiply to get SS
  }
  return(y);            // Return the SS matrix
}


//[[Rcpp::export]]
Rcpp::List sample_mu_tau_c(const arma::mat &data,  const arma::vec &z, const int &L, const int&nu, const arma::mat &W, const arma::colvec &prior_mean, const double&prior_prec){
  int p = data.n_cols;  // Defining the number of variables
  Rcpp::IntegerVector n_h(L); // Defining the vector n_h;
  
  arma::mat x_h(L, p, arma::fill::zeros); // Define a matrix to store cluster specific data means. It is a Lxp matrix
  arma::cube ss_h(p, p, L, arma::fill::zeros); // Define an array to store cluster specific data SS matrix. It is a pxpxL array
  // To store Wishart distribution parameters
  arma::cube W_h(p, p, L);  // Define an array to store cluster specific data scale matrix for posterior Wishart distribution. It is a pxpxL array
  IntegerVector nu_h(L);    // Define an integer vector to store cluster specific degrees of freedom for posterior Wishart distribution.
  
  // To store MVN distribution parameters
  arma::mat mu_hat(p, L); // Define a matrix to store cluster specific mean matrix for posterior MVN distribution.
  arma::mat diffr(p, L);  // Define a matrix to store cluster specific (x_h - prior_mean) vector
  
  for(int h = 0; h < L; h++){
    n_h(h) = count_my(z, (h + 1)); // Calculate the cluster h sample sizes
    arma::mat temp = data_subset(data, z, (h + 1)); // Calculate the subsetted data matrix for cluster h
    x_h.row(h) = mean_my_vec_c(temp); // Calculate the mean vector for cluster h with subsetted data
    ss_h.slice(h) = ss_my_c(temp); // Calculate the SS matrix for cluster h with subsetted data
  }
  
  NumericVector n_h_double = as<NumericVector>(n_h); // As n_h is an Integer vector, convert it to numeric vector
  Rcpp::NumericVector prior_prec_hat = (prior_prec + n_h_double); // Updated precision of posterior MVN distribution
  
  for(int h = 0; h < L; h++){
    nu_h(h) = nu + n_h(h); // Calculate the updated cluster specific degrees of freedom for posterior Wishart distribution
    
    // nu_h[h] = nu + n_h(h); // Calculate the updated cluster specific degrees of freedom for posterior Wishart distribution
    diffr.col(h) = x_h.row(h).t() - prior_mean;  // Calculate cluster specific (x_h - prior_mean) vector
    // Calculate the updated scale matrix for posterior Wishart distribution
    W_h.slice(h) = arma::pinv(arma::pinv(W, 0.000000001) + ss_h.slice(h) +  ((n_h(h) * prior_prec)/(n_h(h) + prior_prec)) * (diffr.col(h) * diffr.col(h).t()), 0.0000001) ;
    mu_hat.col(h) = (prior_prec*prior_mean + n_h(h)*x_h.row(h).t())/(prior_prec_hat(h)); // Updated mean for the posterior Normal distribution
  }
  
  arma::cube tau_draw(p, p, L); // Define an array to store draws from posterior Wishart distribution
  arma::mat mu_draw(p, L);      // Define a matrix to store draws from posterior MVN distribution
  arma::cube var2(p, p, L);
  arma::cube var(p, p, L);
  for(int h = 0; h < L; h++){
    // Draw Precision matrix from a Wishart distribution
    tau_draw.slice(h) = arma::wishrnd(arma::symmatu(W_h.slice(h)), nu_h(h));
    // Draw mean vector from a Multivariate Normal distribution
    var2.slice(h) = arma::pinv(tau_draw.slice(h), 0.00000000000000001);
    var.slice(h) = arma::pinv(tau_draw.slice(h), 0.00000000000000001)/prior_prec_hat(h);
    mu_draw.col(h) = mvnrnd(mu_hat.col(h), arma::symmatu(var.slice(h)));
  }
  return Rcpp::List::create(Rcpp::Named("var.hat") = var, Rcpp::Named("var2.hat") = var2, Rcpp::Named("prior_prec") = prior_prec_hat, Rcpp::Named("tau.draw") = tau_draw, Rcpp::Named("mu.draw") = mu_draw,  Rcpp::Named("n_h") = n_h);
  
}

// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov){
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);    
}

// [[Rcpp::export]]
arma::vec dmvnorm_my( arma::mat x,  arma::mat mean,  arma::mat sigma, bool log){ 
  
  arma::vec distval = Mahalanobis(x,  mean, sigma);
  
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  double log2pi = 1.8378770664093454835606594728112352797227949472755668;
  arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
  
  if(log){ 
    return(logretval);
    
  }else { 
    return(exp(logretval));
  }
}
// Function to calculate the density for multivariate normal distribution
// [[Rcpp::export]]
arma::mat multi_normal(const arma::mat& data, const int& L, const arma::mat& mu, const arma::cube& tau, const bool& log_prob) {
  int n = data.n_rows; // Defining the number of observations
  // Initiate a matrix of dimension n x L, filled with 0 s
  arma::mat density(n, L, arma::fill::zeros);
  for(int h = 0; h < L; h++){
    density.col(h) = dmvnorm_my(data, mu.row(h), arma::symmatu(arma::pinv(tau.slice(h), 0.00000000000000001)), log_prob); // Calculate the density for multivariate normal distribution
  }
  return(density);
}


// [[Rcpp::export]]
Rcpp::IntegerVector sample_my(const arma::mat& prob) {
  int n = prob.n_rows; // Defining the number of observations
  Rcpp::IntegerVector z(n);  // Defining the vector z
  int L = prob.n_cols; // Define the number of columns which corresponds to the number of clusters
  arma::vec cluster = arma::regspace(1,  L); //Define a vector to store labels 1, 2, ... , L
  
  for(int i = 0; i < n; i++){
    z(i) = Rcpp::RcppArmadillo::sample(cluster, 1, true, prob.row(i).t())[0];
  }
  return(z);
}


// [[Rcpp::export]]
Rcpp::IntegerVector sample_z_c(const arma::mat& data, const arma::vec& beta, const int& L, const arma::mat& mu, const arma::cube& tau) {
  int n = data.n_rows; // Defining the number of observations
  arma::mat probability(n, L, arma::fill::zeros); // Initiate a probability matrix of dimension n x L, filled with 0 s
  Rcpp::IntegerVector z(n); // Defining the vector z, to return a vector of cluster labels
  arma::mat normal_density(n, L, arma::fill::zeros); // Matrix to store MVN density at the data points
  
  normal_density = multi_normal(data, L, mu, tau, true); // Calculate the log of the MVN distribution density
  
  for(int i = 0; i < n; i++){
    for(int h = 0; h < L; h++){
      probability(i, h) = log(beta(h)) + normal_density(i, h) ; // Calculate the probability in log-scale
    }
    probability.row(i) = exp(probability.row(i) - max(probability.row(i)));  // Log-sum-exp trick
    probability.row(i) = probability.row(i) + 0.0000000001; // Add a very small number to stabilize computation
    probability.row(i) = probability.row(i)/sum(probability.row(i)); // Normalize to make it a probability
  }
  z = sample_my(probability); // Sample the labels with the probability matrix generated
  
  return(z);
}

// [[Rcpp::export]]
arma::mat prob(const arma::mat& data, const arma::vec& beta, const int& L, const arma::mat& mu, const arma::cube& tau) {
  int n = data.n_rows; // Defining the number of observations
  arma::mat probability(n, L, arma::fill::zeros); // Initiate a probability matrix of dimension n x L, filled with 0 s
  Rcpp::IntegerVector z(n); // Defining the vector z, to return a vector of cluster labels
  arma::mat normal_density(n, L, arma::fill::zeros); // Matrix to store MVN density at the data points
  
  normal_density = multi_normal(data, L, mu, tau, true); // Calculate the log of the MVN distribution density
  
  for(int i = 0; i < n; i++){
    for(int h = 0; h < L; h++){
      probability(i, h) = log(beta(h)) + normal_density(i, h) ; // Calculate the probability in log-scale
    }
    probability.row(i) = exp(probability.row(i) - max(probability.row(i)));  // Log-sum-exp trick
    probability.row(i) = probability.row(i) + 0.0000000001; // Add a very small number to stabilize computation
    probability.row(i) = probability.row(i)/sum(probability.row(i)); // Normalize to make it a probability
  }
  z = sample_my(probability); // Sample the labels with the probability matrix generated
  
  return(probability);
}

//[[Rcpp::export]]
Rcpp::List sample_mu_tau_gdp_c(const arma::mat &data1, const arma::mat &data2, const arma::mat &data3, const arma::mat &data4, const arma::mat &data5, const arma::mat &data6, const arma::mat &data7, const arma::mat &data8,        
                               const arma::vec &z1, const arma::vec &z2, const arma::vec &z3, const arma::vec &z4, const arma::vec &z5, const arma::vec &z6, const arma::vec &z7, const arma::vec &z8,     
                               const int &L, const int&nu, const arma::mat &W, const arma::colvec &prior_mean, const double&prior_prec){
  int p = data1.n_cols;  // Defining the number of variables
  Rcpp::IntegerVector n1_h(L); // Defining the vector n1_h;
  Rcpp::IntegerVector n2_h(L); // Defining the vector n2_h;
  Rcpp::IntegerVector n3_h(L); // Defining the vector n3_h;
  Rcpp::IntegerVector n4_h(L); // Defining the vector n4_h;
  Rcpp::IntegerVector n5_h(L); // Defining the vector n5_h;
  Rcpp::IntegerVector n6_h(L); // Defining the vector n6_h;
  Rcpp::IntegerVector n7_h(L); // Defining the vector n7_h;
  Rcpp::IntegerVector n8_h(L); // Defining the vector n8_h;
  
  
  arma::mat x1_h(L, p, arma::fill::zeros); // Define a matrix to store cluster specific data means. It is a Lxp matrix
  arma::mat x2_h(L, p, arma::fill::zeros); // Define a matrix to store cluster specific data means. It is a Lxp matrix
  arma::mat x3_h(L, p, arma::fill::zeros); // Define a matrix to store cluster specific data means. It is a Lxp matrix
  arma::mat x4_h(L, p, arma::fill::zeros); // Define a matrix to store cluster specific data means. It is a Lxp matrix
  arma::mat x5_h(L, p, arma::fill::zeros); // Define a matrix to store cluster specific data means. It is a Lxp matrix
  arma::mat x6_h(L, p, arma::fill::zeros); // Define a matrix to store cluster specific data means. It is a Lxp matrix
  arma::mat x7_h(L, p, arma::fill::zeros); // Define a matrix to store cluster specific data means. It is a Lxp matrix
  arma::mat x8_h(L, p, arma::fill::zeros); // Define a matrix to store cluster specific data means. It is a Lxp matrix
  
  arma::cube ss1_h(p, p, L, arma::fill::zeros); // Define an array to store cluster specific data SS matrix. It is a pxpxL array
  arma::cube ss2_h(p, p, L, arma::fill::zeros); // Define an array to store cluster specific data SS matrix. It is a pxpxL array
  arma::cube ss3_h(p, p, L, arma::fill::zeros); // Define an array to store cluster specific data SS matrix. It is a pxpxL array
  arma::cube ss4_h(p, p, L, arma::fill::zeros); // Define an array to store cluster specific data SS matrix. It is a pxpxL array
  arma::cube ss5_h(p, p, L, arma::fill::zeros); // Define an array to store cluster specific data SS matrix. It is a pxpxL array
  arma::cube ss6_h(p, p, L, arma::fill::zeros); // Define an array to store cluster specific data SS matrix. It is a pxpxL array
  arma::cube ss7_h(p, p, L, arma::fill::zeros); // Define an array to store cluster specific data SS matrix. It is a pxpxL array
  arma::cube ss8_h(p, p, L, arma::fill::zeros); // Define an array to store cluster specific data SS matrix. It is a pxpxL array
  
  // To store Wishart distribution parameters
  arma::cube W1_h(p, p, L);  // Define an array to store cluster specific data scale matrix (first part) for posterior Wishart distribution. It is a pxpxL array
  arma::cube W2_h(p, p, L);  // Define an array to store cluster specific data scale matrix (second part) for posterior Wishart distribution. It is a pxpxL array
  arma::cube W3_h(p, p, L);  // Define an array to store cluster specific data scale matrix (third part) for posterior Wishart distribution. It is a pxpxL array
  IntegerVector nu_hat(L);    // Define an integer vector to store cluster specific degrees of freedom for posterior Wishart distribution.
  
  // To store MVN distribution parameters
  arma::mat mu_hat(p, L); // Define a matrix to store cluster specific mean matrix for posterior MVN distribution.
  arma::mat diffr1(p, L);  // Define a matrix to store cluster specific (x1_h - prior_mean) vector
  arma::mat diffr2(p, L);  // Define a matrix to store cluster specific (x2_h - prior_mean) vector
  arma::mat diffr3(p, L);  // Define a matrix to store cluster specific (x3_h - prior_mean) vector
  arma::mat diffr4(p, L);  // Define a matrix to store cluster specific (x4_h - prior_mean) vector
  arma::mat diffr5(p, L);  // Define a matrix to store cluster specific (x5_h - prior_mean) vector
  arma::mat diffr6(p, L);  // Define a matrix to store cluster specific (x6_h - prior_mean) vector
  arma::mat diffr7(p, L);  // Define a matrix to store cluster specific (x7_h - prior_mean) vector
  arma::mat diffr8(p, L);  // Define a matrix to store cluster specific (x8_h - prior_mean) vector
  
  for(int h = 0; h < L; h++){
    //Group 1
    n1_h(h) = count_my(z1, (h + 1)); // Calculate the cluster h sample sizes
    arma::mat temp1 = data_subset(data1, z1, (h + 1)); // Calculate the subsetted data matrix for cluster h
    x1_h.row(h) = mean_my_vec_c(temp1); // Calculate the mean vector for cluster h with subsetted data
    ss1_h.slice(h) = ss_my_c(temp1); // Calculate the SS matrix for cluster h with subsetted data
    //Group 2
    n2_h(h) = count_my(z2, (h + 1)); // Calculate the cluster h sample sizes
    arma::mat temp2 = data_subset(data2, z2, (h + 1)); // Calculate the subsetted data matrix for cluster h
    x2_h.row(h) = mean_my_vec_c(temp2); // Calculate the mean vector for cluster h with subsetted data
    ss2_h.slice(h) = ss_my_c(temp2); // Calculate the SS matrix for cluster h with subsetted data
    //Group 3
    n3_h(h) = count_my(z3, (h + 1)); // Calculate the cluster h sample sizes
    arma::mat temp3 = data_subset(data3, z3, (h + 1)); // Calculate the subsetted data matrix for cluster h
    x3_h.row(h) = mean_my_vec_c(temp3); // Calculate the mean vector for cluster h with subsetted data
    ss3_h.slice(h) = ss_my_c(temp3); // Calculate the SS matrix for cluster h with subsetted data
    //Group 4
    n4_h(h) = count_my(z4, (h + 1)); // Calculate the cluster h sample sizes
    arma::mat temp4 = data_subset(data4, z4, (h + 1)); // Calculate the subsetted data matrix for cluster h
    x4_h.row(h) = mean_my_vec_c(temp4); // Calculate the mean vector for cluster h with subsetted data
    ss4_h.slice(h) = ss_my_c(temp4); // Calculate the SS matrix for cluster h with subsetted data
    //Group 5
    n5_h(h) = count_my(z5, (h + 1)); // Calculate the cluster h sample sizes
    arma::mat temp5 = data_subset(data5, z5, (h + 1)); // Calculate the subsetted data matrix for cluster h
    x5_h.row(h) = mean_my_vec_c(temp5); // Calculate the mean vector for cluster h with subsetted data
    ss5_h.slice(h) = ss_my_c(temp5); // Calculate the SS matrix for cluster h with subsetted data
    //Group 6
    n6_h(h) = count_my(z6, (h + 1)); // Calculate the cluster h sample sizes
    arma::mat temp6 = data_subset(data6, z6, (h + 1)); // Calculate the subsetted data matrix for cluster h
    x6_h.row(h) = mean_my_vec_c(temp6); // Calculate the mean vector for cluster h with subsetted data
    ss6_h.slice(h) = ss_my_c(temp6); // Calculate the SS matrix for cluster h with subsetted data
    //Group 7
    n7_h(h) = count_my(z7, (h + 1)); // Calculate the cluster h sample sizes
    arma::mat temp7 = data_subset(data7, z7, (h + 1)); // Calculate the subsetted data matrix for cluster h
    x7_h.row(h) = mean_my_vec_c(temp7); // Calculate the mean vector for cluster h with subsetted data
    ss7_h.slice(h) = ss_my_c(temp7); // Calculate the SS matrix for cluster h with subsetted data
    //Group 8
    n8_h(h) = count_my(z8, (h + 1)); // Calculate the cluster h sample sizes
    arma::mat temp8 = data_subset(data8, z8, (h + 1)); // Calculate the subsetted data matrix for cluster h
    x8_h.row(h) = mean_my_vec_c(temp8); // Calculate the mean vector for cluster h with subsetted data
    ss8_h.slice(h) = ss_my_c(temp8); // Calculate the SS matrix for cluster h with subsetted data
  }
  
  NumericVector n1_h_double = as<NumericVector>(n1_h); // As n_h is an Integer vector, convert it to numeric vector
  NumericVector n2_h_double = as<NumericVector>(n2_h); // As n_h is an Integer vector, convert it to numeric vector
  NumericVector n3_h_double = as<NumericVector>(n3_h); // As n_h is an Integer vector, convert it to numeric vector
  NumericVector n4_h_double = as<NumericVector>(n4_h); // As n_h is an Integer vector, convert it to numeric vector
  NumericVector n5_h_double = as<NumericVector>(n5_h); // As n_h is an Integer vector, convert it to numeric vector
  NumericVector n6_h_double = as<NumericVector>(n6_h); // As n_h is an Integer vector, convert it to numeric vector
  NumericVector n7_h_double = as<NumericVector>(n7_h); // As n_h is an Integer vector, convert it to numeric vector
  NumericVector n8_h_double = as<NumericVector>(n8_h); // As n_h is an Integer vector, convert it to numeric vector
  
  Rcpp::NumericVector prior_prec_hat = (prior_prec + n1_h_double + n2_h_double + n3_h_double + n4_h_double + n5_h_double + n6_h_double + n7_h_double + n8_h_double); // Updated precision of posterior MVN distribution
  
  for(int h = 0; h < L; h++){
    nu_hat(h) = nu + n1_h(h) + n2_h(h) + n3_h(h) + n4_h(h) + n5_h(h) + n6_h(h) + n7_h(h) + n8_h(h); // Calculate the updated cluster specific degrees of freedom for posterior Wishart distribution
    //Group specific deviations
    diffr1.col(h) = x1_h.row(h).t() - prior_mean;  // Calculate cluster specific (x1_h - prior_mean) vector
    diffr2.col(h) = x2_h.row(h).t() - prior_mean;  // Calculate cluster specific (x2_h - prior_mean) vector
    diffr3.col(h) = x3_h.row(h).t() - prior_mean;  // Calculate cluster specific (x3_h - prior_mean) vector
    diffr4.col(h) = x4_h.row(h).t() - prior_mean;  // Calculate cluster specific (x4_h - prior_mean) vector
    diffr5.col(h) = x5_h.row(h).t() - prior_mean;  // Calculate cluster specific (x5_h - prior_mean) vector
    diffr6.col(h) = x6_h.row(h).t() - prior_mean;  // Calculate cluster specific (x6_h - prior_mean) vector
    diffr7.col(h) = x7_h.row(h).t() - prior_mean;  // Calculate cluster specific (x7_h - prior_mean) vector
    diffr8.col(h) = x8_h.row(h).t() - prior_mean;  // Calculate cluster specific (x8_h - prior_mean) vector
    // Calculate the updated scale matrix for posterior Wishart distribution
    W1_h.slice(h) = ss1_h.slice(h) + ss2_h.slice(h) + ss3_h.slice(h) + ss4_h.slice(h) + ss5_h.slice(h) + ss6_h.slice(h) + ss7_h.slice(h) + ss8_h.slice(h) ; 
    
    W2_h.slice(h) = ((n1_h(h) * prior_prec)/(n1_h(h) + prior_prec)) * (diffr1.col(h) * diffr1.col(h).t()) + 
      ((n2_h(h) * prior_prec)/(n2_h(h) + prior_prec)) * (diffr2.col(h) * diffr2.col(h).t()) + 
      ((n3_h(h) * prior_prec)/(n3_h(h) + prior_prec)) * (diffr3.col(h) * diffr3.col(h).t()) + 
      ((n4_h(h) * prior_prec)/(n4_h(h) + prior_prec)) * (diffr4.col(h) * diffr4.col(h).t()) + 
      ((n5_h(h) * prior_prec)/(n5_h(h) + prior_prec)) * (diffr5.col(h) * diffr5.col(h).t()) + 
      ((n6_h(h) * prior_prec)/(n6_h(h) + prior_prec)) * (diffr6.col(h) * diffr6.col(h).t()) + 
      ((n7_h(h) * prior_prec)/(n7_h(h) + prior_prec)) * (diffr7.col(h) * diffr7.col(h).t()) + 
      ((n8_h(h) * prior_prec)/(n8_h(h) + prior_prec)) * (diffr8.col(h) * diffr8.col(h).t());   
    
    W3_h.slice(h) = n1_h(h)*n2_h(h)*((x1_h.row(h).t() - x2_h.row(h).t()) * (x1_h.row(h).t() - x2_h.row(h).t()).t()) +
      n1_h(h)*n3_h(h)*((x1_h.row(h).t() - x3_h.row(h).t()) * (x1_h.row(h).t() - x3_h.row(h).t()).t()) +
      n1_h(h)*n4_h(h)*((x1_h.row(h).t() - x4_h.row(h).t()) * (x1_h.row(h).t() - x4_h.row(h).t()).t()) +
      n1_h(h)*n5_h(h)*((x1_h.row(h).t() - x5_h.row(h).t()) * (x1_h.row(h).t() - x5_h.row(h).t()).t()) + 
      n1_h(h)*n6_h(h)*((x1_h.row(h).t() - x6_h.row(h).t()) * (x1_h.row(h).t() - x6_h.row(h).t()).t()) +
      n1_h(h)*n7_h(h)*((x1_h.row(h).t() - x7_h.row(h).t()) * (x1_h.row(h).t() - x7_h.row(h).t()).t()) + 
      n1_h(h)*n8_h(h)*((x1_h.row(h).t() - x8_h.row(h).t()) * (x1_h.row(h).t() - x8_h.row(h).t()).t()) +
      n2_h(h)*n3_h(h)*((x2_h.row(h).t() - x3_h.row(h).t()) * (x2_h.row(h).t() - x3_h.row(h).t()).t()) +
      n2_h(h)*n4_h(h)*((x2_h.row(h).t() - x4_h.row(h).t()) * (x2_h.row(h).t() - x4_h.row(h).t()).t()) +
      n2_h(h)*n5_h(h)*((x2_h.row(h).t() - x5_h.row(h).t()) * (x2_h.row(h).t() - x5_h.row(h).t()).t()) + 
      n2_h(h)*n6_h(h)*((x2_h.row(h).t() - x6_h.row(h).t()) * (x2_h.row(h).t() - x6_h.row(h).t()).t()) +
      n2_h(h)*n7_h(h)*((x2_h.row(h).t() - x7_h.row(h).t()) * (x2_h.row(h).t() - x7_h.row(h).t()).t()) + 
      n2_h(h)*n8_h(h)*((x2_h.row(h).t() - x8_h.row(h).t()) * (x2_h.row(h).t() - x8_h.row(h).t()).t()) +
      n3_h(h)*n4_h(h)*((x3_h.row(h).t() - x4_h.row(h).t()) * (x3_h.row(h).t() - x4_h.row(h).t()).t()) +
      n3_h(h)*n5_h(h)*((x3_h.row(h).t() - x5_h.row(h).t()) * (x3_h.row(h).t() - x5_h.row(h).t()).t()) + 
      n3_h(h)*n6_h(h)*((x3_h.row(h).t() - x6_h.row(h).t()) * (x3_h.row(h).t() - x6_h.row(h).t()).t()) +
      n3_h(h)*n7_h(h)*((x3_h.row(h).t() - x7_h.row(h).t()) * (x3_h.row(h).t() - x7_h.row(h).t()).t()) + 
      n3_h(h)*n8_h(h)*((x3_h.row(h).t() - x8_h.row(h).t()) * (x3_h.row(h).t() - x8_h.row(h).t()).t()) +
      n4_h(h)*n5_h(h)*((x4_h.row(h).t() - x5_h.row(h).t()) * (x4_h.row(h).t() - x5_h.row(h).t()).t()) + 
      n4_h(h)*n6_h(h)*((x4_h.row(h).t() - x6_h.row(h).t()) * (x4_h.row(h).t() - x6_h.row(h).t()).t()) +
      n4_h(h)*n7_h(h)*((x4_h.row(h).t() - x7_h.row(h).t()) * (x4_h.row(h).t() - x7_h.row(h).t()).t()) + 
      n4_h(h)*n8_h(h)*((x4_h.row(h).t() - x8_h.row(h).t()) * (x4_h.row(h).t() - x8_h.row(h).t()).t()) +
      n5_h(h)*n6_h(h)*((x5_h.row(h).t() - x6_h.row(h).t()) * (x5_h.row(h).t() - x6_h.row(h).t()).t()) +
      n5_h(h)*n7_h(h)*((x5_h.row(h).t() - x7_h.row(h).t()) * (x5_h.row(h).t() - x7_h.row(h).t()).t()) + 
      n5_h(h)*n8_h(h)*((x5_h.row(h).t() - x8_h.row(h).t()) * (x5_h.row(h).t() - x8_h.row(h).t()).t()) +
      n6_h(h)*n7_h(h)*((x6_h.row(h).t() - x7_h.row(h).t()) * (x6_h.row(h).t() - x7_h.row(h).t()).t()) + 
      n6_h(h)*n8_h(h)*((x6_h.row(h).t() - x8_h.row(h).t()) * (x6_h.row(h).t() - x8_h.row(h).t()).t()) +
      n7_h(h)*n8_h(h)*((x7_h.row(h).t() - x8_h.row(h).t()) * (x7_h.row(h).t() - x8_h.row(h).t()).t());
    
    mu_hat.col(h) = (prior_prec*prior_mean + n1_h(h)*x1_h.row(h).t() + n2_h(h)*x2_h.row(h).t() + n3_h(h)*x3_h.row(h).t() + 
      n4_h(h)*x4_h.row(h).t() + n5_h(h)*x5_h.row(h).t() + n6_h(h)*x6_h.row(h).t() + n7_h(h)*x7_h.row(h).t() + n8_h(h)*x8_h.row(h).t())/(prior_prec_hat(h)); // Updated mean for the posterior Normal distribution
  }
  
  arma::cube W_h(p, p, L);
  
  for(int h = 0; h < L; h++){
    W_h.slice(h) = arma::pinv(arma::pinv(W, 0.000000001) + W1_h.slice(h) + (W2_h.slice(h)/prior_prec_hat(h)) + (W3_h.slice(h)/prior_prec_hat(h)), 0.000000001);
  }
  arma::cube tau_draw(p, p, L); // Define an array to store draws from posterior Wishart distribution
  arma::mat mu_draw(p, L);      // Define a matrix to store draws from posterior MVN distribution
  arma::cube var(p, p, L);
  for(int h = 0; h < L; h++){
    // Draw Precision matrix from a Wishart distribution
    tau_draw.slice(h) = arma::wishrnd(W_h.slice(h), nu_hat(h));
    // Draw mean vector from a Multivariate Normal distribution
    var.slice(h) = arma::pinv(tau_draw.slice(h), 0.000000001)/prior_prec_hat(h);
    mu_draw.col(h) = mvnrnd(mu_hat.col(h), arma::symmatu(var.slice(h)));
  }
  return Rcpp::List::create(Rcpp::Named("tau.draw") = tau_draw, Rcpp::Named("mu.draw") = mu_draw, 
                            Rcpp::Named("n1_h") = n1_h,
                            Rcpp::Named("n2_h") = n2_h,
                            Rcpp::Named("n3_h") = n3_h,
                            Rcpp::Named("n4_h") = n4_h,
                            Rcpp::Named("n5_h") = n5_h,
                            Rcpp::Named("n6_h") = n6_h,
                            Rcpp::Named("n7_h") = n7_h,
                            Rcpp::Named("n8_h") = n8_h);
  
}

// [[Rcpp::export]]
double log_l(const arma::mat& data, const arma::vec& beta, const Rcpp::IntegerVector &z, const arma::mat& mu, const arma::cube& tau){
  int n = data.n_rows; //Calculate the number of observations
  int L = beta.size(); //Calculate the number of clusters
  arma::colvec Beta = beta;
  arma::mat log_like_beta(1, n);
  arma::mat log_like_normal(1, n);
  int index;
  arma::mat log_like_full(1, n);
  
  for(int i = 0; i < n; i++){
    index = z(i) - 1;
    log_like_beta.col(i) = log(Beta.row(index));
    log_like_normal.col(i) = dmvnorm(data.row(i), mu.row(index).t(), arma::symmatu(arma::pinv(tau.slice(index), 0.000000000001)), true);
    log_like_full.col(i) = log_like_beta.col(i) + log_like_normal.col(i);
  }
  return arma::accu(log_like_full);
}
